function Turner_generateImages(newIDs)

    %%% Created jpegs (montages/slices) for display %%%

    %%% Set the Matlab environment variable path
    setMatlabEnv;
    
    niiDir = '/Volumes/Projects/TURNER/Data/nifti/';
    qaDir = '/Users/harshads/Sites/Turner/QA/';

    outPrefix = 'f';

    conn = mysqlConnect;
    
    %%% If no subjectIDs specified, select all
    if ~exist('newIDs', 'var')
        results = mysqlQuery(conn, 'select distinct subjectID, sessionID, examID from imaging order by examID desc');
        subjectIDs = results(:, 1);
        sessionIDs = results(:, 2);
        examIDs = results(:, 3);
    %%% else get sessionIDs, examIDs for specified subjects
    else
        subjectIDs = {};
        sessionIDs = {};
        examIDs = {};
        for n = 1:length(newIDs)
            results = mysqlQuery(conn, sprintf('select distinct subjectID, sessionID, examID from imaging where subjectID = ''%s'' order by examID desc',  newIDs{n}));
            subjectIDs = [subjectIDs; results(:, 1)];
            sessionIDs = [sessionIDs; results(:, 2)];
            examIDs = [examIDs; results(:, 3)];
        end
        
    end
    
    for s = 1:length(subjectIDs)

        imgDirs = mysqlQuery(conn, sprintf('select imagecentralLocation from imaging where subjectID=%d and sessionID=''%s'' and examID=%d', subjectIDs{s}, sessionIDs{s}, examIDs{s}));
        
        for i = 1:length(imgDirs)
            
            fprintf('%s : ', imgDirs{i});
            
            dataPaths = strsplit(imgDirs{i}, '/');
            
            %%% For anatomical images
            if strcmp(dataPaths{3}, 'Structural')
                if any(cellfun(@(x) strcmpi(dataPaths{5}(5:end), x), {'localizer', 'asset'}))
                    fprintf('Skipped\n');
                    continue;
                end
                
                niiPath = fullfile(niiDir, dataPaths{1}, dataPaths{2}, dataPaths{4}, dataPaths{5});
                qaPath = fullfile(qaDir, dataPaths{1}, dataPaths{2}, dataPaths{4}(2:end), dataPaths{5});
                niiFile = fullfile(niiPath, [outPrefix '.nii.gz']);
                
                if exist(niiFile, 'file') 
                    createImages(niiPath, qaPath);
                else
                    fprintf('%s does not exist\n', niiFile);
                end
                
            elseif ismember(dataPaths{3}, {'DTI', 'Functional'})

                examID = num2str(cell2mat(mysqlQuery(conn, ['select examID from imaging where imageCentralLocation=''' imgDirs{i} ''''])));
                
                niiPath = fullfile(niiDir, dataPaths{1}, dataPaths{2}, ['E' examID], strrep(dataPaths{4}, ['_E' examID], ''));
                niiFile = fullfile(niiPath, [outPrefix '.nii.gz']);
                
                mcPath = fullfile(niiPath, 'FSL_MotionCorrect', [outPrefix '_mc']);
                mcFile = [mcPath '.nii.gz'];
                
                qaPath = fullfile(qaDir, dataPaths{1}, dataPaths{2}, examID, strrep(dataPaths{4}, ['_E' examID], ''));
                if exist(niiFile, 'file') && exist(mcFile, 'file')
                    createMovie(niiFile, mcFile, qaPath, dataPaths{3});
                    plotMotion(mcPath, qaPath, imgDirs{i}, conn);
                else
                    fprintf('%s does not exist\n', niiFile);
                end
                
                mrVista_motionFile = fullfile(niiPath, 'mrVista_motion.png');
                if exist(mrVista_motionFile, 'file') && ~exist(fullfile(qaPath, 'mrVista_motion.png'), 'file')
                   system(['cp ' mrVista_motionFile ' ' qaPath]); 
                end
            end 
        end
        
    end
    
    if matlabpool('size')
        matlabpool CLOSE FORCE
    end


    
    
function createImages(imgDir, outPath)

    jpgdir = fullfile(outPath, 'slices');
    
    imgfile = fullfile(imgDir, 'f.nii.gz');
    
    info = niftiRead(imgfile);
    img = info.data;
    if any(diff(info.pixdim))
        maxdim = find(info.pixdim == max(info.pixdim));
    else
        maxdim = 1;
    end
    
    if maxdim == 1
        imgMean = mean(mean(img, 2), 3);
        img = img(find(imgMean > 110), :, :);
    elseif maxdim == 2
        imgMean = mean(mean(img, 1), 3);
        img = img(:, find(imgMean > 110), :);
    elseif maxdim == 3
        imgMean = mean(mean(img, 1), 2);
        img = img(:, :, find(imgMean > 110));
    end
        
    %%% exclude slices without brain voxels
    imgscale=3.5;
    
    nImages = size(img, maxdim);
    
    jpgfile = arrayfun(@(i) fullfile(jpgdir, sprintf('slice%03d.jpg', i)), 1:nImages, 'UniformOutput', false);
    if ~any(~cellfun(@exist, jpgfile))
        fprintf('Already created\n');
        return;
    end
    
    if ~exist(jpgdir, 'dir')
        mkdir(jpgdir);
    end
   
    fprintf('Creating slices in %s\n', jpgdir);
    
    for i = 1:nImages
        close all
        hn=figure(1);
        if maxdim == 1
            tempimg = squeeze(img(i,:,:));
        elseif maxdim == 2
            tempimg = squeeze(img(:,i,:));
        elseif maxdim == 3
            tempimg = squeeze(img(:,:,i));
        end
        tempimg = fliplr(flipud(tempimg'));
        imagesc(tempimg), colormap gray
        set(gcf, 'Units', 'pixels', 'PaperPosition', [0 0 imgscale imgscale], 'Visible', 'off', 'Color', [0 0 0], 'inverthardcopy', 'off');
        set(gca, 'Units', 'inches', 'xtick',[],'ytick',[], 'Position', [0 0 imgscale imgscale]);
        text(5, size(tempimg,1)-10, num2str(i), 'color', 'yellow', 'FontSize', 7);
        print(hn,'-djpeg', jpgfile{i});
    end
    
    
    
    
function createMovie(niiFile, mcfile, outPath, imageType)


    noCols = 6;
    niiInfo = mrReadNifti(niiFile);
    img_raw = niiInfo.data;
    niiInfo = niiInfo.hdr;
    noSlices = niiInfo.dim(3);
    noFrames = niiInfo.dim(4);

    if strcmp(imageType, 'DTI')
        scH = 2.25;
        scW1 = 0.375;
        noSlabs = floor((noSlices-2)/ noCols);
        remSlices = mod(noSlices-2, noCols);
        slabs(:,1) = [2:noCols:noSlices-1]';
    else
        scH = 1.5;
        scW1 = 0.25;
        noSlabs = floor(noSlices / noCols);
        remSlices = mod(noSlices, noCols);
        slabs(:,1) = [1:noCols:noSlices-remSlices]';
    end
    
    if ~nnz(arrayfun(@(s) ~exist(fullfile(outPath, ['movie_raw_' num2str(s) '.gif']), 'file'), 1:noSlabs)) %% ~nnz(arrayfun(@(s) ~exist(fullfile(outpath, ['movie_mc_' num2str(s) '.gif']), 'file'), 1:noSlabs))
        fprintf('Already created\n');
        return;
    end    
    
%     info2 = niftiRead(mcfile);
%     img_mc = info2.data;
%     
%     if info1.sto_xyz(3,3) < 0
%         for fr = 1:size(img_raw, 4)
%             for z = 1:size(img_mc,3)
%                 img_mc(:,:,z,fr) = flipud(img_mc(:,:,z,fr));
%             end
%         end
%     end
    
    loops = 65535; %% loop forever
    delay = 0.25; %% seconds
    
    for s = 1:noSlabs-1
        slabs(s, 2) = slabs(s+1, 1) - 1;
    end
    slabs(end, 2) = slabs(end-1,2)+noCols + remSlices;

    for i = 1%:2
        if i == 1
            img = img_raw;
            token = 'raw';
        else
            img = img_mc;
            token = 'mc';
        end
        
        if ~matlabpool('size')
            matlabpool OPEN 6
        end
        
        parfor s = 1:noSlabs

            outDir = fullfile(outPath, ['montage' num2str(s) '_' token]);
            outFile = fullfile(outPath, ['movie_' token '_' num2str(s) '.gif']);
            if exist(outFile, 'file')
                continue;
            end

            if ~exist(outDir, 'dir')
                mkdir(outDir);
            end
            fprintf('Creating %s\n', outFile);

            jpgFile = cell(noFrames, 1);
            for fr = 1:noFrames

                jpgFile{fr} = fullfile(outDir, sprintf('Montage%d_%03d.jpg', s, fr));
                if ~exist(jpgFile{fr}, 'file') || 1
                    mntg = makeMontage(squeeze(img(:,:,slabs(s,1):slabs(s,2),fr)));
                    close all
                    hn=figure(1); 
                    imagesc(mntg), colormap gray
                    
                    scW = scW1 * (slabs(s,2)-slabs(s,1)+1);
                    set(gcf, 'Units', 'pixels', 'PaperPosition', [0 0 noCols*scW scH], 'Visible', 'off', 'Color', [0 0 0], 'inverthardcopy', 'off');
                    set(gca, 'Units', 'inches', 'xtick',[],'ytick',[], 'Position', [0 0 noCols*scW scH]);
                    text(2, size(mntg,1)-10, num2str(fr), 'color', 'yellow', 'FontSize', 12);
                    print(hn,'-djpeg', jpgFile{fr});
                end
            end

            for fr = 1:noFrames
                [M, cMap] = rgb2ind(imread(jpgFile{fr}),256);

                if fr==1
                    imwrite(M,cMap, outFile,'gif','LoopCount',loops,'DelayTime',delay);
                else
                    imwrite(M,cMap, outFile,'gif','WriteMode','append','DelayTime',delay);
                end
            end

            system(['rm -r ' outDir]);
        end
    end
    
    
    
    
function plotMotion(mcpath, outpath, imgdir, conn)

    outfile = fullfile(outpath, 'motionplot.html');

    if exist(outfile, 'file')
        return;
    end
    
    fprintf('Creating FSL Motion Plot\n');

    mot = textread([mcpath '.par']);
    mot_absrms = textread([mcpath '_abs.rms']);
    mot_relrms = textread([mcpath '_rel.rms']);

    meanrms_abs = textread([mcpath '_abs_mean.rms']);
    meanrms_rel = textread([mcpath '_rel_mean.rms']);
    
    datapaths = strsplit(imgdir, '/');
    if strcmp(datapaths{3}, 'DTI')
        count = cell2mat(mysqlQuery(conn, ['select count(*) from dti where imageCentralLocation = ''' imgdir '''']));

        if count
            mysqlQuery(conn, sprintf('UPDATE dti SET meanRMS_Absolute=%f, meanRMS_Relative=%f where imageCentralLocation=''%s'';', meanrms_abs, meanrms_rel, imgdir));
        else
        mysqlQuery(conn, sprintf('insert into dti (imageCentralLocation, meanRMS_Absolute, meanRMS_Relative, notes, include) values (''%s'',%f,%f, '''', 1);', imgdir, meanrms_abs, meanrms_rel));
        end
    end
    
    noVols = size(mot, 1);

    fid = fopen(outfile, 'w');
    fprintf(fid, '<html>\n');
    fprintf(fid, '<head>\n');
    fprintf(fid, '<script type="text/javascript" src="https://www.google.com/jsapi"></script>\n');
    fprintf(fid, '<script type="text/javascript">\n');
    fprintf(fid, 'google.load("visualization", "1", {packages:["corechart"]});\n');
    fprintf(fid, 'google.setOnLoadCallback(drawChart);\n');
    fprintf(fid, 'function drawChart() {\n');

    fprintf(fid, 'var data1 = google.visualization.arrayToDataTable([\n');
    fprintf(fid, '[''Volumes'', ''Trans x'', ''Trans y'', ''Trans z''],\n');
    for fr = 1:noVols-1
        fprintf(fid, '[%d,%f,%f,%f],\n', fr, mot(fr, 4), mot(fr, 5), mot(fr, 6));
    end
    fprintf(fid, '[%d,%f,%f,%f]\n]);\n\n', noVols, mot(end, 4), mot(end, 5), mot(end, 6));

    fprintf(fid, 'var data2 = google.visualization.arrayToDataTable([\n');
    fprintf(fid, '[''Volumes'', ''Rot x'', ''Rot y'', ''Rot z''],\n');
    for fr = 1:noVols-1
        fprintf(fid, '[%d,%f,%f,%f],\n', fr, mot(fr, 1), mot(fr, 2), mot(fr, 3));
    end
    fprintf(fid, '[%d,%f,%f,%f]\n]);\n\n', noVols, mot(noVols, 1), mot(noVols, 2), mot(noVols, 3));

    fprintf(fid, 'var data3 = google.visualization.arrayToDataTable([\n');
    fprintf(fid, '[''Volumes'', ''Absolute RMS'', ''Relative RMS''],\n');
    for fr = 1:noVols-2
        fprintf(fid, '[%d,%f,%f],\n', fr, mot_absrms(fr, 1), mot_relrms(fr, 1));
    end
    fprintf(fid, '[%d,%f,%f]\n]);\n\n', noVols, mot_absrms(noVols-1, 1), mot_relrms(noVols-1, 1));

    fprintf(fid, 'var options1 = {\n');
    fprintf(fid, 'title: ''Translation'', \n');
    fprintf(fid, 'titleTextStyle : {fontSize:25},\n');
    fprintf(fid, 'colors: [''red'', ''blue'', ''green''],\n');
    fprintf(fid, 'vAxis: {title: ''%s''}, \n', 'mm');
    fprintf(fid, 'hAxis: {title: ''%s''}\n', 'Volumes');
    fprintf(fid, '};\n\n');

    fprintf(fid, 'var options2 = {\n');
    fprintf(fid, 'title: ''Rotation'', \n');
    fprintf(fid, 'titleTextStyle : {fontSize:25},\n');
    fprintf(fid, 'colors: [''red'', ''blue'', ''green''],\n');
    fprintf(fid, 'vAxis: {title: ''%s''}, \n', 'rad');
    fprintf(fid, 'hAxis: {title: ''%s''}\n', 'Volumes');
    fprintf(fid, '};\n\n');

    fprintf(fid, 'var options3 = {\n');
    fprintf(fid, 'title: ''RMS'', \n');
    fprintf(fid, 'titleTextStyle : {fontSize:25},\n');
    fprintf(fid, 'colors: [''red'', ''blue'', ''green''],\n');
    fprintf(fid, 'hAxis: {title: ''%s''}\n', 'Volumes');
    fprintf(fid, '};\n\n');

    fprintf(fid, 'var chart1 = new google.visualization.LineChart(document.getElementById(''chart_div1''));\n');
    fprintf(fid, 'chart1.draw(data1, options1);\n');

    fprintf(fid, 'var chart2 = new google.visualization.LineChart(document.getElementById(''chart_div2''));\n');
    fprintf(fid, 'chart2.draw(data2, options2);\n');

    fprintf(fid, 'var chart3 = new google.visualization.LineChart(document.getElementById(''chart_div3''));\n');
    fprintf(fid, 'chart3.draw(data3, options3);\n');

    fprintf(fid, '}\n');

    fprintf(fid, '</script>\n');
    fprintf(fid, '</head>\n');
    fprintf(fid, '<body>\n');
    fprintf(fid, '<div id="chart_div3" style="width: 800px; height: 400px;"></div>\n');
    fprintf(fid, ['Mean RMS - Absolute : ' num2str(meanrms_abs)]);
    fprintf(fid, ['<br>Mean RMS - Relative : ' num2str(meanrms_rel) '<br>']);
    fprintf(fid, '<div id="chart_div1" style="width: 800px; height: 400px;"></div>\n');
    fprintf(fid, '<div id="chart_div2" style="width: 800px; height: 400px;"></div>\n');
    fprintf(fid, '</body>\n');
    fprintf(fid, '</html>\n');
    fclose(fid);
    
    
 
    
function [error, outfile] = fsl_motioncorrection(inpath, token)

    error = 0;
    
    infile = fullfile(inpath, [token '.nii.gz']);
    
    outdirname = 'FSL_MotionCorrect';
    outfile = fullfile(inpath, outdirname, [token, '_mc.nii.gz']);
    
    if exist(strrep(outfile, '.gz',''), 'file')
        error = system(strrep(outfile, '.gz',''));
        if error
            system(['rm ' strrep(outfile, '.gz','') ]);
        end
    elseif ~exist(outfile, 'file') || ...
       ~exist(strrep(outfile, 'nii.gz','par'), 'file') || ~exist(strrep(outfile, '.nii.gz','_rel.rms'), 'file') || ~exist(strrep(outfile, '.nii.gz', '_abs.rms'), 'file')
   
        fprintf('FSL Motion Correction : %s\n', inpath);
        if ~exist(fileparts(outfile), 'dir')
            mkdir(fileparts(outfile));
        end
        error = system(['/usr/local/fsl/bin/mcflirt -stats -plots -in ' infile ' -o ' strrep(outfile, '.nii.gz','') ' -rmsrel -rmsabs -refvol 0']);
        if error
            return;
        end
        system(['rm -r ' strrep(outfile, 'nii.gz','mat') ]);
    end

    
    
    
function [isDicom dcmfiles] = isDicomData(imgdir)    

    dcmfiles = dir(fullfile(imgdir, '*.dcm'));
    if isempty(dcmfiles)
        dcmfiles = dir(fullfile(imgdir, '*.DCM'));
    end

    if isempty(dcmfiles)
        isDicom = 0;
    else
        isDicom = 1;
    end
    

    
    
function montage_image = makeMontage(img)

    montage_image = [];
    for sl = 1:size(img,3)
        montage_image = [montage_image fliplr(flipud(img(:,:,sl)'))];
    end    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
% function createMontage(imgdir, outpath)
% 
%     if any(cellfun(@(x) ~isempty(strfind(lower(imgdir), x)), {'localizer', 'asset'}))
%         return;
%     end
% 
%     outfile = fullfile(outpath, 'montage.jpg');
%     
%     if exist(outfile, 'file')
%         return;
%     end
%     
%     anatfiles = dir(fullfile(imgdir, '*.dcm'));
%     if isempty(anatfiles)
%         anatfiles = dir(fullfile(imgdir, '*.DCM'));
%         if isempty(anatfiles)
%             return;
%         end
%     end
%     
%     fprintf('Creating Montage\n%s\n', outfile);
% 
%     anatfiles = {anatfiles.name}';
%     nImages = length(anatfiles);
%     info = dicominfo(fullfile(imgdir, anatfiles{1}));
%     
%     img = zeros(info.Width, info.Height, nImages);
%     imgMean = zeros(1, nImages);
%     for i = 1:nImages
%         img(:,:,i) = dicomread(fullfile(imgdir, anatfiles{i}));
%         tempimg = img(:,:,i);
%         imgMean(i) = mean(tempimg(:));
%     end
% 
%     %%% exclude slices without brain voxels
%     img = img(:, :, find(imgMean > 110));
% 
%     montage_image = makeMontage(img);
%     
%     nrows = size(montage_image,1) / info.Height;
%     imgscale = 1.5;
% 
%     close all
%     h=figure(1); 
%     imagesc(montage_image), colormap gray
%     set(gcf, 'Units', 'pixels', 'PaperPosition', [0 0 10 nrows*imgscale], 'Visible', 'off');
% 
%     set(gca, 'Units', 'inches', 'xtick',[],'ytick',[], 'Position', [0 0 10 nrows*imgscale]);
%     
%     if ~exist(fileparts(outfile), 'dir')
%         mkdir(fileparts(outfile));
%     end
%     print(h,'-djpeg', outfile);
%     close(h);
        









    
    
    
    
    
    
    
    
    
% %     imgscale = 0.5;
% %     jpgfile = cell(noFrames, 1);
% %     for fr = 1:noFrames
% % 
% %         jpgfile{fr} = fullfile(outpath, 'montage', sprintf('Montage%03d.jpg', fr)); 
% %         if ~exist(fileparts(jpgfile{fr}), 'dir')
% %             mkdir(fileparts(jpgfile{fr}));
% %         end
% %         
% %         jpgfile_raw = fullfile(outpath, 'montage_raw', sprintf('Montage%03d_raw.jpg', fr));
% %         if ~exist(fileparts(jpgfile_raw), 'dir')
% %             mkdir(fileparts(jpgfile_raw));
% %         end
% %         
% %         jpgfile_mc = fullfile(outpath, 'montage_mc', sprintf('Montage%03d_mc.jpg', fr));
% %         if ~exist(fileparts(jpgfile_mc), 'dir')
% %             mkdir(fileparts(jpgfile_mc));
% %         end
% %         
% %         if ~exist(jpgfile_raw, 'file') || ~exist(jpgfile_mc, 'file')
% %             
% %             mntg = makeMontage(squeeze(img_raw(:,:,:,fr)), 6);
% %             
% %             nrows = size(mntg, 1) / size(img_raw, 1);
% %             ncols = size(mntg, 2) / size(img_raw, 2);
% % 
% %             close all
% %             hn=figure(1); 
% %             imagesc(mntg), colormap gray
% %             set(gcf, 'Units', 'pixels', 'PaperPosition', [0 0 ncols*imgscale nrows*imgscale], 'Visible', 'off', 'Color', [0 0 0], 'inverthardcopy', 'off');
% %             set(gca, 'Units', 'inches', 'xtick',[],'ytick',[], 'Position', [0 0 ncols*imgscale nrows*imgscale]);
% %             text(5, size(mntg,1)-10, num2str(fr), 'color', 'yellow', 'FontSize', 7);
% %             print(hn,'-djpeg', jpgfile_raw);
% % 
% %             montage_mc = makeMontage(squeeze(img_mc(:,:,:,fr)), 6);
% % 
% %             close all
% %             hn=figure(1); 
% %             imagesc(montage_mc), colormap gray
% %             set(gcf, 'Units', 'pixels', 'PaperPosition', [0 0 ncols*imgscale nrows*imgscale], 'Visible', 'off', 'Color', [0 0 0], 'inverthardcopy', 'off');
% %             set(gca, 'Units', 'inches', 'xtick',[],'ytick',[], 'Position', [0 0 ncols*imgscale nrows*imgscale]);
% %             text(5, size(montage_mc,1)-10, num2str(fr), 'color', 'yellow', 'FontSize', 7);
% %             print(hn,'-djpeg', jpgfile_mc);
% %             
% %             montage_image = [mntg zeros(size(mntg,1), 50), montage_mc];
% % 
% %             nrows = size(montage_image, 1) / size(img_raw, 1);
% %             ncols = size(montage_image, 2) / size(img_raw, 2);
% % 
% %             close all
% %             hn=figure(1); 
% %             imagesc(montage_image), colormap gray
% %             set(gcf, 'Units', 'pixels', 'PaperPosition', [0 0 ncols*imgscale nrows*imgscale], 'Visible', 'off', 'Color', [0 0 0], 'inverthardcopy', 'off');
% %             set(gca, 'Units', 'inches', 'xtick',[],'ytick',[], 'Position', [0 0 ncols*imgscale nrows*imgscale]);
% %             text(5, size(montage_image,1)-10, num2str(fr), 'color', 'yellow', 'FontSize', 7);
% % 
% %             print(hn,'-djpeg', jpgfile{fr});
% %         end
% %     end
% % 
% %     for fr = 1:noFrames
% %         [M, cmap] = rgb2ind(imread(jpgfile{fr}),256);
% % 
% %         if fr==1
% %             imwrite(M,cmap,jpgfile,'gif','LoopCount',loops,'DelayTime',delay);
% %         else
% %             imwrite(M,cmap,jpgfile,'gif','WriteMode','append','DelayTime',delay);
% %         end
% %     end
% %     
% %     system(['rm -r ' fileparts(jpgfile{1})]);
%     
%     
% 
% %     jpgfile = cell(nFrames, 1);
% %     for fr = 1:nFrames
% % 
% %         jpgfile{fr} = fullfile(outpath, 'montage', sprintf('Montage%03d.jpg', fr)); 
% %         if ~exist(fileparts(jpgfile{fr}), 'dir')
% %             mkdir(fileparts(jpgfile{fr}));
% %         end
% %         
% %         jpgfile_raw = fullfile(outpath, 'montage_raw', sprintf('Montage%03d_raw.jpg', fr));
% %         if ~exist(fileparts(jpgfile_raw), 'dir')
% %             mkdir(fileparts(jpgfile_raw));
% %         end
% %         
% %         jpgfile_mc = fullfile(outpath, 'montage_mc', sprintf('Montage%03d_mc.jpg', fr));
% %         if ~exist(fileparts(jpgfile_mc), 'dir')
% %             mkdir(fileparts(jpgfile_mc));
% %         end
% %         
% %         if ~exist(jpgfile_raw, 'file') || ~exist(jpgfile_mc, 'file')
% %             
% %             montage_raw = makeMontage(squeeze(img_raw(:,:,:,fr)), 6);
% %             
% %             nrows = size(montage_raw, 1) / size(img_raw, 1);
% %             ncols = size(montage_raw, 2) / size(img_raw, 2);
% % 
% %             close all
% %             hn=figure(1); 
% %             imagesc(montage_raw), colormap gray
% %             set(gcf, 'Units', 'pixels', 'PaperPosition', [0 0 ncols*imgscale nrows*imgscale], 'Visible', 'off', 'Color', [0 0 0], 'inverthardcopy', 'off');
% %             set(gca, 'Units', 'inches', 'xtick',[],'ytick',[], 'Position', [0 0 ncols*imgscale nrows*imgscale]);
% %             text(5, size(montage_raw,1)-10, num2str(fr), 'color', 'yellow', 'FontSize', 7);
% %             print(hn,'-djpeg', jpgfile_raw);
% % 
% %             montage_mc = makeMontage(squeeze(img_mc(:,:,:,fr)), 6);
% % 
% %             close all
% %             hn=figure(1); 
% %             imagesc(montage_mc), colormap gray
% %             set(gcf, 'Units', 'pixels', 'PaperPosition', [0 0 ncols*imgscale nrows*imgscale], 'Visible', 'off', 'Color', [0 0 0], 'inverthardcopy', 'off');
% %             set(gca, 'Units', 'inches', 'xtick',[],'ytick',[], 'Position', [0 0 ncols*imgscale nrows*imgscale]);
% %             text(5, size(montage_mc,1)-10, num2str(fr), 'color', 'yellow', 'FontSize', 7);
% %             print(hn,'-djpeg', jpgfile_mc);
% %             
% %             montage_image = [montage_raw zeros(size(montage_raw,1), 50), montage_mc];
% % 
% %             nrows = size(montage_image, 1) / size(img_raw, 1);
% %             ncols = size(montage_image, 2) / size(img_raw, 2);
% % 
% %             close all
% %             hn=figure(1); 
% %             imagesc(montage_image), colormap gray
% %             set(gcf, 'Units', 'pixels', 'PaperPosition', [0 0 ncols*imgscale nrows*imgscale], 'Visible', 'off', 'Color', [0 0 0], 'inverthardcopy', 'off');
% %             set(gca, 'Units', 'inches', 'xtick',[],'ytick',[], 'Position', [0 0 ncols*imgscale nrows*imgscale]);
% %             text(5, size(montage_image,1)-10, num2str(fr), 'color', 'yellow', 'FontSize', 7);
% % 
% %             print(hn,'-djpeg', jpgfile{fr});
% %         end
% %     end
% % 
% %     for fr = 1:nFrames
% %         [M, cmap] = rgb2ind(imread(jpgfile{fr}),256);
% % 
% %         if fr==1
% %             imwrite(M,cmap,outfile,'gif','LoopCount',loops,'DelayTime',delay);
% %         else
% %             imwrite(M,cmap,outfile,'gif','WriteMode','append','DelayTime',delay);
% %         end
% %     end
% %     
% %     system(['rm -r ' fileparts(jpgfile{1})]);
% 
% 
% % function montage_image = makeMontage(img, ncols)
% % 
% % %     if ~exist('isDicom', 'var') || isempty(isDicom)
% % %         isDicom = 1;
% % %     end
% % 
% %     if ~exist('ncols', 'var') || isempty(ncols)
% %         ncols = 6;
% %     end
% % % 
% % %     if ~exist('flip', 'var')
% % %         flip = 0;
% % %     end
% % 
% %     nImages = size(img, 3);
% %     width = size(img, 1);
% %     height = size(img, 2);
% %     nrows = ceil(nImages / ncols);
% % 
% %     montage_image = [];
% %     cnt = nImages;
% %     for r = 1:nrows
% %         rowIm = [];
% %         for c = 1:ncols
% %             n = (r-1)*ncols + c; % index into image
% %             if n<=nImages
% %                 rowIm = [rowIm fliplr(flipud(img(:,:,n)'))];
% %             else
% %                 rowIm = [rowIm zeros(width, height)];
% %             end
% %         end
% %         montage_image = [montage_image; rowIm];
% %     end    
% %     
%     
% %     if isDicom
% %         for r = 1:nrows
% %             rowIm = [];
% %             for c = 1:ncols
% %                 n = (r-1)*ncols + c; % index into image
% %                 if n<=nImages
% % %                     if flip
% % %                         rowIm = [rowIm flipud(img(:,:,n)')];
% % %                     else
% %                         rowIm = [rowIm img(:,:,n)];
% % %                     end
% %                 else
% %                     rowIm = [rowIm zeros(width, height)];
% %                 end
% %             end
% %             montage_image = [montage_image; rowIm];
% %         end
% % 
% %     elseif ~isDicom
% %         cnt = nImages;
% %         for r = nrows:-1:1
% %             rowIm = [];
% %             for c = ncols:-1:1
% %                 if cnt > 0
% %                     rowIm = [rowIm flipud(img(:,:, cnt)')];
% %                 else
% %                     rowIm = [rowIm zeros(width, height)];
% %                 end
% %                 cnt = cnt - 1;
% %             end
% %             montage_image = [montage_image; rowIm];
% %         end
% %     end


    
% % %     anatfiles = dir(fullfile(imgdir, '*.dcm'));
% % %     if isempty(anatfiles)
% % %         anatfiles = dir(fullfile(imgdir, '*.DCM'));
% % %         if isempty(anatfiles)
% % %             return;
% % %         end
% % %     end
% % %     
% % %     anatfiles = {anatfiles.name}';
% % %     info = dicominfo(fullfile(imgdir, anatfiles{1}));
% % % 
% % %     nImages = length(anatfiles);
% % %     
% % %     img = zeros(info.Width, info.Height, nImages);
% % %     imgMean = zeros(1, nImages);
% % %     for i = 1:nImages
% % %         img(:,:,i) = dicomread(fullfile(imgdir, anatfiles{i}));
% % %         tempimg = img(:,:,i);
% % %         imgMean(i) = mean(tempimg(:));
% % %     end
% % % 
% % %     %%% exclude slices without brain voxels
% % %     img = img(:, :, find(imgMean > 110));
% % %     imgscale=3.5;
% % %     
% % %     nImages = size(img,3);
% % %     
% % %     jpgfile = arrayfun(@(i) fullfile(jpgdir, sprintf('slice%03d.jpg', i)), 1:nImages, 'UniformOutput', false);
% % %     if ~any(~cellfun(@exist, jpgfile))
% % %         return;
% % %     end
% % %     
% % %     fprintf('Creating slices in %s\n', jpgdir);
% % %     
% % %     for i = 1:nImages
% % %         close all
% % %         hn=figure(1); 
% % %         imagesc(img(:,:,i)), colormap gray
% % %         set(gcf, 'Units', 'pixels', 'PaperPosition', [0 0 imgscale imgscale], 'Visible', 'off', 'Color', [0 0 0], 'inverthardcopy', 'off');
% % %         set(gca, 'Units', 'inches', 'xtick',[],'ytick',[], 'Position', [0 0 imgscale imgscale]);
% % %         text(5, size(img,1)-10, num2str(i), 'color', 'yellow', 'FontSize', 7);
% % %         print(hn,'-djpeg', jpgfile{i});
% % %     end
    

