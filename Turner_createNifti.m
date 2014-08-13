function Turner_createNifti(newIDs)

    clc;

    %%% Input directory (data can be in either location)
    dataDir1 = '/Volumes/Projects/TURNER/Data/sortedData';
    dataDir2 = '/Volumes/ImageCentral';
    
    %%% Output directory
    niiDir = '/Volumes/Projects/TURNER/Data/nifti/';
    
    %%% Output directory name to store FSL motion corrected data
    fslDir = 'FSL_MotionCorrect';
    
    %%% Set the matlab environment variables(fsl location, toolsmac )
    setMatlabEnv;

    %%% Output file prefix (e.g. f.nii.gz)
    outToken = 'f';

    %%% Connect to mysql database
    conn = mysqlConnect;
    
    %%% If no IDs specified (newIDs), select all
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
    
    %%% Loop through each subject
    for s = 1:length(subjectIDs)

        %%% Prints to output
        fprintf('\n%s\nSubject ID : %d\nSession ID : %s\nExam ID : E%d\n%s\n\n', repmat('*', 1, 50), subjectIDs{s}, sessionIDs{s}, examIDs{s}, repmat('*', 1, 50));

        %%% Get the location of data for current subject
        imgDirs = mysqlQuery(conn, sprintf('select imagecentralLocation from imaging where subjectID=%d and sessionID=''%s'' and examID=%d', subjectIDs{s}, sessionIDs{s}, examIDs{s}));
        
        %%% Loop through
        for i = 1:length(imgDirs)
            
            %%% Check if input is Structural, Functional or DTI
            dataPaths = strsplit(imgDirs{i}, '/');
            
            dataType = dataPaths{3}(1);
            
            %%% For Anatomical images
            if strcmp(dataPaths{3}, 'Structural')
                
                %%% Skip creating nifti if Localizer or Asset
                if any(cellfun(@(x) strcmpi(dataPaths{5}(5:end), x), {'localizer', 'asset'}))
                    continue;
                end
                
                %%% Output directory to save nifti
                outDir = fullfile(dataPaths{1}, dataPaths{2}, dataPaths{4}, dataPaths{5});
             
            %%% For DTI/Functional images
            elseif ismember(dataPaths{3}, {'DTI', 'Functional'})

                %%%
                examID = num2str(cell2mat(mysqlQuery(conn, sprintf('select examID from imaging where imageCentralLocation=''%s''', imgDirs{i}))));
                outDir = fullfile(dataPaths{1}, dataPaths{2}, ['E' examID], strrep(dataPaths{4}, ['_E' examID], ''));
                
            end
            
            if exist(fullfile(dataDir1, num2str(subjectIDs{s})), 'dir')
                inDir = dataDir1;
            else
                inDir = fullfile(dataDir2, ['Subjs_' num2str(floor(subjectIDs{s}/100)*100)]);
            end
            
            [error, result, niiFile] = createNifti(imgDirs{i}, inDir, fullfile(niiDir, outDir), outToken, dataType, conn);
            
            if ~result
                fprintf('%s\n', error);
                continue;
            end
            
            %%% Run FSL Motion Correction (not for structural)
            if ismember(dataType, {'F', 'D'}) && length(dir(fullfile(niiDir, outDir, fslDir, '*.rms'))) ~= 4
                fslMotionCorrection(fullfile(niiDir, outDir), fslDir, outToken);
            end

            %%% Plot mrVista motion parameters (functional only)
            if strcmp(dataType, 'F') && nnz(cellfun(@(x) ~isempty(strfind(niiFile, x)), {'CategoryLocalizer', 'MeridianMapping', 'Eccentricity', 'MTLocalizer'}))
                mrVistaMotionCorrection(niiFile, imgDirs{i}, conn);
            end
        end
    end

    
    

    
%%% Function to create nifti files    
function [error, result, zniiFile] = createNifti(icLoc, rawDir, niiDir, outToken, dataType, conn)

    error = 1;
    result = 0;
    
    %%% Input directory to dicom/spiral files
    rawDir = fullfile(rawDir, icLoc);

    %%% Output file
    niiFile = fullfile(niiDir, [outToken '.nii']);
    zniiFile = [niiFile '.gz'];
    
    %%% If the uncompressed file exists, compress it (f.nii -> f.nii.gz)
    if exist(niiFile, 'file')
        system(['gzip ' niiFile]);
        result = 1;
        error = 0;
    end
    
    %%% Create output nifti directory
    if ~exist(niiDir, 'dir')
        mkdir(niiDir);
    end
    
    if ~exist(zniiFile, 'file')

        fprintf('\nCreating nifti : %s\n', zniiFile);
        
        
        %%% Depending on input data, call respective function
        switch dataType
            
            %%% Use niftiFromDicom to convert structural dicom files
            case 'S'
                tempfile = niftiFromDicom(rawDir, niiDir);
                system(['mv ' tempfile ' ' zniiFile]);
                system(['cp ' fullfile(rawDir, 'cibsrDcm2fmp.csv') ' ' niiDir]);
                
            %%% Use dcm2nii to convert DTI files
            case 'D'
                [~, dcmFiles] = isDicomData(rawDir);
                [error, result] = system(['dcm2nii -o ' niiDir ' ' fullfile(rawDir, dcmFiles(1).name)]);

                if error || ~isempty(strfind(result, 'ERROR'))
                    return;
                end

                tempfile = regexp(result, '->.*.nii\nGZip', 'match');
                system(['mv ' fullfile(niiDir, [tempfile{1}(3:end-5) '.gz']) ' ' niiFile '.gz']);
                system(['rm ' fullfile(niiDir, ['*o' tempfile{1}(3:end-5) '.gz']) ' 2> /dev/null ']);

            %%% Use Gary's scripts to convert spiral files
            case 'F'
            
                eFile = mysqlQuery(conn, sprintf('select pfileName from imaging where imageCentralLocation=''%s''', icLoc));
                eFile = dir(fullfile(rawDir, ['E*' eFile{1}]));

                %%% Use Gary's script to convert spiral data to nifti (default options)
                [error, ~] = system(['makenifti ' fullfile(rawDir, eFile.name) ' ' strrep(niiFile, '.nii', '')]);

                if error
                    return;
                end

                [~, nVols] = system(['fslnvols ' niiFile]);

                %%% GZip the nifti file
                system(['gzip -f ' niiFile]);

                %%% Remove first and last 4 timepoints from ECC, MM
                curDir = pwd;
                cd(niiDir);
                if str2num(nVols) == 104 && nnz(cellfun(@(x) ~isempty(strfind(niiDir, x)), {'MeridianMapping', 'Eccentricity'}))
                    system(['fslsplit ' outToken]);
                    for v = [0:3 100:103]
                        system(sprintf('rm vol%04d.nii.gz', v));
                    end

                    for v = 0:95
                        system(sprintf('mv vol%04d.nii.gz vol%04d.nii.gz', v+4, v));
                    end
                    system(['fslmerge -t ' outToken ' vol*']);
                    system('rm vol*');
                end
                cd(curDir);
        end
    else
        fprintf('%s exists\n', zniiFile);
        error = 0;
    end

    result = 1;



    
    
%%% Checks if input if of type DICOM       
function [isDicom dcmfiles] = isDicomData(inDir)    

    dcmfiles = dir(fullfile(inDir, '*.dcm'));
    if isempty(dcmfiles)
        dcmfiles = dir(fullfile(inDir, '*.DCM'));
    end

    if isempty(dcmfiles)
        isDicom = 0;
    else
        isDicom = 1;
    end
    
    

    
%%% Runs within run motion correction, using first volume as reference    
function [error, outFile] = fslMotionCorrection(inDir, fslDir, inToken)

    error = 0;
    
    inFile = fullfile(inDir, [inToken '.nii.gz']);
    outFile = fullfile(inDir, fslDir, [inToken, '_mc.nii.gz']);
    
    if exist(outFile(1:end-3), 'file')
        error = system(strrep(outFile, '.gz',''));
        if error
            system(['rm ' strrep(outFile, '.gz','') ]);
        end
    elseif ~exist(outFile, 'file') || ...
       ~exist(strrep(outFile, 'nii.gz','par'), 'file') || ~exist(strrep(outFile, '.nii.gz','_rel.rms'), 'file') || ~exist(strrep(outFile, '.nii.gz', '_abs.rms'), 'file')
   
        fprintf('FSL Motion Correction : %s\n', inDir);
        if ~exist(fileparts(outFile), 'dir')
            mkdir(fileparts(outFile));
        end
        [error, ~] = system(['mcflirt -stats -plots -in ' inFile ' -o ' strrep(outFile, '.nii.gz','') ' -rmsrel -rmsabs -refvol 0 >/dev/null']);
        if error
            return;
        end
        system(['rm -r ' strrep(outFile, 'nii.gz','mat') ]);
    end

    
    
    
function mrVistaMotionCorrection(niiFile, icLoc, conn)
    
    plotFile = fullfile(fileparts(niiFile), 'mrVista_motion.png');
     
    if exist(plotFile, 'file')
        return;
    end
    
    fprintf('Plotting %s\n', plotFile);
    
    motionParams = Turner_mrVistaWithinMotionCorrection(niiFile);
    
    %% report on motion; record motion estimates in dataTYPES
    hFig = figure;
    t = 1:size(motionParams, 2);
    plot(t, motionParams(1, :), t, motionParams(2, :), t, sqrt(sum(motionParams(1:2, :).^2))); drawnow;
    legend('Translational', 'Rotational', 'Total', 'Location', 'Best');

    results = mysqlQuery(conn, sprintf('select seriesNo, seriesName from imaging where imageCentralLocation=''%s''', icLoc));
    seriesNo = results{1};
    seriesName = results{2};

    title(sprintf('Motion : %03d %s', seriesNo, strrep(seriesName, '_', ' ')));
    xlabel('Time (frames)');
    ylabel('Motion (voxels)');
    
    print(hFig, '-dpng', plotFile);
    close(hFig);
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
% %     
% %     
% % %%% NOT USED    
% % function outlierDetection(infile)
% % 
% % %     infile = '/Users/harshads/Projects/Turner/Analysis/nifti/18501/13-09-05.3_3T2/E12571/005_HARDI/f.nii.gz'; 
% %     %18226/13-03-04.1_3T3/E2524/006_HARDI/f.nii.gz'; 
% %     %13135/13-09-17.1_3T2/E12639/010_HARDI/f.nii.gz';
% %     info = readFileNifti(infile);
% %     
% %     dir_img = info.data;
% %     
% %     b0_info = info;
% %     b0_info.dim = [128 128 50 1];
% %     b0_info.fname = fullfile(fileparts(infile), 'B0_1.nii.gz');
% %     b0_info.data = dir_img(:,:,:, 1:26:156);
% %     writeFileNifti(b0_info);
% %     
% %     setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
% %     system(['bet ' fullfile(fileparts(infile), 'B0_1.nii.gz') ' ' fullfile(fileparts(infile), 'B0_1_brain.nii.gz')]);
% %     
% %     b0_img = readFileNifti(fullfile(fileparts(infile), 'B0_1_brain.nii.gz'));
% %     b0_img = b0_img.data;
% %     slices = find(mean(mean(b0_img)) > 150);
% %     
% %     %% using only gradient direction images, and not B0 images
% %     dir_img = dir_img(:,:,:,setdiff(1:156, 1:26:156));
% %     
% % %     temp_img = zeros(128,128,1);
% %     directions = cell(150, 1);
% %     
% %     for sl = 1:length(slices)
% %        meanval = squeeze(mean(mean(dir_img(:,:,slices(sl),:)))); 
% %        
% %        q1 = quantile(meanval, .25);
% %        q3 = quantile(meanval, .75);
% %        
% %        ll = q1 - 1.5*iqr(meanval);
% %        ul = q3 + 1.5*iqr(meanval);
% %        
% %        outlierNos = union(find(meanval < ll), find(meanval > ul));
% %        if ~isempty(outlierNos)
% % %           fprintf('Outliers detected for slice %d at %s\n', slices(sl), sprintf('%-4d', outlierNos)); 
% % %           temp_img = cat(3, temp_img, squeeze(dir_img(:,:,slices(sl),outlierNos)));
% % %           directions = unique([directions(:); outlierNos(:)]);
% %           for o = 1:length(outlierNos)
% %               if isempty(directions{outlierNos(o)})
% %                   directions{outlierNos(o)} = sl;
% %               else
% %                   directions{outlierNos(o)} = [directions{outlierNos(o)} sl];
% %               end
% %           end
% %        end       
% %     end
% %     
% %     for d = 1:150
% %         if length(directions{d})
% %            fprintf('Direction %d\tSlices : %s\n', d + ceil(d/25), sprintf(' %d', directions{d}));
% %         end
% %     end
% %     
% %     
% %     
    
    
    
    
