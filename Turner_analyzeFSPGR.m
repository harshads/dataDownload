%%% Runs the mrVista ACPC realign module, saves output as f_acpc.nii.gz (in
%%% the same folder as f.nii.gz)
%%% Copies acpc aligned file and the corresponding cibsrDcm2fmp.csv file to the
%%% Freesurfer/gbb2/forStart folder.
%%% Runs SPM8 Very Light Bias Regularization on the acpc realigned images
%%% and saves bias corrected nifti


function Turner_analyzeFSPGR(subjectID)

    %%% Check the user running the script
    userName = char(java.lang.System.getProperty('user.name'));
    if ~strcmp(userName, 'harshads')
        fprintf('Not added for other users\n');
        return;
    end
    
    %%% Add path to SPM8 revision5236
    addpath /Volumes/ToolsMac/SPM/spm8_r5236/
    
    %%% Initialize spm_jobman gui
    spm_jobman('initcfg');

    %%% nifti directory
    dataDir = '/Volumes/Projects/TURNER/Data/nifti';
    
    %%% SPM8 template file to run bias regularization
    spmTemplateFile = '/Volumes/Projects/TURNER/Data/Files/SPM8_BiasRegularization_Template.mat';
    
    %%% output directory (freesurfer)
    fsDir = '/Volumes/Xspace/Freesurfer/gbb2linux';
    
    conn = mysqlConnect;

    %%% If no subjectID was specified, run for all subjects, check if FSPGR
    %%% was marked usable (>0)
    if ~exist('subjectID', 'var') 
        imgDirs = mysqlQuery(conn, 'select imageCentralLocation from anatomical where usability > 0');
    else
        imgDirs = mysqlQuery(conn, ['select imageCentralLocation from anatomical where usability > 0 and imageCentralLocation like ''' subjectID '%''']);
    end
    
%     %%% all subject IDs (cibsrIDs), can have multiple if multiple FSPGRs
%     %%% were marked usable
%     all_subjectIDs = strtok(imgDirs, '/');
    
    for i = 1:length(imgDirs)

        inFile = fullfile(dataDir, strrep(imgDirs{i}, 'Structural/', ''), 'f.nii.gz');
        subjectID = strtok(imgDirs{i}, '/');
        acpcFile = fullfile(fileparts(inFile), 'f_acpc.nii.gz');
        biasFile = fullfile(fileparts(inFile), 'mf_acpc.nii.gz');

        %%% Run ACPC realign
        if ~exist(acpcFile, 'file')
            fprintf('\n\nACPC Realign - %s\n', inFile);
            mrAnatAverageAcpcNifti(inFile, acpcFile);
            
            %%% Close ACPC realign windows
            close(findobj('Type', 'figure', 'Name',  'Average'));
            close(findobj('Type', 'figure', 'Name',  'Average Aligned'));
            close(findobj('Type', 'figure', 'Name',  [inFile ' (ref)']));
        end
        
        %%% Copy cibsrDcm2fmp.csv to nifti folder
        csvFile = fullfile(fileparts(inFile), 'cibsrDcm2fmp.csv');
        if ~exist(csvFile, 'file')
            in_csvFile1 = fullfile('/Volumes/Projects/TURNER/Data/sortedData/', imgDirs{i}, 'cibsrDcm2fmp.csv');
            in_csvFile2 = fullfile(['/Volumes/ImageCentral/Subjs_' num2str(floor((str2num(subjectID)/100))*100)], imgDirs{i}, 'cibsrDcm2fmp.csv');
            
            if ~exist(in_csvFile1, 'file') && ~exist(in_csvFile2, 'file')
                Turner_copyFMPcsv(subjectID);
            end
            
            if exist(in_csvFile1, 'file') 
                system(['cp ' in_csvFile1 ' ' csvFile]);
            elseif exist(in_csvFile2, 'file') 
                system(['cp ' in_csvFile2 ' ' csvFile]);
            end
        end
        
        %%% Run SPM8 New Segment Bias Regularization
        if ~exist(biasFile, 'file')
            system(['gunzip ' acpcFile]);
            load(spmTemplateFile);
            matlabbatch{1}.spm.tools.preproc8.channel.vols = {strrep(acpcFile, '.gz', '')};
            spm_jobman('run_nogui', matlabbatch);
            system(['gzip ' strrep(acpcFile, '.gz', '')]);
            system(['gzip ' strrep(biasFile, '.gz', '')]);
        end
        
%         info1 = niftiRead(biasfile);
%         info2 = niftiRead(['/Volumes/Xspace/Freesurfer/gbb2linux/rawdata/fromLinux/m' subjectID '.nii']);
%         if mean(info1.data(:) - info2.data(:)) < 0.01
%             continue;
%         end

        timePoint = cell2mat(mysqlQuery(conn, sprintf('select Timepoint from imaging where imageCentralLocation=''%s''', imgDirs{i})));
        allFiles = mysqlQuery(conn, sprintf('select imageCentralLocation from anatomical where usability <> 0 and imageCentralLocation like ''%s%%'' order by imageCentralLocation', subjectID));
        tempFiles = mysqlQuery(conn, sprintf('select imageCentralLocation from imaging where Timepoint=%d and imageCentralLocation like ''%s%%''', timePoint, subjectID));
        allFiles =  intersect(allFiles, tempFiles);
        
        outToken = sprintf('m%s.%d', subjectID, timePoint);
        if length(allFiles) > 1
           outToken = sprintf('%s_%d',  outToken, strmatch(imgDirs{i}, allFiles));
        end
        
        outData{1} = fullfile(fsDir, 'forStart', [outToken '.nii.gz']);
        outData{2} = fullfile(fsDir, 'rawdataDone', [outToken '.nii.gz']);
        outData{3} = fullfile(fsDir, 'finishSuccess', [outToken '_5.3l']);
        outData{4} = fullfile(fsDir, 'finishSuccess', [outToken '_5.3l3T']);
        
        if ~nnz(cellfun(@exist, outData))            
            fprintf('\ncp %s %s\n',biasFile,  outData{1});
            system(['cp ' biasFile ' ' outData{1}]);
            system(['cp ' csvFile ' ' strrep(outData{1}, '.nii.gz', '_cibsrDcm2fmp.csv')]);
        end
    end



