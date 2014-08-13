clear
clc

addpath /Volumes/ToolsMac/SPM/spm8_r5236/    
addpath /Volumes/ToolsMac/SPM/fMRItools/

setMatlabEnv;

% spm('defaults', 'fmri');
% spm_jobman('initcfg');

conn = mysqlConnect;

datadir = '/Volumes/Projects/TURNER/Data/nifti';

if ~isempty(conn.Message)
    fprintf('Could not connect to database\n%s\n', conn.Message);
    return;
end

%%% Motion Statistics for Functional Data (CatLoc, MM, ECC, MTLoc)
% [imgCentralLocations, examIDs] = mysql('select functional.imageCentralLocation, imaging.examID from functional join imaging on functional.imageCentralLocation=imaging.imageCentralLocation order by functional.imageCentralLocation desc');
queryResult = mysqlQuery(conn, 'select functional.imageCentralLocation, imaging.examID from functional join imaging on functional.imageCentralLocation=imaging.imageCentralLocation order by functional.imageCentralLocation desc');
imgCentralLocations = queryResult(:, 1);
examIDs = cell2mat(queryResult(:, 2));

for i = 1:length(imgCentralLocations)
    
    queryResult = mysqlQuery(conn, sprintf('select * from functional where imageCentralLocation = ''%s''', imgCentralLocations{i}));

    if ~nnz(cellfun(@isnan, queryResult(:,end-10:end)))
        continue;
    end

    indir = fullfile(datadir, strrep(imgCentralLocations{i}, 'Functional', ['E' num2str(examIDs(i))]));
    %%% Sometimes, functional names have examIDs as suffix, if multiple
    %%% examIDs have functional runs with same series No.
    indir = strrep(indir, ['_E' num2str(examIDs(i))], '');
    infile = fullfile(indir, 'f.nii.gz');
    
    fprintf('\n\nAnalyzing %s\n', infile);
    
    %%% FSL Motion Outliers
    fsldir = fullfile(indir, 'FSL_MotionCorrect');
    
    %%% Temp Code
    reAnalyze = 0;
    if exist(fullfile(fsldir, 'f_mo_refrms'), 'file')
        temp=load(fullfile(fsldir, 'f_mo_refrms'));
        if length(temp) == 104 && nnz(cellfun(@(x) ~isempty(strfind(indir, x)), {'MeridianMapping', 'Eccentricity'}))
           reAnalyze = 1;
        end
    end
    
    if ~exist(fullfile(fsldir, 'f_mo_refrms'), 'file') || reAnalyze
        fprintf('FSL Motion Outliers : %s\n', infile);
        system(['fsl_motion_outliers -i ' infile ' -o ' fullfile(fsldir, 'f_mo_refrms > /dev/null')]);
    end
    
    noVols = size(textread(fullfile(fsldir, 'f_mo_refrms')), 1);
    fsl_pctmo = size(textread(fullfile(fsldir, 'f_mo_refrms')), 2) / noVols;
    
    fsl_absrms_max = max(textread(fullfile(fsldir, 'f_mc_abs.rms')));
    fsl_absrms_mean = mean(textread(fullfile(fsldir, 'f_mc_abs.rms')));
    fsl_relrms_max = max(textread(fullfile(fsldir, 'f_mc_rel.rms')));
    fsl_relrms_mean = mean(textread(fullfile(fsldir, 'f_mc_rel.rms')));

    %%% ART Motion params
    spmdir = fullfile(indir, 'SPM_MotionCorrect');

    if ~exist(fullfile(spmdir, 'rf.nii.gz'), 'file') || reAnalyze
        
        if ~exist(spmdir, 'dir')
            mkdir(spmdir);
        end
        
        system(['cp ' infile ' ' spmdir]);
        infile = strrep(infile, indir, spmdir);
        system(['gunzip -f ' infile]);
        load /Volumes/Projects/TURNER/Data/Files/SPM8_MotionCorrect.mat
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {arrayfun(@(x) fullfile([strrep(infile, '.gz', '') ',' num2str(x)]), 1:noVols, 'UniformOutput', false)'};
        spm_jobman('run', matlabbatch);
        system(['rm ' strrep(infile, '.gz','')]);
        system(['gzip -f ' fullfile(spmdir, '*.nii')]);
    end

    stats = art_motionstats({spmdir}, spmdir, spmdir);
    system(['rm ' fullfile(spmdir, 'MotionChars*.txt')]);
    stats(5) = stats(5)/noVols;
    stats(6) = stats(6)/noVols;
    
    sqlstr = sprintf('FSL_motionOutliersPct = %.4f, FSL_absoluteRMSMax = %.4f, FSL_absoluteRMSMean = %.4f, FSL_relativeRMSMax = %.4f, FSL_relativeRMSMean = %.4f%s', ...
        fsl_pctmo, fsl_absrms_max, fsl_absrms_mean, fsl_relrms_max, fsl_relrms_mean, ...
        cell2mat(arrayfun(@(s) sprintf(', ART_MotionStats%d = %.4f', s, stats(s)), 1:6, 'UniformOutput', false)));
    
    sqlstr = sprintf('UPDATE functional SET %s where imageCentralLocation=''%s''', sqlstr, imgCentralLocations{i});
    
    fprintf('%s\n', sqlstr);
    mysqlQuery(conn, sqlstr); 

end
















% 
% 
% 
% %%% Motion Statistics for Resting State Data (CatLoc, MM, ECC, MTLoc)
% [imgCentralLocations, examIDs] = mysql('select resting.imageCentralLocation, imaging.examID from resting join imaging on resting.imageCentralLocation=imaging.imageCentralLocation order by resting.id');
% 
% conn = database('turner', 'root', 'toortoor', '/Library/Java/Extensions/mysql-connector-java-5.1.27-bin.jar', 'jdbc:mysql://localhost/turner');
% 
% for i = 1:length(imgCentralLocations)
%     
%     curs = exec(conn, sprintf('select * from resting where imageCentralLocation = ''%s''', imgCentralLocations{i}));
%     results = fetch(curs);
% 
%     if ~nnz(cellfun(@isnan, results.Data(:,end-10:end)))
%         continue;
%     end
% 
%     indir = fullfile(datadir, strrep(imgCentralLocations{i}, 'Functional', ['E' num2str(examIDs(i))]));
%     indir = strrep(indir, ['_E' num2str(examIDs(i))],'');
%     infile = fullfile(indir, 'f.nii.gz');
%     
%     fprintf('Analyzing %s\n', infile);
%     
%     %%% FSL Motion Outliers
%     fsldir = fullfile(indir, 'FSL_MotionCorrect');
%     if ~exist(fullfile(fsldir, 'f_mo_refrms'), 'file')
%         fprintf('FSL Motion Outliers : %s\n', infile);
%         system(['fsl_motion_outliers -i ' infile ' -o ' fullfile(fsldir, 'f_mo_refrms')]);
%     end
%     
%     noVols = size(textread(fullfile(fsldir, 'f_mo_refrms')), 1);
%     fsl_pctmo = size(textread(fullfile(fsldir, 'f_mo_refrms')), 2) / noVols;
%     
%     fsl_absrms_max = max(textread(fullfile(fsldir, 'f_mc_abs.rms')));
%     fsl_absrms_mean = mean(textread(fullfile(fsldir, 'f_mc_abs.rms')));
%     fsl_relrms_max = max(textread(fullfile(fsldir, 'f_mc_rel.rms')));
%     fsl_relrms_mean = mean(textread(fullfile(fsldir, 'f_mc_rel.rms')));
% 
%     %%% ART Motion params
%     spmdir = fullfile(indir, 'SPM_MotionCorrect');
%     if ~exist(fullfile(spmdir, 'rf.nii.gz'), 'file')
%         
%         if ~exist(spmdir, 'dir')
%             mkdir(spmdir);
%         end
%         
%         system(['cp ' infile ' ' spmdir]);
%         infile = strrep(infile, indir, spmdir);
%         system(['gunzip -f ' infile]);
%         load /Volumes/Projects/TURNER/Data/Files/SPM8_MotionCorrect.mat
%         matlabbatch{1}.spm.spatial.realign.estwrite.data = {arrayfun(@(x) fullfile([strrep(infile, '.gz', '') ',' num2str(x)]), 1:noVols, 'UniformOutput', false)'};
%         spm_jobman('run', matlabbatch);
%         system(['rm ' strrep(infile, '.gz','')]);
%         system(['gzip -f ' fullfile(spmdir, 'rf.nii')]);
%         system(['gzip -f ' fullfile(spmdir, 'meanf.nii')]);
%     end
% 
%     stats = art_motionstats({spmdir}, spmdir, spmdir);
%     system(['rm ' fullfile(spmdir, 'MotionChars*.txt')]);
%     stats(5) = stats(5)/noVols;
%     stats(6) = stats(6)/noVols;
%     
%     sqlstr = sprintf('FSL_pctMotionOutliers = %.4f, FSL_AbsoluteRMS_Max = %.4f, FSL_AbsoluteRMS_Mean = %.4f, FSL_RelativeRMS_Max = %.4f, FSL_RelativeRMS_Mean = %.4f%s', ...
%         fsl_pctmo, fsl_absrms_max, fsl_absrms_mean, fsl_relrms_max, fsl_relrms_mean, ...
%         cell2mat(arrayfun(@(s) sprintf(', ART_MotionStats%d = %.4f', s, stats(s)), 1:6, 'UniformOutput', false)));
%     
%     mysql(sprintf('UPDATE resting SET %s where imageCentralLocation=''%s''', sqlstr, imgCentralLocations{i}));
% 
% end
