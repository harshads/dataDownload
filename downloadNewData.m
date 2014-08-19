function downloadNewData

    %%% Download Turner data from ScanParking to local folders %%%

    clc;
    
    %%% global variables 
    global MOVEFLAG createdDirs taskNames taskTokens

    %%% Set the matlab environment
    setMatlabEnv;

    %%% Check if connection to server exists
    try
        sortdir = '/Volumes/Projects/TURNER/Data/sortedData'; %'/Users/harshads/Projects/Turner/Data/';%
        cd(sortdir);

        parkingdir = '/Volumes/scanParking/STAGING/Turner/';
        cd(parkingdir);

    catch
        fprintf('Cannot cd into network drives\n');
        return;
    end

    %%% get the cibsrIDs in scanParking
    subjectIDs = getSubjectIDs(parkingdir);
    
    %%% Set some variables
    
    %%% Input session folder should match the pattern
    scanner = {'3T1', '3T2', '3T3'};
    studyNo = {'\.[A-Z0-9][A-Z0-9][A-Z0-9]_', '\.[A-Z][A-Z]_', '\.[0-9]_'};
    pattern = cellfun(@(y) cellfun(@(x) [y x], scanner, 'UniformOutput', false), studyNo, 'UniformOutput', false);
    pattern = [pattern{:}];
    
    %%% Output folder name for behavioral data
    behDir = 'Functional/Behavioral';
    
    %%% Retinotopy task names
    taskNames = {'CategoryLocalizer', 'MeridianMapping', 'Eccentricity', 'MTLocalizer'};
    
    %%% Corresponding names for behavioral mat files
    taskTokens = {'CatLoc', 'mm_cw', 'ecc_out', 'mtloc'};
    
    %%% Flag not used
    MOVEFLAG = 1;
    
    all_outNames = {};
    
    %%% Loop through the IDs
    for i = 1:length(subjectIDs)
        
        %%% Current output subject directory
        subjectDir = fullfile(sortdir, subjectIDs{i});
        
        %%% Current input directory in scan parking
        cur_parkingDir = fullfile(parkingdir, subjectIDs{i});
        
        %%% Check if existing subject directory already has sessions 
        excludedDirs = getDirs(subjectDir);
        
        %%% Only need to copy new session data
        newDirs = setdiff(getDirs(cur_parkingDir), excludedDirs);
        
        if isempty(newDirs)
            continue;
        end
        
        fprintf('\n%s\nSUBJECT ID : %s\n%s\n', repmat('*', 1, 50), subjectIDs{i}, repmat('*', 1, 50));
        
        curDir = pwd;
        cd(cur_parkingDir);
        newDirs = checkBehavioral(newDirs);
        cd(curDir);
        
        if isempty(newDirs)
            fprintf('MISSING BEHAVIORAL DATA\n');
            continue;
        end
        
        %%% Check if any of the folders match the required pattern
        %%% names (eg. YY-MM-DD.NO_3T2 etc)
        sessionDirs = {};
        for p = pattern;
            sessionDirs = [sessionDirs newDirs(~cellfun(@isempty, regexp(newDirs, p{1})))];
            newDirs = setdiff(newDirs, sessionDirs);
            if isempty(newDirs)
                break;
            end
        end
        
        %%% Check if folder did not match required format
        if ~isempty(newDirs)
            fprintf('THE FOLLOWING SESSIONS DO NOT MATCH REQUIRED NAMING FORMAT\n');
            fprintf('CANNOT MOVE THE FOLLOWING FOLDERS\n');
            cellfun(@(x) fprintf('%s\n', fullfile(subjectDir, x)), newDirs);
        end
        
        if isempty(sessionDirs)
            fprintf('NO DATA TO SORT\n');
            continue;
        end

        %%% Copy the scan data from scanParking to sorting folder
        if ~exist(subjectDir, 'dir')
            mkdir(subjectDir);
        end
        
        copyData = cellfun(@(x) sprintf('cp -r %s %s', fullfile(cur_parkingDir, x), subjectDir), sessionDirs, 'UniformOutput', false);
        fprintf('%s\n', copyData{:});
        tic, cellfun(@(x) system(x), copyData), toc;
        
        %%% Helps to remove .DS_Store files
        system(['find ' subjectDir ' -type f -name ".DS_Store" -exec rm {} \;']);

        %%% Remove empty cells
%         sessionDirs(cellfun(@isempty, sessionDirs)) = [];

        %%% For every new session 
        for s = 1:length(sessionDirs)
        
            %%% CD into new session directory
            cd(fullfile(subjectDir, sessionDirs{s}));
            
            %%% Get the scan information for all dicom and spiral files
            scanInfo = getInformation(sessionDirs{s});
            
            %%% Variables to hold information about data that will be
            %%% sorted, screen saves or unsorted (for duplicate or
            %%% incomplete data)
            unsortedInfo = {};
            sortedInfo = {};
            screensaveInfo = {};
            
            %%% Keeps track of folders the script creates
            createdDirs = {};
                        
            %%% Loop through the information
            for sc = 1:length(scanInfo)
                    
                if scanInfo{sc}.sort
                    if scanInfo{sc}.isduplicate
                        scanInfo{sc}.dest = 'Unsorted/Duplicate';
                        unsortedInfo{end+1} = scanInfo{sc};
                    elseif strcmp(scanInfo{sc}.scantype, 'Screen Save')
                        screensaveInfo{end+1} = scanInfo{sc};
                    elseif ~scanInfo{sc}.iscomplete
                        scanInfo{sc}.dest = 'Unsorted/Incomplete';
                        unsortedInfo{end+1} = scanInfo{sc};
                    else
                        sortedInfo{end+1} = scanInfo{sc};
                    end
                else
                    scanInfo{sc}.dest = '../Unsorted';
                    unsortedInfo{end+1} = scanInfo{sc};
                end
            end
            
            %%% If we found data that needs to be sorted
            if length(sortedInfo)
                
                sortedInfo = sortData(sortedInfo);
                moveFiles('Screenshots', {'.png', '.jpg', '.tiff'});
                moveOtherScreenshots(screensaveInfo);
                moveFiles(behDir, {'.edat', '.edat2'});
                moveFiles(behDir, {'.mat'}, sortedInfo);
                sortedInfo = movePsychBehavioral(sortedInfo, behDir);
                sortEyeTrack(sortedInfo, behDir);
            end
               
            %%% Sort the incomplete/duplicate/other data
            if length(unsortedInfo)
                moveUnsorted(unsortedInfo);
            end
            
            %%% Clean up empty folders
            dataCleanUp;
            
            %%% Paths to sorted data
            outNames = cellfun(@(x) fullfile(subjectIDs{i}, sessionDirs{s}, x.filename), sortedInfo, 'UniformOutput', false);
            
            %%% Create cibsrDcm2FMP.csv file to upload to FileMakePro
            cd(sortdir);
            sqlInfo = createFMPcsv(outNames);

            %%% Update sql database   
            if ~isempty(sqlInfo)
                update_sqlDB(sqlInfo);
            end

            %%% Marks the folder as sorted
            system(['touch ' fullfile(sortdir, subjectIDs{i}, sessionDirs{s}, '.sorted')]);

            %%% Cell array for further steps (creating nifti, creating jpegs).
            all_outNames = [all_outNames; outNames'];
            
        end
        
        %%% Writes the individual csv files
        cd(sortdir);
        Turner_copyFMPcsv(subjectIDs{i});

    end

    if isempty(all_outNames)
        return;
    end
    
    %%% Get a list of subjectIDs analyzed in this run
    subjectIDs = unique(strtok(all_outNames, '/'));
    
    %%% Create nifti data for these subjects
    Turner_createNifti(subjectIDs);
    
    %%% Update study information for these subjects
    Turner_updateSQL;

    %%% Runs for user 'harshads'
    curUser = char(java.lang.System.getProperty('user.name'));
    if strcmp(curUser, 'harshads')
        Turner_generateImages(subjectIDs);
    end

    
    
function newDirs = checkBehavioral(newDirs)

    global taskTokens
    
    goodIndex = [];
    names = {'catloc', 'mm', 'ecc', 'mtloc'};

    for n = 1:length(newDirs)
        
        [~, funcFiles] = system(['find ' newDirs{n} ' -type f -name "E*.7"']);
        funcFiles = strsplit(funcFiles, char(10));
        if ~length(funcFiles)
            continue;
        end
        
        examNos = cellfun(@(f) strtok(getSpiralField(f, 'examnum/seriesnum'), '/'), funcFiles, 'UniformOutput', false);
        seriesNames = cellfun(@(f) strtok(getSpiralField(f, 'series description'), '/'), funcFiles, 'UniformOutput', false);

        matFlag = 1;
        
        for exNo = sort(unique(examNos))
            
            exNo = num2str(str2num(regexprep(exNo{1}, '[^0-9]', '')));
            
            for tNm = names

                curIndex = find(cellfun(@(x) ~isempty(strfind(x(1:6), tNm{1})), seriesNames));

                noFuncFiles = length(curIndex);
                if noFuncFiles

                    matToken = taskTokens{find(cellfun(@(x) ~isempty(strfind(lower(x), tNm{1})), taskTokens))};

                    [~, behFiles] = system(['find ' newDirs{n} ' -type f -name "' exNo '*' matToken '*.mat"']);
                    noBehFiles = nnz(~cellfun(@isempty, strsplit(behFiles, char(10))));
                    
                    if noBehFiles ~= noFuncFiles
                        matFlag = 0;
                    end
                end
            end
        end
        
        if matFlag 
            goodIndex(end+1) = n;
        end
    end
    
    newDirs = arrayfun(@(n) newDirs(n), goodIndex);
    

%%% Read required line (token) from the spiral E*.7 text file    
function match_str = getSpiralField(efile, token)

    token = strrep(token, ' ', '\ ');
    [~, match_str] = system(['grep -r ' token ' ' efile '| cut -d ''='' -f2 | sed ''s/^[ \t]*//;s/[ \t]*$//''']);
    match_str = strtrim(match_str);

    

    
%%% Returns information about each dicom and spiral file found
%%% Fields include filename, scanner, scandate, pfile, examno
%%% seriesno, scantype, sessionid iscomplete
function out_scanInfo = getInformation(sessiondir)

    %%% List of first dicom file and Spiral E*.7 file
    [~, alldirs] = system('find . -depth -type d | grep -v Unsorted | sed ''s|./||'' | sed ''s|$|linebreak|''');
    alldirs = strsplit(alldirs, 'linebreak')';
    alldirs = alldirs(~cellfun(@isempty, alldirs));

    datafiles = cell(length(alldirs), 1);
    for d = 1:length(alldirs)
        tmpfiles = [dir(fullfile(alldirs{d}, '*.dcm')); dir(fullfile(alldirs{d}, '*.DCM'))];
        tmpfiles = {tmpfiles.name}';
        if ~isempty(tmpfiles)
            datafiles{d} = fullfile(alldirs{d}, tmpfiles{1});
        end
    end

    datafiles(cellfun(@isempty, datafiles)) = [];

    [~, tmpfiles] = system('find . -type f -iname ''E*.7'' | sort | sed ''s|./||''');

    if ~isempty(tmpfiles)
        tmpfiles = strsplit(tmpfiles, char(10));
        [~, idx] = sort(cellfun(@(x) x(max(strfind(x, '/'))+1:end), tmpfiles, 'UniformOutput',false));
        datafiles = [datafiles; tmpfiles(idx)'];
    end

    fprintf('\n%s Found these files %s\n', repmat('-', 1, 15), repmat('-', 1, 15));
    for d = 1:length(datafiles)
        fprintf('%s\n', datafiles{d});
    end
        
    session_scandate = datestr(datenum(sessiondir(1:8), 'yy-mm-dd'), 'yyyymmdd');
    session_scanner = sessiondir(max(strfind(sessiondir, '_'))+1:end);
        
    out_scanInfo = {};
    no_ExamSeries = {};
    
    for d = 1:length(datafiles)
        
        curfile = datafiles{d};
        
        scanInfo.filename = curfile;
        scanInfo.scandate = '';
        scanInfo.pfile = 'none';
        scanInfo.examno = '';
        scanInfo.seriesno = '';
        scanInfo.scantype = '';
        scanInfo.sessionid = sessiondir;
        scanInfo.sort = 0;
        scanInfo.iscomplete = 0;
        scanInfo.isduplicate = 0;
        scanInfo.imagetype = 'dicom';
        scanInfo.timestamp = '';
        
        %%%%%%%%%%%%%% Spiral file %%%%%%%%%%%%%%
        if strcmp(curfile(end-1:end), '.7')
            
            scanInfo.imagetype = 'spiral';
            
            [examNo, seriesNo] = strtok(getSpiralField(curfile, 'examnum/seriesnum'), '/');
            seriesNo = seriesNo(2:end);
            
            %%% Set the ExamNo, SeriesNo & ScanType
            scanInfo.examno = ['E' num2str(str2num(examNo))];
            scanInfo.seriesno = str2num(seriesNo(2:end));
            scanInfo.scantype = 'Functional';
            
            %%% Set the ScanDate
            scanDate = getSpiralField(curfile, 'date of scan');
            scanDate = scanDate([1:6 8:end]);
            scanInfo.scandate = datestr(datenum(scanDate, 'mm/dd/yy'), 'yyyymmdd');
            
            %%% Split into path & file name
            [oPath oFile oExt] = fileparts(curfile);
            oFile = [oFile oExt];
            
            %%% Set the P filename
            scanInfo.pfile = cell2mat(regexp(oFile, 'P\d+.7', 'match'));
            
            %%% Complete path to P file
            pFile = fullfile(oPath, scanInfo.pfile);
            
            %%% If mag/hdr files exist, check the size matches the required size
            if exist([pFile '.hdr'], 'file') && exist([pFile '.mag'], 'file')
                try
                    noSlices = str2num(getSpiralField(curfile, 'slquant'));
                    noTimePoints = str2num(getSpiralField(curfile, 'num time frames'));
                    imgSize = str2num(getSpiralField(curfile, 'imgsize'));

                    magSize = dir([pFile '.mag']);
                    magSize = magSize.bytes;
                    if magSize/noSlices/noTimePoints/imgSize^2 == 4
                        scanInfo.iscomplete = 1;
                    end
                catch err
                    fprintf('ERROR reading %s\n', curfile);
                end

                try
                    scannerID = getSpiralField(curfile, 'scanner id');
    
                    if ismember(scannerID, {'No1', 'No2', 'No3'})
                        scannerID = strrep(scannerID, 'No', '3T');
                    end
                catch err
                    scannerID = 'Unknown';
                end
            end
            
            scanInfo.sort = strcmp(scanInfo.scandate, session_scandate) && strcmp(scannerID, session_scanner);
            
            scanInfo.timestamp = getFuncTimestamp(curfile);
            
        %%%%%%%%%%%%%% Dicom file %%%%%%%%%%%%%%
        elseif strcmpi(curfile(end-3:end), '.dcm')
            
            %%% Read the dicom header
            dcmInfo = dicominfo(curfile);
            
            %%% Dicoms for screensaves are missing information, will be caught by exception
            try
                
                if ismember(dcmInfo.InstitutionName, {'Stanford Lucas Center', 'Lucas Center'})
                    magnetNo = strrep(dcmInfo.PerformedStationName, 'No', '');
                    scannerID = sprintf('%dT%s', dcmInfo.MagneticFieldStrength, magnetNo);
                    siteName = 'Lucas';
                elseif strcmp(dcmInfo.InstitutionName, 'CNI')
                    scannerID = sprintf('CNI%dT', dcmInfo.MagneticFieldStrength);
                    siteName = 'CNI';
                end
                
                scanInfo.scantype = getDicomScanType(curfile, siteName);
                
                [dcmPath, ~, dcmExt] = fileparts(curfile);
                numFiles = length(dir(fullfile(dcmPath, ['*' dcmExt])));
                
                if numFiles == dcmInfo.ImagesInAcquisition
                    scanInfo.iscomplete = 1;
                end    
                
            catch err
                
                %%% Dicoms are screen saves
                scanInfo.scantype = 'Screen Save';
                scannerID = session_scanner;
                
                %%% Mark as complete
                scanInfo.iscomplete = 1;
            end
                
            scanInfo.scandate = dcmInfo.StudyDate;
            scanInfo.examno = sprintf('E%d', str2num(dcmInfo.StudyID));
            scanInfo.seriesno = dcmInfo.SeriesNumber;
                  
            scanInfo.sort = strcmp(scanInfo.scandate, session_scandate) && strcmp(scannerID, session_scanner);
            
        end    
        
        no_ExSe = sprintf('%s_%d', scanInfo.examno, scanInfo.seriesno);
        
        if ismember(no_ExSe, no_ExamSeries)
            scanInfo.isduplicate = 1;
        elseif scanInfo.iscomplete
            no_ExamSeries{end+1} = no_ExSe;
        end

        out_scanInfo{end+1} = scanInfo;    
        
        
    end 
        

    
    
function curDirs = getDirs(indir)

    %%% Find all files in current directory
    curDirs = dir(indir);
    
    %%% Exclude first 2 ('.', '..'), and check if others are directories
    try
        curDirs = curDirs([false false arrayfun(@(x) curDirs(x).isdir, 3:length(curDirs))]);
        curDirs = sort({curDirs.name})';
    catch err
        curDirs = {};
    end
    
    
    
%%% Find CIBSR IDs in Turner scanParking folder    
function subjectIDs = getSubjectIDs(datadir)

    subjectIDs = getDirs(datadir);
    subjectIDs = subjectIDs(cellfun(@(x) ~isnan(str2double(x)) && isdir(fullfile(datadir, x)), subjectIDs));

    
    
    
%%% Called by Dicom Files only
function scanType = getDicomScanType(infile, Site)

    seriesName = getSeriesName(infile);
    scanType = 'Structural';
    dcmInfo = dicominfo(infile);
    
    if strcmp(seriesName, 'HARDI') || ~isempty(strfind(lower(dcmInfo.SeriesDescription), 'dti'))
        scanType = 'DTI';
    elseif strcmp(Site, 'Lucas')
       if ~isempty(strfind(lower(dcmInfo.Private_0019_109c), 'dti')) || ~isempty(strfind(lower(dcmInfo.Private_0019_109e), 'dti_epi'))
           scanType = 'DTI';
       end
    end
   
    
    

%%% Get the seriesName
function out_seriesName = getSeriesName(infile)

    study_seriesNames = setSeriesNames;
    
    try
        info = dicominfo(infile);
        in_seriesName = info.SeriesDescription;
    catch err
        in_seriesName = getSpiralField(infile, 'series description');
    end
    
    try
        out_seriesName = study_seriesNames(in_seriesName);
    catch err
        out_seriesName = regexprep(in_seriesName, '[^a-zA-Z0-9]', '');
    end
    
    

    
%%% Set the seriesNames based on seriesNames in Pfile or Dicom files
function study_seriesNames = setSeriesNames
        
    study_seriesNames = containers.Map;
    study_seriesNames('loc') = 'Localizer';
    study_seriesNames('3plnloc') = 'Localizer';
    study_seriesNames('ASSET calibration') = 'ASSET';
    study_seriesNames('FSPGR_sagittal') = 'FSPGR';
    study_seriesNames('HARDI 150 R=2') = 'HARDI';
    study_seriesNames('InplaneAnatomy') = 'InplaneAnatomy';
    study_seriesNames('mtloc_39sl_3mm_19.2_tr2.5') = 'MTLocalizer';
    study_seriesNames('Resting_30sl_Tr2') = 'Resting';
    study_seriesNames('catloc_39sl_3mm_19.2_tr2.4') = 'CategoryLocalizer';
    study_seriesNames('mm1_39sl_3mm_19.2') = 'MeridianMapping';
    study_seriesNames('ecc1_39sl_3mm_19.2') = 'Eccentricity';
    
    
    
%%% Sorts all the imaging files
function sortInfo = sortData(sortInfo)

    global MOVEFLAG createdDirs

    fprintf('\n%s Moving data %s\n', repmat('-', 1, 15), repmat('-', 1, 15));

    for sc = 1:length(sortInfo)

        filename = sortInfo{sc}.filename;
        seriesNo = sortInfo{sc}.seriesno;
        examNo = sortInfo{sc}.examno;
        
        seriesName = getSeriesName(filename);
        seriesFolder = sprintf('%03d_%s', seriesNo, seriesName);
        
        if ~isempty(strfind(filename, '/'))
            indir = filename(1:max(strfind(filename, '/'))-1);
            oFile = filename(max(strfind(filename, '/'))+1:end);
        else
            indir = '.';
            oFile = filename;
        end
            
        if strcmp(sortInfo{sc}.scantype, 'Structural')
            tempdir = examNo;
        else
            tempdir = '';
        end
        
        outdir = fullfile(sortInfo{sc}.scantype, tempdir, seriesFolder);
        
        %%% If input and output are same, 
        if strcmp(indir, outdir)
            createdDirs{end+1} = outdir;
            fprintf('Source and Destination are the same %s\n', outdir);
            continue;
        end
        
        %%% Check if output folder exists; will add a suffix of examNo
        if exist(outdir, 'dir')
            outdir = [outdir '_' examNo];
        end
        
        %%% This should not happen, but in case it does
        %%% Duplicate directory
        tempdir = outdir;
        count = 1;
        moveflag = 1;
        while exist(tempdir, 'dir')
            tempdir = sprintf('%s_%d', outdir, count);
            if strcmp(tempdir, indir)
                moveflag = 0;
                break;
            count = count + 1;
            end
        end
            
        outdir = tempdir;
            
        %%% Dicom file
        if strcmp(sortInfo{sc}.imagetype, 'dicom')
            
            fprintf('MOVING %s TO %s\n', indir, outdir);
            if MOVEFLAG
                outPath = fileparts(outdir);
                if ~exist(outPath, 'dir')
                    mkdir(outPath);
                end
                system(sprintf('mv "%s" %s', indir, outdir));
            end
            
        %%% Spiral file
        elseif strcmp(sortInfo{sc}.imagetype, 'spiral')

            %%% Get the P filename
            pfile = sortInfo{sc}.pfile(1:end-2);
            
            infiles = fullfile(indir, ['*' pfile '*']);
                        
            fprintf('MOVING %s TO %s\n', infiles, outdir);
            if MOVEFLAG
                if ~exist(outdir, 'dir')
                    mkdir(outdir);  
                end
                system(sprintf('mv %s %s', infiles, outdir));
            end
        end
        
        createdDirs{end+1} = outdir;
        
        sortInfo{sc}.filename = fullfile(outdir, oFile);
        
    end
    
    
    
    
function moveFiles(outdir, pattern, scanInfo)
    
    global MOVEFLAG createdDirs
    
    fprintf('\n%s Moving %s %s\n', repmat('-', 1, 15), outdir, repmat('-', 1, 15));
    
    examNos = [];
    if exist('scanInfo', 'var')
        examNos = sort(unique(cellfun(@(s) s.examno, scanInfo, 'UniformOutput', false)));
    end
    
    files = {};
    for p = pattern
        if length(examNos)
            for e = 1:length(examNos)
                [~, tempfiles] = system(['find . -type f -name "*' examNos{e}(2:end) '*' p{1} '" | sed ''s|./||'' | sed ''s|$|linebreak|''']);
                tempfiles = strsplit(tempfiles, 'linebreak');
                tempfiles = tempfiles(~cellfun(@isempty, tempfiles));
                files = [files tempfiles];
            end
        else
            [~, tempfiles] = system(['find . -type f -name "*' p{1} '" | sed ''s|./||'' | sed ''s|$|linebreak|''']);
            tempfiles = strsplit(tempfiles, 'linebreak');
            tempfiles = tempfiles(~cellfun(@isempty, tempfiles));
            files = [files tempfiles];
        end
    end
    
    %%% For every matching file
    for f = 1:length(files)
        
        infile = files{f};
        
        if strcmp(strtok(infile(3:end), '/'), outdir)
            continue;
        end
        
        outfile = renameFile(outdir, infile);
    
        fprintf('MOVING %s TO %s\n', infile, outfile);
        if MOVEFLAG
            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end
            
            createdDirs{end+1} = outdir;
            system(sprintf('mv "%s" "%s"\n', infile, outfile));
        end
        
        infile = [infile(1:max(strfind(infile, '.'))-1) '.txt'];
        if exist(infile, 'file')
            outfile = [outfile(1:max(strfind(outfile, '.'))-1) '.txt'];
            fprintf('MOVING %s TO %s\n', infile, outfile);
            if MOVEFLAG
                system(sprintf('mv "%s" "%s"', infile, outfile));
            end
        end
    end
    
    
    
    
function outfile = renameFile(outdir, infile)
    
    [~, fName, fExt] = fileparts(infile);
    fName = [fName fExt];
    outfile = fullfile(outdir, fName);
    count = 1;
    while exist(outfile, 'file')
        idx = max(strfind(fName, '.'));
        outfile = fullfile(outdir, sprintf('%s_%d%s', fName(1:idx-1), count, fName(idx:end)));
        count = count + 1;
    end

    
    
    
function moveOtherScreenshots(screensaveInfo)
    
    global MOVEFLAG createdDirs
    
    fprintf('\n%s Moving Screenshots %s\n', repmat('-', 1, 15), repmat('-', 1, 15));
    
    picdir = 'Screenshots';
    
    %%% Screenshots folder
    for ss = 1:length(screensaveInfo)
        
        ssInfo = screensaveInfo{ss};
        
        %%% Input directory
        infile = ssInfo.filename;
        indir = fileparts(infile);
    
        %%% Dicom folder name
        inFolder = strsplit(indir, '/');
        inFolder = inFolder{end};
        
        if exist(fullfile(picdir, inFolder), 'dir')
            inFolder = [inFolder '_' ssInfo.examno];
        end
        
        outdir = checkDir(inFolder, picdir);
        
        fprintf('MOVING %s TO %s\n', indir, outdir);
        if MOVEFLAG
            if ~exist(picdir, 'dir')
                mkdir(picdir);
            end
            createdDirs{end+1} = outdir;
            system(sprintf('mv "%s" "%s"', indir, outdir));
        end
    end

    
    
    
function outdir = checkDir(inFolder, indir)

    cnt = 1;
    while exist(fullfile(indir, inFolder), 'dir')
        inFolder = [inFolder '_' num2str(cnt)];
        cnt = cnt + 1;
    end
        
    outdir = fullfile(indir, inFolder);
    
    
    
    
function scanInfo = movePsychBehavioral(scanInfo, behDir)

    global MOVEFLAG createdDirs taskNames taskTokens
    
    bExt = '.mat';
    
    examNos = cellfun(@(s) s.examno, scanInfo, 'UniformOutput', false);
    infiles = cellfun(@(s) s.filename, scanInfo, 'UniformOutput', false);
    timestamps = cellfun(@(s) s.timestamp, scanInfo, 'UniformOutput', false);
    
    for t = 1:length(taskNames)
        
        %%% Current task
        tName = taskNames{t};
        
        %%% Get a list of unique examIDs
        exNos = sort(unique(examNos));
        
        %%% Find all files for current task
        idx1 = find(cellfun(@(x) ~isempty(strfind(x, tName)), infiles));
        
        pattern = taskTokens{t};
        
        %%% For every examID
        for e = 1:length(exNos)
            
            %%% Find the index for current examID
            idx2 = strmatch(exNos{e}, examNos);
            
            %%% Final index of current task for current examID
            idx = intersect(idx1, idx2);
        
            funcFiles = infiles(idx);   
            funcTS = cell2mat(timestamps(idx));
            
            behFiles = dir(fullfile(behDir, [exNos{e}(2:end) '_*' pattern '*' bExt]));
            
            if ~isempty(behFiles) && length(funcFiles) == length(behFiles)
                
                behFiles = cellfun(@(x) fullfile(behDir, x), {behFiles.name}', 'UniformOutput', false);
                
                behTS = cellfun(@(x) datenum(regexp(x, '\d{8}_\d{4}', 'match'), 'mmddyyyy_HHMM'), behFiles);
                [~, idx] = sort(behTS);
                behFiles = behFiles(idx);
                
                [funcTS, idx] = sort(funcTS);
                funcFiles = funcFiles(idx);
                
                for b = 1:length(behFiles)
                    
                    [indir, fName, fExt] = fileparts(funcFiles{b});
                    fName = [fName fExt];
                    outdir = indir;
                    
                    if ~strcmp(tName, 'MTLocalizer')
                        
                        [~, bName, ~] = fileparts(behFiles{b});
                        runNo = strtok(bName(length(exNos{e})+1:end), '_');
                        outdir = [outdir '_' runNo];
                        
                        fprintf('MOVING %s TO %s\n', indir, outdir);
                        if MOVEFLAG
                            system(sprintf('mv %s %s', indir, outdir));
                            createdDirs{strmatch(indir, createdDirs)} = outdir;
                        end
                        
                        scanInfo{strmatch(funcFiles{b}, infiles)}.filename = fullfile(outdir, fName);
                    end
                    
                    infiles{find(cellfun(@(x) ~isempty(strfind(x, indir)), infiles))} = fullfile(outdir, fName);
                    
                    outToken = outdir(max(strfind(outdir, '/'))+1:end);
                    [bPath, bName] = fileparts(behFiles{b});
                    bName = [outToken '__' bName bExt];
                    outfile = fullfile(bPath, bName);
                    fprintf('MOVING %s TO %s\n', behFiles{b}, outfile);
                    if MOVEFLAG
                        system(sprintf('mv %s %s', behFiles{b}, outfile));
                    end
                end
            end
        end
    end


    
    
function ftimestamp = getFuncTimestamp(func_file)

    ftime = regexprep(getSpiralField(func_file, 'time of scan'), '[^0-9:]', '');
    fdate = regexprep(getSpiralField(func_file, 'date of scan'), '[^0-9//]', '');
    
    ftimestamp = datenum([fdate([1:6 8:end]) ':' ftime], 'mm/dd/yy:HH:MM');
    

    
    
function dataCleanUp    
    
    global createdDirs
        
    tempDirs = {};
    for c = 1:length(createdDirs)
        cDir = createdDirs{c};
        idx = max(strfind(cDir, '/'));
        while ~isempty(idx)
            cDir = cDir(1:idx-1);
            tempDirs{end+1} = cDir;
            idx = max(strfind(cDir, '/'));
        end
    end
    
    createdDirs = unique([createdDirs tempDirs]);
    
    system('find . -type f -name ".DS_Store" -exec rm {} \; 2>/dev/null');
    system('find . -type f -name ANAT_XFER_DONE -exec rm {} \; 2>/dev/null');
    system('find . -type f -name SUMMARY -exec rm {} \; 2>/dev/null');
    system('find . -type d -depth -exec rmdir {} \; 2>/dev/null');
    
    grepstr = sprintf(' | grep -v %s', createdDirs{:});
    [~, move_dirs] = system(sprintf('find . -maxdepth 1 -type d %s | sed ''s|./||'' | sed ''s|$|linebreak|''', grepstr));
    move_dirs = strsplit(move_dirs, 'linebreak')';
    move_dirs = move_dirs(~cellfun(@isempty, move_dirs));
    
    outdir = 'Unsorted';
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    if ~isempty(move_dirs)
        for m = 1:length(move_dirs)
            if ~strcmp(move_dirs{m}, '.')
                system(sprintf('mv "%s" %s\n', move_dirs{m}, outdir));
            end
        end
    end
            
    [~, move_files] = system('find . -maxdepth 1 -type f | sed ''s|./||'' | sed ''s|$|linebreak|''');
    move_files = strsplit(move_files, 'linebreak')';
    move_files = move_files(~cellfun(@isempty, move_files));
    
    if ~isempty(move_files)
        for m = 1:length(move_files)
            system(sprintf('mv "%s" %s\n', move_files{m}, outdir));
        end
    end

    %%% Final clean up in case anything was missed
    system('find . -type f -name ".DS_Store" -exec rm {} \; 2>/dev/null');
    system('find . -type d -depth -exec rmdir {} \; 2>/dev/null');

    system('rmdir Unsorted 2>/dev/null');
    
    
    
    
function sortEyeTrack(scanInfo, behDir)

    global MOVEFLAG createdDirs taskNames
    
    fprintf('\n%s Moving Eyetracking %s\n', repmat('-', 1, 15), repmat('-', 1, 15));
    
    etFiles = {};
    etTimeStamps = [];
    etDurations = [];
    
    [~, tempfiles] = system(['find . -type f -name "*.txt" | grep -v ' behDir ' | sed ''s|./||'' | sed ''s|$|linebreak|''']);
    tempfiles = strsplit(tempfiles, 'linebreak');
    tempfiles = tempfiles(~cellfun(@isempty, tempfiles));
    
    %%% Find the first dicom file in every sub directory
    for t = 1:length(tempfiles)
        
        curfile = tempfiles{t};
        
        fid = fopen(curfile, 'r');
        lines = {};
        while ~feof(fid)
            lines{end+1} = fgetl(fid);
        end
        fclose(fid);
        
        if nnz(cellfun(@(l) ~isempty(strfind(l, 'PupilWidth')), lines))
            etFiles{end+1} = curfile;
            %%% Split string using tab - char(9)
            ts = strsplit(lines{find(cellfun(@(l) ~isempty(strfind(l, 'TimeStamp')), lines))}, char(9));
               
            etTimeStamps(end+1) = datenum([ts{4:end}], 'mmmdd,yyyy,HH:MM:SSAM');
            
            dur = strsplit(lines{end});
            etDurations(end+1) = str2double(dur{2});
        end
    end
            
    %%% Sort the Behavioral files based on Timestamp
    
    [etTimeStamps, idx] = sort(etTimeStamps);
    etFiles = etFiles(idx);
    etDurations = etDurations(idx);
    
    funcFiles = {};
    funcTimeStamps = [];
    funcDurations = [];
    
    for sc = 1:length(scanInfo)
        
        curfile = scanInfo{sc}.filename;
        
        if nnz(cellfun(@(t) ~isempty(strfind(curfile, t)), taskNames))
            
            funcFiles{end+1} = curfile;
        
            funcTimeStamps(end+1) = scanInfo{sc}.timestamp;

            noFrames = str2num(getSpiralField(curfile, 'num time frames'));
            repTime = str2double(strrep(getSpiralField(curfile, 'TR'), 'msec', ''));
            funcDurations(end+1) = noFrames*repTime/1000;
        end
    end
    
    %%% Sort the Behavioral files based on Timestamp
    [funcTimeStamps, idx] = sort(funcTimeStamps);
    funcFiles = funcFiles(idx);
    funcDurations = funcDurations(idx);
    
    if length(etFiles)
        fprintf('Eyetracking files found : \n');
        for e = 1:length(etFiles)
            fprintf('%s %s %f\n', etFiles{e}, etTimeStamps(e), etDurations(e));
        end
        fprintf('\n');
    end
    
    idx_incomplete = find(etDurations < min(funcDurations));
    
    outdir = 'Unsorted/Incomplete';
    for i = idx_incomplete
        
        createdDirs{end+1} = outdir;
        
        fprintf('Moving %s to %s\n', etFiles{i}, outdir);
        if MOVEFLAG
            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end
            system(sprintf('mv "%s" "%s"', etFiles{i}, outdir));
        end
    end
    
    etDurations(idx_incomplete) = [];
    etTimeStamps(idx_incomplete) = [];
    etFiles(idx_incomplete) = [];
    
    %%% Check for duplicates
    %%% if duplicates found, move copies (with shorter duration or all but one) to unsorted/incomplete
    %%% 
    if length(etTimeStamps) ~= length(unique(etTimeStamps))
        
        idx_remove = [];
        for eTS = sort(unique(etTimeStamps))
           
            idx = find(ismember(etTimeStamps, eTS));
            if length(idx) > 1

                idx_shortDur = idx(find(etDurations(idx) < max(etDurations(idx))));
                
                %%% if all duplicates of same size
                if length(idx) == length(idx_shortDur)
                    idx_shortDur(end) = [];
                end
                
                idx_remove = [idx_remove idx_shortDur];
               
                outdir = 'Unsorted/Incomplete';
                for i = idx_shortDur

                    createdDirs{end+1} = outdir;

                    fprintf('Moving %s to %s\n', etFiles{i}, outdir);
                    if MOVEFLAG
                        if ~exist(outdir, 'dir')
                            mkdir(outdir);
                        end
                    system(sprintf('mv "%s" "%s"', etFiles{i}, outdir));
                    end
                end
            end
        end
        
        etFiles(idx_remove) = [];
        etTimeStamps(idx_remove) = [];
        etDurations(idx_remove) = [];
    end
    
    for e = 1:length(etFiles)

        inFile = etFiles{e};
        [~, eName, ~] = fileparts(inFile);
        outName = regexprep([eName '.txt'], '[^A-Za-z0-9-.]', '_');

        %%% If #Eyetracking files == #Functional files
        if length(etFiles) == length(funcFiles) && ...
            etDurations(e) >= funcDurations(e) && abs(etTimeStamps(e) - funcTimeStamps(e)) < 0.0104 %%% 15mins 
        
            outdir = 'Functional/EyeTracking';
            outPrefix = strtok(strrep(funcFiles{e}, 'Functional/', ''), '/');
            outFile = fullfile(outdir, sprintf('%s_%s', outPrefix, outName));

        %%% If #Eyetracking files == #Functional files but file is of
        %%% lesser duration than functional duration, move to unsorted
        elseif length(etFiles) == length(funcFiles)
            
            outdir = 'Unsorted/Unmatched';
            outFile = fullfile(outdir, outName);

        %%% If #Eyetracking files not equal to #Functional files, cannot
        %%% decide which eyetracking files are for which task, move to
        %%% Eyetracking folder
        elseif length(etFiles) ~= length(funcFiles)
            
            outdir = 'Functional/EyeTracking';
            outFile = fullfile(outdir, outName);
            
        end

        fprintf('Moving %s to %s\n', inFile, outFile);
        if MOVEFLAG
            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end
            system(sprintf('mv "%s" %s', inFile, outFile));
        end
        createdDirs{end+1} = outdir;
    end    

          
    
    
function moveUnsorted(unsortedInfo)

    global MOVEFLAG createdDirs
    
    fprintf('\n%s Moving unsorted data %s\n', repmat('-',1, 15), repmat('-',1, 15));
    
    for u = 1:length(unsortedInfo)
        
        scInfo = unsortedInfo{u};
        
        filename = scInfo.filename;
    
        outdir = scInfo.dest;
        
        indir = fileparts(filename);
        seriesName = getSeriesName(filename);
        outFolder = sprintf('%s_%03d_%s', scInfo.examno, scInfo.seriesno, seriesName);
        outdir = fullfile(outdir, outFolder);
        createdDirs{end+1} = outdir;
        
        %%% Dicom file
        if strcmp(scInfo.imagetype, 'dicom')
            fprintf('MOVING %s TO %s\n', indir, outdir);
            if MOVEFLAG
                if ~exist(outdir, 'dir')
                    mkdir(outdir);
                end
                system(sprintf('mv %s %s', indir, outdir));
            end
        
        %%% Spiral file
        elseif strcmp(scInfo.imagetype, 'spiral')
    
            infiles = fullfile(indir, ['*' scInfo.pfile(1:end-2) '*']);
            
            fprintf('MOVING %s TO %s\n', infiles, outdir);
            if MOVEFLAG
                if ~exist(outdir, 'dir')
                    mkdir(outdir);  
                end
                system(sprintf('mv %s %s', infiles, outdir));
            end
        end
    end
    
    
    
    
function update_sqlDB(sqlInfo)

    global taskNames
     
    %%% Connect to mysql database
    conn = mysqlConnect;

    columnNames = fieldnames(sqlInfo{1});
    for s = 1:length(sqlInfo)

        scanDates = cell2mat(mysqlQuery(conn, sprintf('SELECT distinct scanDate FROM imaging where subjectID=%d', sqlInfo{s}.subjectID)));
        timePoint = 1;
        if ~isempty(scanDates)
            if nnz(datenum(sqlInfo{s}.scanDate) - datenum(scanDates) > 180)
                timePoint = 2;
            end
        end
        
        %%% check if data already exists in mysql table
        rowCount = cell2mat(mysqlQuery(conn, sprintf('SELECT COUNT(*) FROM imaging where subjectID=%d and examID=%d and seriesNo=%d and scanDate=''%s''', ...
            sqlInfo{s}.subjectID, sqlInfo{s}.examID, sqlInfo{s}.seriesNo, sqlInfo{s}.scanDate)));
        
        if ~rowCount
            
            curUser = char(java.lang.System.getProperty('user.name'));
            
            queryStr = sprintf('INSERT INTO imaging (%s', sprintf('%s,', columnNames{:}));
            queryStr = sprintf('%s,userName, timePoint) VALUES (', queryStr(1:end-1));
            for col = 1:length(columnNames)
                fieldVal = eval(['sqlInfo{' num2str(s) '}.' columnNames{col}]);
                if strcmp(class(fieldVal), 'char')
                    queryStr = sprintf('%s''%s'',', queryStr, fieldVal);
                else
                    if mod(fieldVal, 1)
                        queryStr = sprintf('%s%.4f,', queryStr, fieldVal);
                    else
                        queryStr = sprintf('%s%d,', queryStr, fieldVal);
                    end
                end
            end

            queryStr = sprintf('%s,''%s'',%d);', queryStr(1:end-1), curUser, timePoint);
            fprintf('%s\n', queryStr);
            mysqlQuery(conn, queryStr);
        end

        %%% insert functional information into functional table
        if nnz(cellfun(@(x) ~isempty(strfind(sqlInfo{s}.seriesName, x)), [taskNames, 'Resting']))
            rowCount = cell2mat(mysqlQuery(conn, sprintf('SELECT COUNT(*) FROM functional where imageCentralLocation=''%s''', sqlInfo{s}.imagecentralLocation)));
            if ~rowCount
                queryStr = sprintf('insert into functional (imageCentralLocation, include, notes) values (''%s'', -2, '''');', sqlInfo{s}.imagecentralLocation);
                fprintf('%s\n', queryStr);
                mysqlQuery(conn, queryStr);
            end
            
            if nnz(cellfun(@(x) ~isempty(strfind(sqlInfo{s}.seriesName, x)), taskNames))
       
                results = mysqlQuery(conn, sprintf('select subjectID, sessionID, examID, seriesNo from imaging where imageCentralLocation=''%s''', sqlInfo{s}.imagecentralLocation));
                subjectID = results{1};
                sessionID = results{2};
                examID = results{3};
                seriesNo = results{4};
                
                results = mysqlQuery(conn, sprintf(['select imageCentralLocation, seriesNo from imaging where seriesName ', ...
                                                               '= ''InplaneAnatomy'' and subjectID=%d and sessionID=''%s'' and examID = %d'], subjectID, sessionID, examID));
                inplaneFile = results(1);
                inplaneSeriesNo = results{2};

                cur_inplaneFile = cell2mat(mysqlQuery(conn, sprintf('select inplane from functional where imageCentralLocation=''%s''', sqlInfo{s}.imagecentralLocation)));
                if length(inplaneFile)>1

                    diffNo = abs(inplaneSeriesNo - seriesNo);
                    idx = find(diffNo == min(diffNo));
                    new_inplaneFile = inplaneFile{idx(1)};

                elseif length(inplaneFile)
                   new_inplaneFile = inplaneFile{1};
                else
                   continue
                end

                %%% Do not overwrite existing inplane files
                if isempty(cur_inplaneFile)%%strcmp(new_inplaneFile, cur_inplaneFile{1})
                    queryStr = sprintf('UPDATE functional set inplane = ''%s'' where imageCentralLocation = ''%s''', new_inplaneFile, sqlInfo{s}.imagecentralLocation);
                    fprintf('%s\n', queryStr);
                    mysqlQuery(conn, queryStr);
                end
            end
        end

        %%% insert FSPGR information into anatomical table
        if strcmp(sqlInfo{s}.seriesName, 'FSPGR')
            rowCount = cell2mat(mysqlQuery(conn, sprintf('SELECT COUNT(*) FROM anatomical where imageCentralLocation=''%s''', sqlInfo{s}.imagecentralLocation)));
            if ~rowCount
                queryStr = sprintf('insert into anatomical (imageCentralLocation, wrapArtifactRating, motionArtifactRating, usability, notes) values (''%s'', -2, -2, -2, '''');', ...
                    sqlInfo{s}.imagecentralLocation);
                fprintf('%s\n', queryStr);
                mysqlQuery(conn, queryStr);
            end

        end
        
        %%% insert HARDI information into dti table
        if strcmp(sqlInfo{s}.seriesName, 'HARDI')
            rowCount = cell2mat(mysqlQuery(conn, sprintf('SELECT COUNT(*) FROM dti where imageCentralLocation=''%s''', sqlInfo{s}.imagecentralLocation)));
            if ~rowCount
                queryStr = sprintf('insert into dti (imageCentralLocation, include) values (''%s'', -2);', sqlInfo{s}.imagecentralLocation);
                fprintf('%s\n', queryStr);
                mysqlQuery(conn, queryStr);
            end

        end
    end
    
    
    
    
    
    
    
    
%     study_seriesNames('T2 Anatomy') = 'T2';
%     study_seriesNames('FSPGR 1x1x1 (w/o BRAVO)') = 'FSPGR';
%     study_seriesNames('HiRes_Hippocampus') = 'HiRes_Hippocampus';
%     study_seriesNames('Sprlio_HappyCafe') = 'HappyCafe';
%     study_seriesNames('Sprlio_HALF ER GoNoGo 1') = 'GoNoGo1';
%     study_seriesNames('Sprlio_HALF ER GoNoGo 2') = 'GoNoGo2';
%     study_seriesNames('Sprlio_Rest') = 'Resting';
%     study_seriesNames('Sprlio-Resting') = 'Resting';
%     study_seriesNames('Sprlio-MID') = 'MK_MID';
%     study_seriesNames('Sprlio4skip.5_HappyCafe') = 'HappyCafes';
%     study_seriesNames('T2Anatomy_fMRI_4skip1') = 'T2';


% idx = find(cellfun(@(x) strcmp(x.seriesName, 'InplaneAnatomy'), sqlInfo));

% inplanes_examIDs = arrayfun(@(i) sqlInfo{idx(i)}.examID, 1:length(idx));
% inplanes_seriesNos = arrayfun(@(i) sqlInfo{idx(i)}.seriesNo, 1:length(idx));
% inplanes_imgCentralLocations = arrayfun(@(i) sqlInfo{idx(i)}.imagecentralLocation, 1:length(idx), 'UniformOutput', false);

%             idx = find(ismember(inplanes_examIDs, sqlInfo{row}.examID));
%             cur_seriesNos = inplanes_seriesNos(idx);
%             cur_imgCentralLocations = inplanes_imgCentralLocations(idx);
%             if length(idx)==1
%                 cur_inplane = cur_imgCentralLocations{1};
%             elseif length(idx) > 1
%                 cur_inplane = cur_imgCentralLocations{max(find(sqlInfo{row}.seriesNo - cur_seriesNos > 0))};
%             else
%                 cur_inplane = '';
%             end
%             if rowcount
%                 query_str = sprintf('update functional set subjectID=%d, seriesName=''%s'', inplane=''%s'', include=1, notes='''' where imageCentralLocation=''%s'';', ...
%                     sqlInfo{row}.subjectID, sqlInfo{row}.seriesName, cur_inplane, sqlInfo{row}.imagecentralLocation);
%             else
%                 query_str = sprintf('insert into functional (subjectID, seriesName, imageCentralLocation, inplane, include, notes, timePoint) values (%d, ''%s'', ''%s'', ''%s'', 1, '''', %d);', ...
%                     sqlInfo{row}.subjectID, sqlInfo{row}.seriesName, sqlInfo{row}.imagecentralLocation, cur_inplane, timepoint);
%             end

% %%% insert Resting information into resting table
% if strcmp(sqlInfo{row}.seriesName, 'Resting')
%     query_str = sprintf('SELECT COUNT(*) FROM resting where imageCentralLocation=''%s''', sqlInfo{row}.imagecentralLocation);
%     rowcount = mysql(query_str);
% 
%     if ~rowcount
%         query_str = sprintf('insert into resting (imageCentralLocation, usability, notes) values (''%s'', 0, '''');', ...
%             sqlInfo{row}.imagecentralLocation);
%         mysql(query_str);
%     end
% 
% end



%         if rowcount
% 
%             query_str = sprintf('UPDATE imaging SET ');
%             for col = 1:length(columnNames)
%                 fieldval = eval(['sqlInfo{' num2str(row) '}.' columnNames{col}]);
%                 if strcmp(class(fieldval), 'char')
%                     query_str = sprintf('%s%s=''%s'',', query_str, columnNames{col}, fieldval);
%                 else
%                     if mod(fieldval, 1)
%                         query_str = sprintf('%s%s=%.4f,', query_str, columnNames{col}, fieldval);
%                     else
%                         query_str = sprintf('%s%s=%d,', query_str, columnNames{col}, fieldval);
%                     end
%                 end
%             end
%             query_str = sprintf('%s where subjectID=%d and examID=%d and seriesNo=%d and scanDate=''%s'';', query_str(1:end-1), ...
%                 sqlInfo{row}.subjectID, sqlInfo{row}.examID, ...
%                 sqlInfo{row}.seriesNo, sqlInfo{row}.scanDate);
%         else

