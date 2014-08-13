%%%% update mysql "tracking" table using information from general tracking excel file
%%%% update mysql "functionalQuality" table using information from general tracking excel file
% function Turner_updateSQL()

function Turner_updateSQL

    conn = mysqlConnect;

    dataDir1 = '/Volumes/Bryce/Projects/Turner/Data/';
    dataDir2 = '/Volumes/Projects/TURNER/Data/sortedData';

    %%% Find the last modified file (if multiple files exist with the name
    %%% GeneralTracking
    xlsDir = '/Volumes/Projects/TURNER/K23Mosaic/Tracking';
    files = dir(fullfile(xlsDir, 'GeneralTracking*'));
    [~, idxMatched] = sort(arrayfun(@(x) datenum(files(x).date, 'dd-mmm-yyyy HH:MM:SS'), 1:length(files)));
    xlsFile = fullfile(xlsDir, files(idxMatched(end)).name);

    %%% Read the General Tracking excel file
    [~, ~, xlsInfo] = xlsread(xlsFile);

    %%% Find the row where the ID's start (search for token "ID")
    idxStart = find(cellfun(@(x) strcmp('ID', x), xlsInfo(:,1)))+1;

    %%% Get the header names
    headerNames = xlsInfo(idxStart - 1, :);

    %%% Find columns with header names
    idxCol = find(cellfun(@isstr, headerNames));
    headerNames = headerNames(idxCol);

    %%% Exclude everything above
    xlsInfo = xlsInfo(idxStart:end, idxCol);

    %%% Find subjectIDs from the IMAGING data table
    mysql_subjectIDs1 = cell2mat(mysqlQuery(conn, 'select distinct(subjectID) from imaging order by subjectID'));

    %%% Find subjectIDs from the TRACKING data table
    mysql_subjectIDs2 = cell2mat(mysqlQuery(conn, 'select subjectID from tracking order by subjectID'));

    %%% Check to make sure no duplicate IDs in tracking
    if length(unique(mysql_subjectIDs2)) ~= length(mysql_subjectIDs2)
        fprintf('ERROR: Duplicate IDs in mysql "tracking" table\n');
        return;
    end

    %%% Runs either all the IDs again and updates; or runs only new IDs
    runAll = 1; 

    %%% If not RunAll, run for newer IDs
    if runAll
        new_subjectIDs = mysql_subjectIDs1; 
    else
        new_subjectIDs = setdiff(mysql_subjectIDs1, mysql_subjectIDs2);
    end

    %%% Find the subjectIDs row index (check that the cell contains non empty
    %%% numeric string
    idxSub = find(cellfun(@(x) isnumeric(x) && ~isnan(x), xlsInfo(:,1)));%%find(cellfun(@(x) (~isstr(x) && ~isnan(x)), info(:,1)));

    fprintf('\n%s\nUpdating DATA IMAGING table\n%s\n', repmat('*', 1, 40), repmat('*', 1, 40));

    for n = 1:length(new_subjectIDs)

        results = mysqlQuery(conn, sprintf('select distinct timepoint, sessionID from imaging where subjectID = %d', new_subjectIDs(n)));

        timePoints = cell2mat(results(:, 1));
        sessionIDs = results(:, 2);

        for t = 1:length(timePoints)

            hasVisionData = 0;
            cnt = cell2mat(mysqlQuery(conn, sprintf('select count(*) from functional where imageCentralLocation not like ''%%_Resting'' and imageCentralLocation like ''%%%s%%''', sessionIDs{t})));
            if cnt > 0
                hasVisionData = 1;
            end

            noRows = cell2mat(mysqlQuery(conn, sprintf('select hasVisionData from dataImaging where subjectID = %d and timepoint = %d and sessionID = ''%s''', new_subjectIDs(n), timePoints(t), sessionIDs{t})));

            %%% If data EXISTS, and has the SAME VALUE
            if ~isempty(noRows) && noRows == hasVisionData
                continue;
            else
                fprintf('\nUpdating IMAGING table for SubjectID %d \n', new_subjectIDs(n));

                %%% If data DOES NOT EXIST in table
                if isempty(noRows)
                    queryStr = sprintf('INSERT INTO dataImaging (subjectID, timepoint, sessionID, hasVisionData) values (%d, %d, ''%s'', %d)', ...
                        new_subjectIDs(n), timePoints(t), sessionIDs{t}, hasVisionData);
                %%% If data has DIFFERENT VALUE in table
                else
                    queryStr = sprintf('UPDATE dataImaging set hasVisionData = %d where subjectID = %d and timepoint = %d and sessionID = ''%s''', ...
                            hasVisionData, new_subjectIDs(n), timePoints(t), sessionIDs{t});
                end
                fprintf('%s\n', queryStr);
                mysqlQuery(conn, queryStr);
            end            
        end


        %%% Find the first dicom file to read header
        iclocDir = mysqlQuery(conn, sprintf('select imageCentralLocation from imaging where subjectID = %d', new_subjectIDs(n)));
        cnt = 1;
        dcmInfo = [];
        while cnt <= length(iclocDir)
            inFile = fullfile(dataDir1, iclocDir{cnt}, 'I0001.dcm');
            if ~exist(inFile, 'file')
                inFile = fullfile(dataDir2, iclocDir{cnt}, 'I0001.dcm');
            end
            if exist(inFile, 'file')
                dcmInfo = dicominfo(inFile);
                break;
            end
            cnt = cnt + 1;
        end

        if isempty(dcmInfo)
            continue;
        end

        %%% Check if subject exists in tracing folder
        %%% Then Check the name from the dicom with the name in excel
        %%% Needed sometimes when CIBSR IDs have not been entered
        noRows = cell2mat(mysqlQuery(conn, sprintf('select count(*) from tracking where subjectID = %d', new_subjectIDs(n))));
        idxMatched = min(strmatch(new_subjectIDs(n), cell2mat(xlsInfo(idxSub,1))));
        if isempty(idxMatched)
            fname = dcmInfo.PatientName.GivenName;
            lname = dcmInfo.PatientName.FamilyName;

            %%% Check the First Name and Last Name in the dicom header with the
            %%% excel file header
            idxMatched = find(cellfun(@(x) strcmp(x, lname), xlsInfo(:, strmatch('Last Name', headerNames))));
            idxMatched = idxMatched(ismember(xlsInfo(idxMatched, strmatch('First Name', headerNames)), fname));
        else
            idxMatched = idxSub(idxMatched);
        end

        %%% If a match has been found
        if ~isempty(idxMatched)

            groupName = xlsInfo{idxMatched, strmatch('Group', headerNames)};
            subAge = xlsInfo{idxMatched, strmatch('Age at Participation', headerNames)};
            subGender = xlsInfo{idxMatched, strmatch('Sex', headerNames)};
            if ~strcmp(subGender, 'F') && ~strcmp(subGender, 'M')
                subGender = '';
            end
            studyK23 = strcmp(xlsInfo{idxMatched, strmatch('K23', headerNames)}, 'x');
            studyGBB2 = strcmp(xlsInfo{idxMatched, strmatch('GBB2', headerNames)}, 'x');
            studyKSTS = strcmp(xlsInfo{idxMatched, strmatch('KS/TS', headerNames)}, 'x');

            if isnan(subAge)
                subAge = datevec(datenum(dcmInfo.StudyDate, 'yyyymmdd') - datenum(dcmInfo.PatientBirthDate, 'yyyyddmm'));
                subAge = subAge(1);
            end

            %%% Insert new information into TRACKING
            if ~noRows
                queryStr = sprintf('INSERT INTO tracking (subjectID,groupName,age,gender,study_K23,study_GBB2,study_KSTS) values (%d,''%s'',%d,''%s'', %d, %d, %d)', ...
                    new_subjectIDs(n), groupName, subAge, subGender, studyK23, studyGBB2, studyKSTS);
                fprintf('%s\n', queryStr);
                mysqlQuery(conn, queryStr);

            %%% Check if information has changed in TRACKING
            else

                results = mysqlQuery(conn, sprintf('select groupName,age,gender,study_K23,study_GBB2, study_KSTS from tracking where subjectID = %d', new_subjectIDs(n)));

                str1 = sprintf('%s,%d,%s,%d,%d,%d', results{1}, results{2}, results{3}, results{4}, results{5}, results{6});
                str2 = sprintf('%s,%d,%s,%d,%d,%d', groupName, subAge, subGender, studyK23, studyGBB2, studyKSTS);

                if ~strcmpi(str1, str2)
                    queryStr = sprintf('UPDATE tracking set groupName = ''%s'', age = %d, gender = ''%s'', study_K23 = %d, study_GBB2 = %d, study_KSTS = %d where subjectID = %d', ...
                        groupName, subAge, subGender, studyK23, studyGBB2, studyKSTS, new_subjectIDs(n));
                    fprintf('%s\n', queryStr);
                    mysqlQuery(conn, queryStr);
                end
            end            

        %%% If not matches were found    
        else
            if ~noRows
                queryStr = sprintf('INSERT INTO tracking (subjectID,age,study_K23,study_GBB2,study_KSTS) values (%d, 0, 0, 0, 0)', new_subjectIDs(n));
                fprintf('%s\n', queryStr);
                mysqlQuery(conn, queryStr);
            else
                %% do not update sql table; values might have been updated manually
                %mysql(sprintf('UPDATE tracking set age = 0, study_K23 = 0, study_GBB2 = 0, study_KSTS = 0, hasImagingData = %d, hasVisionData = %d where subjectID = %d', hasImagingData, hasVisionData, new_subjectIDs(n)));
            end
        end
    end


    %%% Updates FUNCTIONAL QUALITY with new IDs found in the DATA IMAGING table
    %%% which have VISION data
    clear;

    conn = mysqlConnect;

    fprintf('\n%s\nUpdating FUNCTIONAL QUALITY table\n%s\n', repmat('*', 1, 40), repmat('*', 1, 40));

    %%% updates functionalQuality table using IDs from tracking subjectIDs
    %%% who have vision data
    subIDs1 = mysqlQuery(conn, 'select concat(subjectID, ",", timepoint) from dataImaging where hasVisionData = 1');
    subIDs2 = mysqlQuery(conn, 'select concat(subjectID, ",", timepoint) from functionalQuality');

    new_subjectIDs = setdiff(subIDs1, subIDs2);

    for n = 1:length(new_subjectIDs)
        idxMatched = strfind(new_subjectIDs{n}, ',');
        queryStr = sprintf('insert into functionalQuality (subjectID, timepoint) values (%s, %s)', new_subjectIDs{n}(1:idxMatched-1), new_subjectIDs{n}(idxMatched+1:end));

        fprintf('%s\n', queryStr);
        mysqlQuery(conn, queryStr);
    end



    %%% Updates IMAGING table (sets the study information using TRACKING table
    %%% information)
    clear;

    conn = mysqlConnect;

    %%% updates study column in imaging table using IDs from tracking subjectIDs
    %%% who have vision data
    subIDs = cell2mat(mysqlQuery(conn, 'select distinct subjectID from imaging where study = '''''));
    studyNames = {'K23Mosaic', 'GBB2', 'KS/TS'};

    fprintf('\n%s\nUpdating IMAGING table\n', repmat('*', 1, 40));
    fprintf('Setting study information for subjects using TRACKING table information\n%s\n\n', repmat('*', 1, 40));

    for n = 1:length(subIDs)


        study = cell2mat(mysqlQuery(conn, sprintf('select study_K23, study_GBB2, study_KSTS from tracking where subjectID=%d', subIDs(n))));
        if ~any(subIDs)
            continue;
        end

        str = sprintf('%s;', studyNames{find(study)});
        queryStr = sprintf('update imaging set study=''%s'' where subjectID = %d;', str(1:end-1), subIDs(n));
        fprintf('%s\n', queryStr);
        mysqlQuery(conn, queryStr);
    end





























% % % for n = 1:length(new_subjectIDs)
% % %     
% % %     fprintf('\nUpdating SubjectID %d \n', new_subjectIDs(n));
% % %     
% % %     hasImagingData = 0;
% % %     cnt = mysql(sprintf('select count(subjectID) from imaging where subjectID = %d and seriesName NOT IN (''ASSET'', ''Localizer'')', new_subjectIDs(n)));
% % %     if cnt ~= 0
% % %         hasImagingData = 1;
% % %     end
% % % 
% % %     hasVisionData = 0;
% % %     cnt = mysql(sprintf('select count(subjectID) from functional where subjectID = %d and seriesName NOT IN (''ASSET'', ''Localizer'')', new_subjectIDs(n)));
% % %     if cnt ~= 0
% % %         hasVisionData = 1;
% % %     end
% % %     
% % %     noRows = mysql(sprintf('select count(subjectID) from tracking where subjectID = %d', new_subjectIDs(n)));
% % %      
% % %     indir = mysql(sprintf('select imageCentralLocation from imaging where subjectID = %d', new_subjectIDs(n)));
% % %     cnt = 1;
% % %     dcm_info = [];
% % %     while cnt <= length(indir)
% % %         infile = fullfile('/Volumes/Bryce/Projects/Turner/Data/', indir{cnt}, 'I0001.dcm');
% % %         if ~exist(infile, 'file')
% % %             infile = fullfile('/Users/harshads/Projects/Turner/Data/', indir{cnt}, 'I0001.dcm');
% % %         end
% % %         if exist(infile, 'file')
% % %             dcm_info = dicominfo(infile);
% % %             break;
% % %         end
% % %         cnt = cnt + 1;
% % %     end
% % %     
% % %     idx = min(find(cell2mat(info(xl_idx,1)) - new_subjectIDs(n) == 0));
% % %     if isempty(idx)
% % %         if isempty(dcm_info)
% % %             continue;
% % %         end
% % %         fname = dcm_info.PatientName.GivenName;
% % %         lname = dcm_info.PatientName.FamilyName;
% % % 
% % %         idx = find(cellfun(@(x) strcmp(x, lname), info(:,3)));
% % %         idx = idx(ismember(info(idx, 4), fname));
% % %     else
% % %         idx = xl_idx(idx);
% % %     end
% % % 
% % %     if ~isempty(idx)
% % %         groupName = info{idx,7};
% % %         age = info{idx, 5};
% % %         gender = info{idx, 6};
% % %         study_K23 = strcmp(info{idx, 8}, 'x');
% % %         study_GBB2 = strcmp(info{idx, 9}, 'x');
% % %         study_KSTS = strcmp(info{idx, 10}, 'x');
% % %         
% % %         if isnan(age)
% % %             age = datevec(datenum(dcm_info.StudyDate, 'yyyymmdd') - datenum(dcm_info.PatientBirthDate, 'yyyyddmm'));
% % %             age = age(1);
% % %         end
% % %         
% % %         if ~noRows
% % %             mysql(sprintf('INSERT INTO tracking (subjectID,groupName,age,gender,study_K23,study_GBB2,study_KSTS,hasImagingData,hasVisionData) values (%d,''%s'',%d,''%s'', %d, %d, %d, %d, %d)', ...
% % %                 new_subjectIDs(n), groupName, age, gender, study_K23, study_GBB2, study_KSTS, hasImagingData, hasVisionData));
% % %         else
% % %             mysql(sprintf('UPDATE tracking set groupName = ''%s'', age = %d, gender = ''%s'', study_K23 = %d, study_GBB2 = %d, study_KSTS = %d, hasImagingData = %d, hasVisionData = %d where subjectID = %d', ...
% % %                 groupName, age, gender, study_K23, study_GBB2, study_KSTS, hasImagingData, hasVisionData, new_subjectIDs(n)));
% % %         end            
% % %     else
% % %         if ~noRows
% % %             mysql(sprintf('INSERT INTO tracking (subjectID,age,study_K23,study_GBB2,study_KSTS,hasImagingData,hasVisionData) values (%d, 0, 0, 0, 0, %d, %d)', new_subjectIDs(n), hasImagingData, hasVisionData));
% % %         else
% % %             %% do not update sql table; values might have been updated manually
% % %             %mysql(sprintf('UPDATE tracking set age = 0, study_K23 = 0, study_GBB2 = 0, study_KSTS = 0, hasImagingData = %d, hasVisionData = %d where subjectID = %d', hasImagingData, hasVisionData, new_subjectIDs(n)));
% % %         end
% % %     end
% % % 
% % % end





%     if ~study_K23 && ~study_GBB2 && ~study_KSTS
%         continue;
%     end


% % for f = 1:length(info)
% %     
% %     subjectID = info{f,1};
% %     if (isstr(subjectID) && isempty(str2num(subjectID))) || isnan(subjectID) || any(cellfun(@(f) isnan(f), info(f, [5 6 8:10]))) || ~ismember(info{f, 7}, {'Control', 'Monosomic', 'Mosaic'})
% %         continue;
% %     end
% %     
% %     noRows = mysql(['select count(subjectID) from tracking where subjectID = ' num2str(subjectID)]);
% %     
% %     if ~noRows && ismember(subjectID, mysql_subjectIDs)
% %         
% %         groupName = info{f,7};
% %         age = info{f, 5};
% %         gender = info{f, 6};
% %         study_K23 = strcmp(info{f, 8}, 'x');
% %         study_GBB2 = strcmp(info{f, 9}, 'x');
% %         study_KSTS = strcmp(info{f, 10}, 'x');
% %         
% %         if ~study_K23 && ~study_GBB2 && ~study_KSTS
% %             continue;
% %         end
% %         
% %         hasImagingData = 0;
% % 
% %         cnt = mysql(sprintf('select count(subjectID) from imaging where subjectID = %d and seriesName NOT IN (''ASSET'', ''Localizer'')', subjectID));
% %         if cnt ~= 0
% %             hasImagingData = 1;
% %         end
% %         
% %         hasVisionData = 0;
% % 
% %         cnt = mysql(sprintf('select count(subjectID) from functional where subjectID = %d and seriesName NOT IN (''ASSET'', ''Localizer'')', subjectID));
% %         if cnt ~= 0
% %             hasVisionData = 1;
% %         end
% %         
% %         mysql(sprintf('INSERT INTO tracking (subjectID,groupName,age,gender,study_K23,study_GBB2,study_KSTS,hasImagingData,hasVisionData) values (%d,''%s'',%d,''%s'', %d, %d, %d, %d, %d)', ...
% %             subjectID, groupName, age, gender, study_K23, study_GBB2, study_KSTS, hasImagingData, hasVisionData));
% %     end
% %     
% %     
% % end
% % 
% % %%% updates functionalQuality table using IDs from tracking subjectIDs
% % %%% who have vision data
% % tracking_subjectIDs = mysql('select subjectID from tracking where hasVisionData = 1');
% % 
% % functionalQ_subjectIDs = mysql('select subjectID from functionalQuality');
% % 
% % new_subjectIDs = setdiff(tracking_subjectIDs, functionalQ_subjectIDs);
% % 
% % if ~isempty(new_subjectIDs)
% %     arrayfun(@(x) mysql(sprintf('insert into functionalQuality (subjectID) values (%d)', x)), new_subjectIDs, 'UniformOutput', false);
% % end








%     cnt = mysql(sprintf('select count(subjectID) from functional where subjectID = %d', subjectID));
%     hasVisionData = 0;
%     if cnt ~= 0
%         hasVisionData = 1;
%     end
%     
%     mysql(['update tracking set hasVisionData = ' num2str(hasVisionData) ' where subjectID = ' num2str(subjectID)]);

