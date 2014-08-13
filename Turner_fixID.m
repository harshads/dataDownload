function Turner_fixID(oldID, newID)

conn = mysqlConnect;

%%% Update dataImaging table
mysqlQuery(conn, sprintf('UPDATE dataImaging SET subjectID = %s WHERE subjectID = %s', newID, oldID));


%%% Update functional table
mysqlQuery(conn, sprintf('UPDATE functional SET subjectID=%s WHERE subjectID=%s', newID, oldID));
mysqlQuery(conn, sprintf('UPDATE functional SET imageCentralLocation = REPLACE(imageCentralLocation, %s, %s) WHERE imageCentralLocation LIKE ''%s%%''', oldID, newID, oldID));
results = unique(mysql(sprintf('select inplane from functional where inplane like ''%s/%%''', oldID)));
for r = 1:length(results)
    new_inplane = sprintf('%s%s', newID, results{r}(length(oldID)+1:end));
    mysqlQuery(conn, sprintf('UPDATE functional SET inplane = ''%s'' WHERE inplane = ''%s''', new_inplane, results{r}));
end

%%% Update anatomical table
results = unique(mysql(sprintf('select imageCentralLocation from anatomical where imageCentralLocation like ''%s/%%''', oldID)));
for r = 1:length(results)
    new_imageCentralLocation = sprintf('%s%s', newID, results{r}(length(oldID)+1:end));
    mysqlQuery(conn, sprintf('UPDATE anatomical SET imageCentralLocation = ''%s'' WHERE imageCentralLocation = ''%s''', new_imageCentralLocation, results{r}));
end

%%% Update dti table
mysqlQuery(conn, sprintf('UPDATE dti SET imageCentralLocation = REPLACE(imageCentralLocation, %s, %s) WHERE imageCentralLocation LIKE ''%s%%''', oldID, newID, oldID));

%%% Update functionalQuality table
mysqlQuery(conn, sprintf('UPDATE functionalQuality SET subjectID=%s WHERE subjectID=%s', newID, oldID));

%%% Update imaging table
mysqlQuery(conn, sprintf('UPDATE imaging SET subjectID=%s WHERE subjectID=%s', newID, oldID));
results = unique(mysql(sprintf('select imageCentralLocation from imaging where imageCentralLocation like ''%s/%%''', oldID)));
for r = 1:length(results)
    new_imageCentralLocation = sprintf('%s%s', newID, results{r}(length(oldID)+1:end));
    mysqlQuery(conn, sprintf('UPDATE imaging SET imageCentralLocation = ''%s'' WHERE imageCentralLocation = ''%s''', new_imageCentralLocation, results{r}));
end

%%% Update resting table
results = unique(mysql(sprintf('select imageCentralLocation from resting where imageCentralLocation like ''%s/%%''', oldID)));
for r = 1:length(results)
    new_imageCentralLocation = sprintf('%s%s', newID, results{r}(length(oldID)+1:end));
    mysqlQuery(conn, sprintf('UPDATE resting SET imageCentralLocation = ''%s'' WHERE imageCentralLocation = ''%s''', new_imageCentralLocation, results{r}));
end

%%% Update tracking table
mysqlQuery(conn, sprintf('UPDATE tracking SET subjectID=%s WHERE subjectID=%s', newID, oldID));


%%% Modify nifti folder
dirs = {'/Volumes/Projects/TURNER/Data/nifti', '/Volumes/Projects/TURNER/Data/sortedData'};

curUser = char(java.lang.System.getProperty('user.name'));
if strcmp(curUser, 'harshads')
    dirs{end+1} = '/Users/harshads/Sites/Turner/QA';
end
    
for d = 1:length(dirs)

    olddir = fullfile(dirs{d}, oldID);
    if exist(olddir, 'dir')
        newdir = fullfile(dirs{d}, newID);

        if exist(newdir, 'dir')
            system(sprintf('mv %s/* %s', olddir, newdir));
        else
            system(sprintf('mv %s %s', olddir, newdir));
        end
    end
end

%%% Move data to correct folder in ImageCentral
%%% Check if all data has been moved
datadir = '/Volumes/ImageCentral';
olddir = fullfile(datadir, ['Subjs_' num2str(floor(str2num(oldID)/100)*100)] , oldID);
if exist(olddir, 'dir')
    newdir = fullfile(datadir, ['Subjs_' num2str(floor(str2num(newID)/100)*100)] , newID);

    if exist(newdir, 'dir')
        sessdirs = dir(fullfile(olddir, '*'));
        sessdirs = sort({sessdirs(3:end).name})';
        for s = 1:length(sessdirs)
            if isdir(fullfile(olddir, sessdirs{s})) && ~exist(fullfile(newdir, sessdirs{s}), 'dir')
                system(sprintf('mv %s %s', fullfile(olddir, sessdirs{s}), newdir));
                
                
                %%% Update mergedCibsrDcm2fmp.csv
                cur_csvfile = fullfile(newdir, sessdirs{s}, 'mergedCibsrDcm2fmp.csv');
                
                fid = fopen(cur_csvfile, 'r');
                headerLine = fgetl(fid);
                line = {};
                while ~feof(fid)
                    line{end+1} = fgetl(fid);
                end
                fclose(fid);

                headerNames = strsplit(headerLine, ',');
                for l = 1:length(line)

                    cur_line = line{l};
                    cur_values = strsplit(cur_line, ',');
                    icLocation = cur_values{strmatch('ImageCentral_folder_name', headerNames)};

                    %%% Update the csv "notes" field
                    cur_values{strmatch('Subj_Num', headerNames)} = newID;
                    
                    idx = strmatch('ImageCentral_folder_name', headerNames);
                    oldval = cur_values{idx};
                    cur_values{idx} = [newID oldval(length(oldID)+1:end)];

                    new_line = cur_values{1};
                    for n = 2:length(cur_values)
                        new_line = sprintf('%s,%s', new_line, cur_values{n});
                    end

                    %%% Update the line variable
                    line{l} = new_line;
                end

                fid = fopen(cur_csvfile, 'w');
                fprintf(fid, headerLine);
                fprintf(fid, '\n');  
                for l = 1:length(line)
                    fprintf(fid, line{l});
                    fprintf(fid, '\n');
                end
                fclose(fid);
            end
        end
        
        %%% Delete old folder if its empty
        system(sprintf('find %s -type f -name ''.DS_Store'' -exec rm {} \\; rmdir %s;', olddir, olddir));
        if exist(olddir, 'dir')
           fprintf('\n\nCOULD NOT MOVE ALL FOLDERS FROM %s\n\n', olddir);
           return;
        end
        
        %%% Update individual Dcm2fmp.csv's
        Turner_copyFMPcsv(newID);
    else
        system(sprintf('mv %s %s', olddir, newdir));
    end
end


%%% Fix mrVista sessions
datadir = '/Volumes/Projects/TURNER/Data/mrVista';

olddir = fullfile(datadir, oldID);
if exist(olddir, 'dir')
    
    newdir = fullfile(datadir, newID);
    if exist(newdir, 'dir')
        sessdirs = dir(fullfile(olddir, '*'));
        sessdirs = sort({sessdirs(3:end).name})';
        for s = 1:length(sessdirs)
            if isdir(fullfile(olddir, sessdirs{s})) && ~exist(fullfile(newdir, sessdirs{s}), 'dir')
                system(sprintf('mv %s %s', fullfile(olddir, sessdirs{s}), newdir));
            end
        end
        
        %%% Delete old folder if its empty
        system(sprintf('find %s -type f -name ''.DS_Store'' -exec rm {} \\; rmdir %s;', olddir, olddir));
        if exist(olddir, 'dir')
           fprintf('\n\nCOULD NOT MOVE ALL FOLDERS FROM %s\n\n', olddir);
           return;
        end
        
    else
        system(sprintf('mv %s %s', olddir, newdir));
    end
    
    [~, matfiles] = system(sprintf('find %s -type f -name "mrSESSION.mat"', newdir));
    matfiles = strsplit(matfiles, char(10));
    for m = 1:length(matfiles)
        cd(fileparts(matfiles{m}));
        mrVista_updateSessionPath;
    end
end


