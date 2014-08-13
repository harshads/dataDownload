function Turner_copyFMPcsv(subjectID)

    %%% Local drives where data is saved
    dataDir = '/Volumes/Projects/TURNER/Data/sortedData/';

    if ~exist('subjectID', 'var')
        inDir = uigetdir(dataDir, 'Select subject directory');
    else
        inDir = fullfile(dataDir, subjectID);

        if ~exist(inDir, 'dir')
            dataDir = fullfile('/Volumes/ImageCentral/', ['Subjs_' num2str(floor(str2num(subjectID)/100)*100)]);
            inDir = fullfile(dataDir, subjectID);
        end
    end
    
    cd(dataDir);

    fprintf('Reading FileMakerPro csv file from %s', inDir);
    [~, inFiles] = system(['find ' inDir ' -maxdepth 2 -type f -name "mergedCibsrDcm2fmp.csv"']);

    inFiles = strsplit(inFiles, char(10))';

    for f = 1:length(inFiles)

        fprintf('\nReading %s\n\n', inFiles{f});
        fid = fopen(inFiles{f}, 'r');
        headerLine = fgetl(fid);
        line = {};
        while ~feof(fid)
            line{end+1} = fgetl(fid);
        end
        fclose(fid);

        for l = 1:length(line)

            curLine = line{l};
            curValues = strsplit(curLine, ',');
            imgLocation = curValues{12};

            if ~exist(imgLocation, 'dir')
                fprintf('Cannot find the path %s\n', imgLocation);
                continue;
            end

            if isempty(strfind(imgLocation, 'Functional'))
                outFile = fullfile(imgLocation, 'cibsrDcm2fmp.csv');
            else
                outFile = fullfile(fileparts(imgLocation), [curValues{11} '_cibsrDcm2fmp.csv']);
            end

            fprintf('Writing output file %s\n', outFile);
            fid = fopen(outFile, 'w');
            fprintf(fid, '%s\n%s\n', headerLine, curLine);
            fclose(fid);

        end
    end
