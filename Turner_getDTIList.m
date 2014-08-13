%%% Runs as cronjob every night
conn = mysqlConnect;

icLoc = mysqlQuery(conn, 'select imageCentralLocation from dti');

outfile = '/Volumes/Projects/TURNER/Data/Files/list_DTIdata.txt';
fid = fopen(outfile, 'w');
for f = 1:length(icLoc)
    examID = cell2mat(mysqlQuery(conn, sprintf('select examID from imaging where imageCentralLocation=''%s''', icLoc{f})));
    imgloc = strrep(icLoc{f}, 'DTI', ['E' num2str(examID)]);
    fprintf(fid, '%s\n', imgloc);
end
fclose(fid);

system(['cp ' outfile ' /Volumes/Bryce/Projects/Turner/Files/' ]);
