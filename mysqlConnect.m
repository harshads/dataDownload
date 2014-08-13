function conn = mysqlConnect

    %%% Current user
    curUser = char(java.lang.System.getProperty('user.name'));
    
    %%% serverLocal=1 for user harshads, else 0
    serverLocal = strcmp(curUser, 'harshads');
    
    if serverLocal
        host = 'localhost';
        curUser = 'root';
        pwd = 'c!85rsqL';
    else
        host = '171.65.63.28';
        pwd = 'sqluser';
    end

    conn = database('turner', curUser, pwd, '/Library/Java/Extensions/mysql-connector-java-5.1.27-bin.jar', ['jdbc:mysql://' host '/turner']);
    
    
    
    
    
    
% %     %%% Current user
% %     curUser = char(java.lang.System.getProperty('user.name'));
% %     
% %     %%% serverLocal=1 for user harshads, else 0
% %     serverLocal = strcmp(curUser, 'harshads');
% %     
% %     if serverLocal
% %         addpath ~/matlab/mysql/
% %     else
% %         addpath /Volumes/Projects/TURNER/Data/Files/mysql/
% %     end
% %     
% %     try
% %         mysql('use turner');
% %     catch err
% %         if serverLocal
% %             file = '/Volumes/Projects/TURNER/Data/Files/rootpwd';
% %             host = 'localhost';
% %             user = 'root';
% %         else
% %             file = '/Volumes/Projects/TURNER/Data/Files/sqluserpwd';
% %             host = '171.65.63.28';
% %             user = 'sqluser';
% %         end
% %         
% %         [~, pwd] = system(['cat ' file]);
% %         mysql('open', host, user, strtrim(pwd));
% %         mysql('use turner');
% %     end