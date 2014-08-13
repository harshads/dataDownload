function data = mysqlQuery(conn, query_str)

    curs = exec(conn, query_str);
    results = fetch(curs);
    
    data = {};
    if ~strcmpi(results.Data, 'No Data')
        data = results.Data;
    end