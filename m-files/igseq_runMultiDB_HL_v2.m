function table = igseq_runMultiDB_HL_v2( db )

    table_tmp = {};
    
    conn = igseq_db(db);
    
    %tmp = cellstr('Database');
    db_name = regexprep(db, '.*/', '');
    table_tmp = cellstr(db_name);
    
    q = 'select lineageID, count(*) as n from sequence group by lineageID order by n desc';
    lin_sizes = get(fetch(exec(conn,q)),'data');
    
    q1 = 'select pairedChain, count(distinct(chainType)) as n from sequence group by pairedChain having n > 1';
    s = get(fetch(exec(conn, q1)),'data');
    
    lin_size = {};
    for i = 1:length(s)
        q2 = strcat('select lineageID from sequence where pairedChain = ''', s(i,1), '''');
        t = get(fetch(exec(conn, q2{1})),'data');
        tu = unique( [t{:,1}] );
        
        for j = 1:length(tu)
            lin_sizes(cell2mat(lin_sizes(:,1)) == tu(1,1),:);
            tu(2,j) = cell2mat(get(fetch(exec(conn,q3)),'data'));
            lin_size = [lin_size tu'];
        end
            
    end
    lin_size
    
    table = vertcat(header, table);
    cell2csv('HL_db_stats.csv', table, ',', '2016', '.');
end
    
    