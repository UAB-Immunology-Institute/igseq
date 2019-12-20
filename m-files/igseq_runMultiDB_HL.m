function table = igseq_runMultiDB_HL( directory )

list = dir([directory '/' '*.db']);
list.name
table = {};

for j = 1:length(list)
    db = [list(j).folder '/' list(j).name];
    table_tmp = {};
    j
    conn = igseq_db(db);
    
    %tmp = cellstr('Database');
    db_name = regexprep(db, '.*/', '');
    table_tmp = cellstr(db_name);
    
    q1 = 'SELECT chainType, count(*) from sequence group by chainType';
    
    t = get(fetch(exec(conn, q1)), 'data');
    tp = t';
    table_tmp = [table_tmp tp(2,:)];
    
    q2 = 'select pairedChain, chainType, count(*) from sequence group by pairedChain, chainType';
    
    t2 = get(fetch(exec(conn, q2)), 'data');
    tab2 = tabulate(t2(:,1));
    paired = tab2(find(cell2mat(tab2(:,2)) == 2));
    %tmp = cellstr('Paired');
    table_tmp = [table_tmp num2cell(length(paired))];
    
    q3 = 'select pairedChain, chainType, Vgene, Dgene, Jgene, count(*) from sequence where pairedChain != "null" and chainType == "H" group by pairedChain, chainType, Dgene';
    q4 = 'select pairedChain, chainType, Vgene, Dgene, Jgene, count(*) from sequence where pairedChain != "null" and chainType == "H" group by pairedChain, chainType, Jgene';
    q5 = 'select pairedChain, chainType, Vgene, Dgene, Jgene, count(*) from sequence where pairedChain != "null" and chainType == "H" group by pairedChain, chainType, Vgene';
    
    t3 = get(fetch(exec(conn, q3)), 'data');
    tab3 = tabulate(t3(:,1));
    multiH_D = tab3(find(cell2mat(tab3(:,2)) > 1));
    %tmp = cellstr('# Heavy Chain Dgene Mappings');
    table_tmp = [table_tmp num2cell(length(multiH_D))];
    
    
    t4 = get(fetch(exec(conn, q4)), 'data');
    tab4 = tabulate(t4(:,1));
    multiH_J = tab4(find(cell2mat(tab4(:,2)) > 1));
    %tmp = cellstr('# Heavy Chain Jgene Mappings');
    table_tmp = [table_tmp num2cell(length(multiH_J))];
    
    t5 = get(fetch(exec(conn, q5)), 'data');
    tab5 = tabulate(t5(:,1));
    multiH_V = tab5(find(cell2mat(tab5(:,2)) > 1));
    %tmp = cellstr('# Heavy Chain Vgene Mappings');
    table_tmp = [table_tmp num2cell(length(multiH_V))];
    
    q6 = 'select pairedChain, chainType, Vgene, Dgene, Jgene, count(*) from sequence where pairedChain != "null" and chainType == "L" group by pairedChain, chainType, Jgene';
    q7 = 'select pairedChain, chainType, Vgene, Dgene, Jgene, count(*) from sequence where pairedChain != "null" and chainType == "L" group by pairedChain, chainType, Vgene';
    
    t6 = get(fetch(exec(conn, q6)), 'data');
    tab6 = tabulate(t6(:,1));
    multiL_J = tab6(find(cell2mat(tab6(:,2)) > 1));
    %tmp = cellstr('# Light Chain Jgene Mappings');
    table_tmp = [table_tmp num2cell(length(multiL_J))];
    
    t7 = get(fetch(exec(conn, q7)), 'data');
    tab7 = tabulate(t7(:,1));
    multiL_V = tab7(find(cell2mat(tab7(:,2)) > 1));
    %tmp = cellstr('# Light Chain Vgene Mappings');
    table_tmp = [table_tmp num2cell(length(multiL_V))];
    
    table = vertcat(table, table_tmp);
end
    header = {'DB' 'H' 'L' 'Null' 'Paired' '# Heavy Chain Dgene Multi-Mappings' '# Heavy Chain Jgene Multi-Mappings' '# Heavy Chain Vgene Multi-Mappings' '# Light Chain Jgene Multi-Mappings' '# Light Chain Vgene Multi-Mappings'};
    table = vertcat(header, table);
    cell2csv('HL_db_stats.csv', table, ',', '2016', '.');
end
    
    