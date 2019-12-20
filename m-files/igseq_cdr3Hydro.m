function d = igseq_cdr3Hydro(db)


    conn = igseq_db(db);
    a = get(fetch(exec(conn, 'select AA_CDR3, population from sequence')), 'data');
    close(conn);
    
    af = regexprep( a(:,2), "S\d_", '');
    
    pops = unique(a(:,2));
    for i = 1:length(pops)
        figure;
        popAA = a(find(strcmp(a(:,2),pops(i))));
        for j = 1:length(popAA)
            proteinpropplot(popAA{j,1});
            hold on;
        end
        p = regexprep(pops(i), '_', '-');
        title(p);
    end

return