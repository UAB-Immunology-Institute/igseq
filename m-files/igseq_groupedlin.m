function tf = igseq_groupedlin(db)

% IGSEQ_GROUPEDLIN create grouped lineage table for an analysis unit
%
%    t = IGSEQ_GROUPEDLIN(dbname) returns a table as a cell array where
%    each row is a lineage.  The first argument is a sqlite database name
%    (including the full path) as a string, corresponding to the analysis
%    unit.  The columns are: lineage ID, v-gene, j-gene and then the number
%    of sequences within each population for each lineage in the database.
%    There are n+3 columns, where n is the number of populations.
%
%    Requires Matlab Database Toolbox

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % open database and get populations, lineages IDs and isotypes
    conn = igseq_db(db);
    q = 'SELECT population, lineageID, count(sequenceID), vgene, jgene FROM sequence GROUP BY lineageID, population ORDER BY lineageID';
    a = table2cell(get(fetch(exec(conn, q)), 'data'));
    close(conn);

    % get columns
    pop = a(:, 1);
    lin = cell2mat(a(:, 2));
    sz  = cell2mat(a(:, 3));
    v   = a(:, 4);
    j   = a(:, 5);
    
    % unique lineages and populations
    upop = unique(pop);
    [ulin, iu] = unique(lin);
    uv = v(iu);
    uj = j(iu);
    nlin = length(ulin);
    npop = length(upop);
    
    % convert query results to table
    m = zeros(nlin, npop);
    for i = 1:npop
        k = find(strcmp(upop{i}, pop));
        linpop = lin(k);
        [~, ~, ib] = intersect(linpop, ulin);
        m(ib, i) = sz(k);
    end
    
    % make output table
    t = [cellstr(num2str(ulin)), uv, uj, num2cell(m)];
    [~, isort] = sort(sum(m, 2), 'descend');
    ts = t(isort, :);
    th = [{'LineageID'}, {'V-gene'}, {'J-gene'}, upop'];
    tf = [th; ts];
    
    return
    
    
        
    
    
    
    
    
    
