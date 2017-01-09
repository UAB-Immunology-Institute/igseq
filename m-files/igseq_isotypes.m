function t = igseq_isotypes(db)

% IGSEQ_ISOTYPES isotype distributions for populations in an analysis unit.
%
%    t = IGSEQ_ISOTYPES(dbname) returns a table in the form of a cell array
%    where the rows are populations and the columns are isotypes.  The sum
%    of each row should add up to the total number of sequences for the
%    associated population.  The input argument is a sqlite database name
%    (including the full path) as a string.
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_NSEQ, IGSEQ_NLIN

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % open database and run query
    conn = igseq_db(db);
    a = get(fetch(exec(conn, 'select population, isotype, count(*) from sequence group by population, isotype')), 'data');
    close(conn);
    
    % unique popuations and isotypes
    upop = unique(a(:, 1));
    uiso = unique(a(:, 2));
    
    % initialize return table
    t = cell(length(upop) + 1, length(uiso) + 1);
    t(1, 1) = {''};
    t(2:end, 1) = upop;
    t(1, 2:end) = uiso;
    
    % populate table
    for i = 1:size(a, 1);
        r = 1 + find(strcmp(a(i, 1), upop));
        c = 1 + find(strcmp(a(i, 2), uiso));
        t(r, c) = a(i, 3);
    end
    
    % put in zeros where there are empty cells
    i = cellfun('isempty', t);
    i(1, 1) = 0;
    t(i) = num2cell(0);
    
return