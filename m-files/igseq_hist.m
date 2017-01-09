function s = igseq_hist(db, pop)

% IG_HIST returns clonality distribution histogram
%
% This function requires connection to the sequencing database.
%
% INPUT ARGUMENTS:
%
%    1: path to SQLite db
%    2: abbreviation of cell population (as appears in database)
%
% RETURN VALUE:
%
%    an nx2 matrix where column 1 are lineage sizes and column 2 are
%    numbers of lineages with the corresponding sizes
%
% USAGE EXAMPLE:
%
%    >> H = igseqhist('my analysis unit','PB');
%
% by Alex Rosenberg
% University of Rochester Medical Center
% Division of Allergy Immunology and Rheumatology
% 23-APR-2014

    % open database connection
    conn = igseq_db(db);

    % sub query
    subq = ['(SELECT count(sequenceID) AS sz FROM Sequence',...
        ' WHERE population = ' '''' pop '''',...
        ' GROUP BY lineageID)'];
    
    % main query
    q = ['SELECT X.sz, count(X.sz) FROM ' subq ' X GROUP BY X.sz ORDER BY X.sz'];

    % get result
    s = cell2mat(get(fetch(exec(conn, q)), 'data'));
        
    % close connection
    close(conn);

return
