function a = igseq_nUniqueSeq(db)

% IGSEQ_NUNIQUESEQ number of unique sequences in each population of an analysis unit.
%
%    n = IGSEQ_NUNIQUESEQ(dbname) returns two-column cell array.  The first
%    column contains the population names and the second column contains
%    the number of unique sequences in the corresponding population.  The input
%    argument is a sqlite database name (including the full path) as a
%    string.
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_NSEQ, IGSEQ_NLIN, IGSEQ_ISOTYPES

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % open database and run query
    conn = igseq_db(db);
    a = get(fetch(exec(conn, 'select population, count(distinct NT_VDJ) from sequence group by population')), 'data');
    close(conn);
    a = sortrows(a, 1);

return