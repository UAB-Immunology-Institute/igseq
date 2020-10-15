function u = igseq_upop(db)

% IGSEQ_UPOP list available populations in an analysis unit database.
%
%    u = IGSEQ_UPOP(dbname) returns a list of population names in the
%    specified analysis unit database.  The input argument is a sqlite
%    database name (including the full path) as a string.
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_NSEQ, IGSEQ_NLIN, IGSEQ_ISOTYPES, IGSEQ_SUMMARY

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % open database and run query
    conn = igseq_db(db);
    u = sort(get(fetch(exec(conn, 'select distinct(population) from sequence')), 'data'));
    close(conn);

return