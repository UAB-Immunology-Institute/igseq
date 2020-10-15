function a = igseq_summary(db)

% IGSEQ_SUMMARY brief summary of database.
%
%    n = IGSEQ_SUMMARY(dbname) returns basic info about database.
%
%    Requires Matlab Database Toolbox

% Chris Fucile
% Alex Rosenberg
% 9 Oct 2018
% UAB Department of Microbiology and Informatics Institute
% Copyright (c) 2018. All rights reserved.
% This software is offered with no guarantees of any kind.

    % open database and run query
    conn = igseq_db(db);
    a = cell2mat(get(fetch(exec(conn, 'select count(*) from Sequence')), 'data'));
    close(conn);

return