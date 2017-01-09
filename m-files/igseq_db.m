function dbh = igseq_db(db)

% IGSEQ_DB create database object from an igseq sqlite database.
%
%    conn = IGSEQ_DB(dbname) returns a database object given a sqlite
%    database name (including the full path) as a string input argument.
%
%    Requires Matlab Database Toolbox

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    url = strcat('jdbc:sqlite://',db);
    dbh  = database(db, '', '', 'org.sqlite.JDBC', url);
    
return