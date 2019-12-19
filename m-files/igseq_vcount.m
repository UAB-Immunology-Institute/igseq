function X = igseq_vcount(db, type, varargin)

% IGSEQ_VJCOUNT tally v- usage by sequence or by lineage.
%
%    X = IGSEQ_VCOUNT(dbname, type) returns a data structure containing
%    counts for dfferent v-gene usage.  The first
%    argument is a sqlite database name (including the full path) as a
%    string.  The second input argument is either 'sequence' or 'lineage',
%    depending on whether you want to count the sequences or lineages with
%    a particular v- usage.
%
%    The return data structure contains five fields:
%        method:     the method for counting v- usage
%        count:      the population by v-gene matrix
%        population: corresponding population names for 1st dimension
%        vgene:      corresponding v-gene names for 2nd dimension
%
%
%    Requires Matlab Database Toolbox

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % querys to get v-gene anf j-genes from sqlite database for analysis unit
    q1 = 'SELECT population || "|" || Vgene FROM Sequence';
    q2 = 'SELECT population || "|" || Vgene FROM Sequence GROUP BY lineageID';

    % query to use based on user input
    if strcmp(type, 'sequence')
        q = q1;
    elseif strcmp(type, 'lineage')
        q = q2;
    else
        error('Second argument must be either ''sequence'' or ''lineage''');
    end

    % open a database connection, run query
    conn = igseq_db(db);
    a = get(fetch(exec(conn, q)), 'data');
    close(conn);

    % tabulate concatenated v-gene, population output
    t = tabulate(a);

    % split concatenated output into fields
    b = regexp(t(:, 1), '\|', 'split');
    b1 = cell(length(b), 2);
    for i = 1:length(b), b1(i, :) = b{i}; end

    % population, v-gene for usage table
    p  = b1(:, 1);
    v  = b1(:, 2);

    % number unique v-gene population combinations
    n = size(t, 1);

    % annotations for count matrix dimensions
    up = unique(p);
    uv = unique(v);

    % create count matrix
    c = zeros(length(up), length(uv));
    for i = 1:n
        c(strcmp(p{i}, up), strcmp(v{i}, uv)) = t{i, 2};
    end

    % output as either cell array or structure
    if nargin > 2 && varargin{1} == 1
        X = cell(2 + length(uv), 1 + length(up));
        X(3:end, 1) = uv;
        %X(2, 2:end) = repmat(uj', 1, length(up));
        X(1, 2:end) = reshape(repmat(up, 1, length(uj))', length(up) * length(uj), 1)';
        for i = 1:length(up)
            k1 = ((i - 1) * length(uj)) + 1;
            k2 = k1 + length(uj) - 1;
            X(3:end, 1 + (k1:k2)) = num2cell(squeeze(c(i, :, :)));
        end
    else
        X.method = type;
        X.count = c;
        X.population = up;
        X.vgene = uv;
        %X.jgene = uj;
    end

    return
