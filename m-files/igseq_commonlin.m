function r = igseq_commonlin(db, s)

% IGSEQ_COMMONLIN tabulate degree of shared lineages between populations
%
%    t = IGSEQ_COMMONLIN(dbname, type) returns a square matrix (a cell
%    array) that summarizes the number or percent of shared lineages
%    between populations.  The first argument is a sqlite database name
%    (including the full path) as a string, corresponding to the analysis
%    unit.  The second argument is either 'count' or 'percent', specifying
%    how the data should be summarized.  If 'count', then the cells of the
%    return matrix will be the number of lineages that are common between
%    the two populations.  The values along the diagonal are simply the
%    nubmer of lineages in the population.  Only the lower triangle of this
%    matrix is filled.  If 'percent', then the values in the cells
%    correspond to the percent of shared lineages between the two
%    populations relative to the population indicated by the row.  The
%    first row and column of the output matrix indicates the population
%    names.
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_MHOVERLAP, IGSEQ_GROUPEDLIN

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % error check
    if isempty(find(strcmp(s, {'count', 'percent'}), 1))
        error('Second argument must be either ''count'' or ''percent''');
    end
    
    % get grouped lineage table
    t = igseq_groupedlin(db);

    % extract population names
    pop = t(1, 4:end);

    % number of populations
    n = length(pop);

    % count number of lineages common to pairs of populations
    data = cell2mat(t(2:end, 4:end));
    c = zeros(n);
    for i = 1:n
        for j = 1:i
            c(i, j) = length(find(prod(data(:, [i j]), 2)));
        end
    end

    % compute percentages of lineages of source (row) populations
    c1 = tril(c, -1)' + c;
    cp = 100 * c1 ./ repmat(diag(c1), 1, n);

    % prepare output cell array
    r = cell(n + 1);
    r(1, 1) = {''};
    r(1, 2:end) = pop;
    r(2:end, 1) = pop';

    % return type of data requested
    if strcmp(s, 'count')
        r(2:end, 2:end) = num2cell(c);
    else
        r(2:end, 2:end) = num2cell(cp);
    end

    return



