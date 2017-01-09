function a = igseq_commonlintop(db, idx)

% IGSEQ_COMMONLINTOP shared lineages between populations for top lineages
%
%    t = IGSEQ_COMMONLINTOP(dbname, idx) returns a cell array that
%    summarizes percentage of the top n lineages in a "source" population
%    that share sequences with lineages in a "destination" population.  The
%    first argument is a sqlite database name (including the full path) as
%    a string, corresponding to the analysis unit.  The second argument is
%    a list of numbers corresponding to the top numbers (n's) of lineages.
%    For example, if you wanted to characterize the percent shared lineages
%    corresponding to the top 10, 20, 30, 40, 50 lineages, you could enter
%    10:10:50 (start with 10, advance in steps of 10, and end with 50).
%    Alternatively, you could enter a vector of values like [10 50 100].
%    The return value is a cell array where the first column indicates the
%    number of top lineages.  The first row indicates the source
%    population.  The second row indicates the destination population.  For
%    example, if the ith row and first column was "30", the jth column,
%    first row was "MEMORY" and jth column second row was "ASC", then the
%    value for (i, j) would be "the percentage of the top 30 lineages
%    in the MEMORY population that share related sequences to lineages in
%    the ASC population."  You can plot any column agains the first column
%    to see whether the larger lineages tend to be more connected.  Be
%    aware that for polyclonal populations, the ranking of lineage sizes
%    can be more arbitrary compared to oligoclonal populations.
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_MHOVERLAP, IGSEQ_GROUPEDLIN, IGSEQ_COMMONLIN

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.


    % get grouped lineage table
    t = igseq_groupedlin(db);

    % extract population names
    pop = t(1, 4:end);

    % number of populations
    n = length(pop);

    % number of "n" in "top n lineages"
    ni = length(idx);
    
    % initialize output
    r = nan(ni, n * n);
    
    % initialize population headers
    psource = cell(1, n * n);
    pdest   = cell(1, n * n);
    
    % source and destination indices
    k = [reshape(repmat((1:n), n, 1), n^2, 1) repmat((1:n)', n, 1)];
    nk = size(k, 1);
    
    % step through pairs
    for i = 1:nk
        psource(i) = pop(k(i, 1));        
        pdest(i)   = pop(k(i, 2));
        j1 = 3 + find(strcmp(pop, psource(i)));
        j2 = 3 + find(strcmp(pop, pdest(i)));
        d = flipud(sortrows(cell2mat(t(2:end, [j1 j2])), 1));
        for q = 1:ni
            r(q, i) = 100 * length(find(prod(d(1:idx(q), :), 2) > 0)) / idx(q);
        end
    end
    
    % prepare output table
    a = cell(ni + 2, (n * n) + 1);
    a(1, 1) = {'source'};
    a(2, 1) = {'dest'};
    a(1, 2:end) = psource;
    a(2, 2:end) = pdest;
    a(3:end, 1) = num2cell(idx(:));
    a(3:end, 2:end) = num2cell(r);
    
   return