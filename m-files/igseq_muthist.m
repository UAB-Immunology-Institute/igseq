function m = igseq_muthist(db, t, binsize, maxfrq)

% IGSEQ_MUTHIST mutation frequency histograms for an analysis unit.
%
%    h = IGSEQ_MUTHIST(dbname, type, binsize, maxfrq) returns a table of
%    binned mutation frequency histograms for all population/isotype
%    combinations within an analysis unit (database).  The first argument
%    is a sqlite database name (including the full path) as a string.  The
%    second argument is one of: "total", "silent" or "nonsilent" to specify
%    the type of mutations to consider.  The third argument is the
%    histogram bin size (in percentages) - "1" is a typical value to use.
%    The fourth argument is the maximum value for the histogram (e.g. 40).
%    The return value is a cell array.  The data columns are the histograms
%    for particular population/isotype combinations.  The cells are the
%    counts per bin.  The first two columns are the bin start and bin end,
%    respectively (e.g. a bin for the number of sequences that had mutation
%    frequencies between 3% and 4%).  The third column is the bin center
%    (3.5% in the previous example).  The first row indicates the
%    population.  The second row indicates the isotype.  No isotype label
%    indicates that all isotypes (including "U") were considered.
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_ISOTYPES, IGSEQ_MUTHISTPLOT

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    
    % type of mutation
    switch t
        case      'total',     f = 'V_Mutations';
        case      'silent',    f = 'V_SilentMutations';
        case      'nonsilent', f = 'V_NonsilentMutations';
        otherwise
            error('type must be one of: ''total'', ''silent'', ''nonsilent''');
    end
    
    % open database and run query to get all data
    q = ['select population, isotype, V_Nucleotides, ' f ' from sequence'];
    conn = igseq_db(db);
    a = get(fetch(exec(conn, q)), 'data');
    close(conn);

    % get unique populations and isotypes
    upop = unique(a(:, 1));
    uiso = unique(a(:, 2));
    npop = length(upop);
    niso = length(uiso);
    n = npop * niso;
    
    % compute frequencies
    frq  = 100 * cell2mat(a(:, 4)) ./ cell2mat(a(:, 3));
    
    % remove U from isotype if present
    uiso = [{''}; setdiff(uiso, 'U')];
    
    % population, isotype combinations
    poplabel = reshape(repmat(upop, 1, niso)', 1, n);
    isolabel = repmat(uiso, npop, 1)';
    
    % histogram bins
    edges = 0:binsize:maxfrq;
    binctr = edges(1:(end - 1)) + (binsize / 2);
    nbin = length(binctr);
    
    % initialize histogram data
    h = zeros(nbin, n);
    
    % set through population, isotype combinations and make histograms
    for i = 1:n
        if ~isempty(isolabel{i})
            k = intersect(find(strcmp(a(:, 1), poplabel{i})), find(strcmp(a(:, 2), isolabel{i})));
        else
            k = find(strcmp(a(:, 1), poplabel{i}));
        end
        if ~isempty(k)
            h0 = histc(frq(k), edges)';
            h(:, i) = h0(1:(end - 1));
        end
    end
    
    % prepare output
    m = cell(nbin + 2, n + 3);
    m(3:end, 1) = num2cell(edges(1:(end - 1)))';
    m(3:end, 2) = num2cell(edges(2:end))';
    m(3:end, 3) = num2cell(binctr)';
    m(1, 4:end) = poplabel;
    m(2, 4:end) = isolabel;
    m(1, 1:3) = {t, '', ''};
    m(2, 1:3) = {'bin_start', 'bin_end', 'bin_ctr'};
    m(3:end, 4:end) = num2cell(h);

    return
