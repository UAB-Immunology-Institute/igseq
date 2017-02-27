function t = igseq_groupedliniso(db, varargin)

% IGSEQ_GROUPEDLINISO create grouped lineage table for an analysis unit
%
%    t = IGSEQ_GROUPEDLINISO(dbname) returns a table as a cell array where
%    each row is a lineage.  The first argument is a sqlite database name
%    (including the full path) as a string, corresponding to the analysis
%    unit.  The columns are: lineage ID and then the number of sequences
%    within each population/isotype combination for each lineage in the
%    database.  There are n+1 columns, where n is the number of population-
%    isotype combinations.  The contents of this array can be cut and
%    pasted into Excel.
%
%    NOTE: singleton lineages (lineages with only one sequence in only one
%    population) are not included in the tabular output.
%
%    t = IGSEQ_GROUPEDLINISO(dbname, pops, iso) returns the same output as
%    described above, except it only returns rows where there are more than
%    zero sequences in the specified populations AND for the specified
%    isotype.  The list of populations (second argument) must be a cell
%    array of strings.  The isotype (third argument) must be a string.
%    These two arguments can be empty.  If both are empty, then the
%    function behaves as if only the database name was provided.
%
%    to get a table of lineages broken down broken down by population and
%    isootype with no filtering:
%
%       >> t = igseq_groupedliniso(db);
%
%    to get this table filtered to only show lineages where there are
%    sequences in the specified populations for IgE:
%
%       >> t = igseq_groupedliniso(db, {'POP2', 'POP3', 'POP4'}, 'E');
%
%    to get this table filtered to only show lineages where there are
%    sequences in the specified populations for any isotype:
%
%       >> t = igseq_groupedliniso(db, {'POP2', 'POP3', 'POP4'}, '');
%
%    to get this table filtered to only show lineages where there are
%    sequences in any population but for IgE:
%
%       >> t = igseq_groupedliniso(db, {}, 'E');
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_GROUPEDLIN, IGSEQ_ISOTYPES

% Chris Fucile
% Alex Rosenberg
% 15 Aug 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % error checking
    if nargin ~= 1
        if nargin == 3
            popfilter = varargin{1};
            if ~iscellstr(popfilter), error('pop filter must be cell array of strings'); end
            isofilter = varargin{2};
            if ~ischar(isofilter), error('isotype must be a string'); end
        else
            error('must supply either one or three arguments');
        end
    else
        popfilter = {};
        isofilter = '';
    end    

    % open database, tabulate lineages by  population and isotype
    conn = igseq_db(db);
    q = 'SELECT lineageID, population, isotype, COUNT(*) FROM Sequence GROUP BY lineageID, population, isotype';    
    a = get(fetch(exec(conn, q)), 'data');
    qi = 'SELECT DISTINCT(isotype) FROM SEQUENCE';
    iso = get(fetch(exec(conn, qi)), 'data');
    close(conn);
    
    % put query results in a table
    L0.lin = cell2mat(a(:, 1));
    L0.pop = a(:, 2);
    L0.iso = a(:, 3);
    L0.seq = cell2mat(a(:, 4));

    % get rid of singleton lineages (ending in ".1")
    z = L0.lin ./ round(L0.lin);
    L.lin = L0.lin(z == 1);
    L.pop = L0.pop(z == 1);
    L.iso = L0.iso(z == 1);
    L.seq = L0.seq(z == 1);
    n = length(L.lin);

    % unique lineages    
    ulin = unique(L.lin);
    nlin = length(ulin);

    % unique populations
    pop = igseq_upop(db);
    
    % final error check
    if ~isempty(popfilter) && (length(unique(popfilter)) ~= length(intersect(popfilter, pop)))
        error('one or more populations not in database');
    end
    if ~isempty(isofilter) && isempty(find(strcmp(isofilter, iso), 1))
        error('isoype not in database');
    end
    
    % number of populations ans isotypes
    npop = length(pop(:, 1));
    niso = length(iso);

    % initialize tabulation matrix
    x = zeros(npop, nlin, niso);
    for i = 1:n
        kpop = strcmp(L.pop{i}, pop(:, 1));
        kiso = strcmp(L.iso{i}, iso);
        klin = L.lin(i) == ulin;
        x(kpop, klin, kiso) = L.seq(i);
    end
    
    % flatten to a lineage X (population/isotype) 2-D matrix    
    y = [];
    for i = 1:length(pop)
        y = [y squeeze(x(i, :, :))]; %#ok<AGROW>
    end

    % column headings for output table
    poplabel = reshape(repmat(pop', niso, 1), niso * npop, 1)';
    isolabel = repmat(iso', 1, npop);
    
    % combine labels
    clabel = cell(1, niso * npop);
    for i = 1:length(clabel)
        clabel(i) = cellstr([poplabel{i} '_Ig' isolabel{i}]);
    end
    
    % identify populations for filtering based on user input
    if isempty(popfilter) && ~isempty(isofilter)
        ifilt = find(ismember(isolabel, isofilter));
    elseif ~isempty(popfilter) && isempty(isofilter)
        ifilt = find(ismember(poplabel, popfilter));
    elseif ~isempty(popfilter) && ~isempty(isofilter)
        ifilt = intersect(...
            find(ismember(poplabel, popfilter)),...
            find(ismember(isolabel, isofilter)));
    else
        ifilt = 1:(niso * npop);
    end
    
    % identify lineages that pass filter
    ipass = find(sum(y(:, ifilt), 2) > 0);
    yfilt = y(ipass, :);
    ulinfilt = ulin(ipass);
    nlinfilt = length(ulinfilt);
               
    % generate ouptut table
    t = cell(1 + nlinfilt, 1 + (niso * npop));
    t(1, :) = [{'lineageID'} clabel];
    t(2:end, 1) = num2cell(ulinfilt);
    t(2:end, 2:end) = num2cell(yfilt);

    return