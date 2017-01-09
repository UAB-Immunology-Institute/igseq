function d = igseq_hillplot(db, q, ev, poplist, isolist, colors)

% IGSEQ_HILLPLOT Hill diversity profile plots over a range of orders.
%
%    IGSEQ_HILLPLOT(dbname, orders, evenflag, pop, iso, colors) creates
%    a diversity profile plot over the specified range of orders.  The
%    first argument is a sqlite database name (including the full path) as
%    a string.  The second argument is a vector of n orders across which to
%    compute the diversity index.  The third argument is a flag: a value of
%    1 plots evenness (diversity divided by zero-order-diversity); any
%    other value plots diversity.  The next three arguments are all
%    m-element cell vectors of strings corresponding to m desired plots.
%    The fourth argument is a list of population names, the fifth is a list
%    of isotypes (use '' if to include all), and the sixth is a list of
%    line colors.
%
%    h = IGSEQ_HILLPLOT(dbname, orders, evenflag, pop, iso, colors) creates
%    the figure and returns an n-order by m-population matrix of data that
%    corresponds to the plots.
%
%    Example:
%        >> db = 'Users/test/Desktop/myfile.igseqdb';
%        >> pop = {'42_N', '42_N', '42_SM', '42_SM'};
%        >> iso = {'', 'M', '', 'G'};
%        >> col = {'k', [.5 .5 .5], 'r', [1 .5 .5]};
%        >> igseq_hillplot(db, 0:.1:10, 0, pop, iso, col);
%
%    See:
%        - Hill (1973) Diversity and evenness: a unifying notation and its
%          consequences. Ecology 54:427-432.
%        - Stern et al (2014) B Cells populating the multiple sclerosis
%          brain mature in the draining cervical lymph nodes. Sci Trans Med
%          6(248):248ra107, DOI: 10.1126/scitranslmed.3008879.
%        - Grieff et al (2015) A bioinformatic framework for immune 
%          repertoire diversity profiling enables detection of 
%          immunological status. Genome Med 7:49,
%          DOI: 10.1186/s13073-015-0169-8.
%        
%    Requires Matlab Database Toolbox

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % error checking
    if length(poplist) ~= length(isolist), error('Population and isotype lists must be same length'); end
    if length(poplist) ~= length(colors),  error('Population and colors lists must be same length');  end
    if length(isolist) ~= length(colors),  error('Isaotype and colors lists must be same length');    end
    if ~iscellstr(poplist), error('Population list must be a cell array of strings'); end
    if ~iscellstr(isolist), error('Isotype list must be a cell array of strings');    end
    if ~iscell(colors),     error('Colors must be a cell');                           end

    % open database and get populations, lineages IDs and isotypes
    conn = igseq_db(db);
    a = get(fetch(exec(conn, 'select population, lineageID, isotype from sequence')), 'data');
    close(conn);
    
    % get columns from query result
    pop = a(:, 1);
    lin = cell2mat(a(:, 2));
    iso = a(:, 3);
    
    % number of curves to plot
    n = length(colors);        
    
    % initialize diversity data
    d = nan(length(q), n);
    
    % compute diversity profiles for each population
    for i = 1:n
        if ~isempty(isolist{i})
            ikeep = intersect(find(strcmp(poplist{i}, pop)), find(strcmp(isolist{i}, iso)));
        else
            ikeep = find(strcmp(poplist{i}, pop));
        end
        c = lin(ikeep);
        t = tabulate(c);
        t = sort(t(t(:, 2) > 0, 2), 'descend');                
        d(:, i) = hilldiv(t, q);
    end
    
    % evenness or diversity
    ylab = 'Diversity';
    if ev == 1
        d = d ./ repmat(d(1, :), length(q), 1);
        ylab = 'Evenness';
    end
    
    % figure constants and set up figure
    left = 70;
    right = 20;
    bottom = 50;
    top = 10;
    pw = 600;
    ph = 400;
    fw = left + pw + right;
    fh = bottom + ph + top;
    figure('position', [30 30 fw fh]);
    subplot('position', [left / fw, bottom / fh, pw / fw, ph / fh], 'nextplot', 'add');

    % add plots
    h = nan(n, 1);
    for i = 1:n
        h(i) = plot(q, d(:, i), '-', 'color', colors{i}, 'linewidth', 1.5);
    end
    set(gca,...
        'ticklength', [0 0],...
        'fontname', 'arial',...
        'fontsize', 12,...
        'yticklabel', num2str(get(gca, 'ytick')'),...
        'box', 'off');
    xlabel('Order', 'fontname', 'arial', 'fontsize', 16);
    ylabel(ylab, 'fontname', 'arial', 'fontsize', 16);
    
    % make legend
    leg = cell(n, 1);
    for i = 1:n
        if isempty(isolist{i})
            leg(i) = cellstr(['  ' poplist{i}]);
        else
            leg(i) = cellstr(['  ' poplist{i} ': Ig' isolist{i}]);
        end
    end
    hleg = legend(h, leg, 'location', 'northeast');
    set(hleg,...
        'fontname', 'arial',...
        'fontsize', 14,...
        'interpreter', 'none',...
        'box', 'off');
    
    return