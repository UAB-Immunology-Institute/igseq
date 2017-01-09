function igseq_muthistplot(varargin)

% IGSEQ_MUTHISTPLOT mutation frequency histograms for an analysis unit.
%
%    IGSEQ_MUTHISTPLOT(dbname, type, binsize, maxfrq, pops, isos, cols)
%    plots the distributions of V-gene mutation frequencies for the
%    specified population, isotype combinations in the analysis unit at the
%    specified colors.  The first argument is a sqlite database name
%    (including the full path) as a string.  The second argument is one of:
%    "total", "silent" or "nonsilent" to specify the type of mutations to
%    consider.  The third argument is the histogram bin size (in
%    percentages) - "1" is a typical value to use.  The fourth argument is
%    the maximum value for the histogram (e.g. 40).  The return value is a
%    cell array.  The fifth, sixth and seventh arguments are lists of
%    populations, isotypes and line colors, respectively - all cell arrays
%    of strings.
%
%    Example:
%        >> db = 'Users/test/Desktop/myfile.igseqdb';
%        >> pop = {'42_N', '42_N', '42_SM', '42_SM'};
%        >> iso = {'', 'M', '', 'G'};
%        >> col = {'k', [.5 .5 .5], 'r', [1 .5 .5]};
%        >> igseq_muthistplot(db, 'silent', 1, 40, pop, iso, col);
%
%    IGSEQ_MUTHISTPLOT(dbname, pops, isos, cols) creates the same plot as
%    above, except that the first argument is the output of igseq_muthist.
%
%    Example:
%        >> h = igseq_muthist(db, 'silent', 1, 40);
%        >> igseq_muthistplot(h, pop, iso, col);
%
%    Note that the second usage will generate the plots much faster because
%    there is no overhead for database calls.
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_MUTHIST

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    if nargin == 7
        db = varargin{1};
        t = varargin{2};
        binsize = varargin{3};
        maxfrq = varargin{4};
        pop = varargin{5};
        iso = varargin{6};
        col = varargin{7};
    elseif nargin == 4
        a = varargin{1};
        pop = varargin{2};
        iso = varargin{3};
        col = varargin{4};
    else
        error('Requires either 4 or 7 input arguments');
    end

    % error checking
    if length(pop) ~= length(iso), error('Population and isotype lists must be same length'); end
    if length(pop) ~= length(col), error('Population and colors lists must be same length');  end
    if length(iso) ~= length(col), error('Isaotype and colors lists must be same length');    end
    if ~iscellstr(pop), error('Population list must be a cell array of strings'); end
    if ~iscellstr(iso), error('Isotype list must be a cell array of strings');    end
    if ~iscell(col),    error('Colors must be a cell');                           end
    n = length(pop);

    % get histograms for all populations and isotypes if 7 args
    if nargin == 7
        a = igseq_muthist(db, t, binsize, maxfrq);
    end
    t = a{1, 1};
    
    % available pops and isotypes
    allpop = a(1, 4:end);
    alliso = a(2, 4:end);
    
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
    
    % convert into column indices to find which histograms to plot
    x = cell2mat(a(3:end, 3));
    h = nan(n, 1);
    for i = 1:n
        idx = intersect(find(strcmp(pop{i}, allpop)), find(strcmp(iso{i}, alliso)));
        y = cell2mat(a(3:end, idx + 3));
        h(i) = plot(x, y, '-', 'color', col{i}, 'linewidth', 1);      
    end
    set(gca,...
        'ticklength', [0 0],...
        'fontname', 'arial',...
        'fontsize', 12,...
        'yticklabel', num2str(get(gca, 'ytick')'),...
        'box', 'off');
    xlabel(['% V Mutation Frequency (' t ')'], 'fontname', 'arial', 'fontsize', 16);
    ylabel('Number of Sequences', 'fontname', 'arial', 'fontsize', 16);
    
    % make legend
    leg = cell(n, 1);
    for i = 1:n
        if isempty(iso{i})
            leg(i) = cellstr(['  ' pop{i}]);
        else
            leg(i) = cellstr(['  ' pop{i} ': Ig' iso{i}]);
        end
    end
    hleg = legend(h, leg, 'location', 'northeast');
    set(hleg,...
        'fontname', 'arial',...
        'fontsize', 14,...
        'interpreter', 'none',...
        'box', 'off');

    return
