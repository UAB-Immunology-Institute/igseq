function igseq_alluv2(dbname, p, varargin)

% IGSEQ_ALLUV2 Create a 2-population alluvial clonality plot
%
%    IGSEQ_ALLUV2(dbname, pop) creates two clonality plots where common
%    lineages in the two populations are connected by bezier ribbons.  The
%    first argument is a sqlite database name (including the full path)
%    as a string.  The second areument is a 2-element cell vector of
%    population names as appears in the specified database.  An optional
%    third argument (a 2-element cell vector) specifies alternate
%    population names to print on the plot.
%
%    Usage:
%
%        IGSEQ_ALLUV2('/path/to/my.db', {'54_7', '54_11'})
%
%        IGSEQ_ALLUV2('/path/to/my.db', {'54_7', '54_11'}, {'mem', 'PC'})
%
%    Requires Matlab Database Toolbox

% Chris Fucile
% Alex Rosenberg
% 26 Oct 2018
% University of Alabama at Birmingham, Department of Microbiology
% Copyright (c) 2018. All rights reserved.
% This software is offered with no guarantees of any kind.

    % if specified population names different from those in database
    if nargin > 2
        t = varargin{1};
    else
        t = p;
    end

    % establish connection to database and execute query
    q = [...
        'SELECT population, lineageID, count(sequenceID) FROM ',...
        'Sequence WHERE population IN ("' p{1} '", "' p{2} '") ',...
        'GROUP BY population, lineageID'];
    conn = igseq_db(dbname);
    a = get(fetch(exec(conn, q)), 'data');
    close(conn);
    pop = a(:, 1);
    lin = cell2mat(a(:, 2));
    sz = cell2mat(a(:, 3));

    % separate by pop
    i1 = strcmp(pop, p{1});
    i2 = strcmp(pop, p{2});
    lin1 = lin(i1);
    lin2 = lin(i2);
    sz1 = sz(i1);
    sz2 = sz(i2);

    % combine into one grouped lineage table
    T.lin = unique([lin1; lin2]);
    T.p1 = zeros(length(T.lin), 1);
    T.p2 = zeros(length(T.lin), 1);
    [~, j1, k1] = intersect(lin1, T.lin);
    [~, j2, k2] = intersect(lin2, T.lin);
    T.p1(k1) = sz1(j1);
    T.p2(k2) = sz2(j2);
    
    % express as percent total sequences
    T.p1 = 100 * T.p1 / sum(T.p1);
    T.p2 = 100 * T.p2 / sum(T.p2);

    % figure constants
    cw = 60;
    bw = 400;
    left = 100;
    right = 100;
    bot = 50;
    top = 50;
    ph = 880;
    fw = left + cw + bw + cw + right;
    fh = bot + ph + top;
    posname = [0, 0, 1, bot / fh];   
    pos1t = [left             / fw, (bot + ph) / fh, cw / fw, top / fh];
    pos2t = [(left + cw + bw) / fw, (bot + ph) / fh, cw / fw, top / fh];
    pos1  = [left             / fw, bot / fh, cw / fw, ph / fh];    
    pos2  = [(left + cw + bw) / fw, bot / fh, cw / fw, ph / fh];
    posc  = [(left + cw)      / fw, bot / fh, bw / fw, ph / fh];
    pref = struct(...
        'ticklength', [0 0],...
        'xticklabel', {{}},...
        'yticklabel', {{}},...
        'linewidth', 1,...
        'box', 'on',...
        'xlim', [0 1],...
        'ylim', [0 100]);
    textpref = struct(...
        'fontname', 'arial',...
        'fontsize', 20,...
        'fontweight', 'bold',...
        'horizontalalignment', 'center',...
        'interpreter', 'none');
    titpref = struct(...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'box', 'off');
    maxcol = 300;
    cmap = winter(maxcol);

    % make structure array
    L = struct;
    for i = 1:2
        k = find(T.(['p' num2str(i)]) > 0);
        L(i).lin = T.lin(k);
        L(i).sz = T.(['p' num2str(i)])(k);
    end

    % sort each of the clonality datasets
    S = struct;
    for i = 1:length(L)
        [~, isort] = sort(L(i).sz);
        S(i).lin = L(i).lin(isort);
        S(i).sz  = L(i).sz(isort);
        S(i).y   = [0; cumsum(S(i).sz)];
    end

    % add axis positions
    S(1).pos = pos1;
    S(2).pos = pos2;

    % make figure
    figure('position', [5 5 fw fh]);

    % database name
    subplot('position', posname);
    name2 = regexp(dbname, '/', 'split');
    text(.5, 50, name2{end},...
        'fontname', 'courier',...
        'fontsize', 24,...
        'horizontalalignment', 'center',...
        'interpreter', 'none');
    set(gca, pref, titpref);

    % titles
    subplot('position', pos1t);
    text(.5, 70, t{1}, textpref);
    set(gca, pref, titpref);
    subplot('position', pos2t);
    text(.5, 70, t{2}, textpref);
    set(gca, pref, titpref);

    % make arcs
    subplot('position', posc, 'nextplot', 'add');
    clin = intersect(S(1).lin, S(2).lin);
    cmap = cmap(randperm(size(cmap, 1)), :);
    for i = 1:length(clin)
        i1 = find(S(1).lin == clin(i));
        i2 = find(S(2).lin == clin(i));
        [x12, y12] = bezarc(...
            0, S(1).y(i1), S(1).y(i1 + 1),...
            1, S(2).y(i2), S(2).y(i2 + 1), .8);
        icol = mod(i, maxcol) + 1;
        hf = fill(x12, y12, cmap(icol, :),...
            'edgecolor', cmap(icol, :),...
            'linewidth', .25);
        alpha(hf, .5);
    end
    set(gca, pref, 'xcolor', 'w', 'ycolor', 'w', 'box', 'off');
    uistack(gca, 'bottom');    

    % make clonality plot
    yt = 0:5:100;
    ytl = cell(length(yt), 1);
    for i = 1:length(yt), ytl(i) = cellstr(num2str(yt(i))); end
    h = nan(length(S), 1);
    for i = 1:length(S)
        h(i) = subplot('position', S(i).pos, 'nextplot', 'add');
        plot([0 1], [S(i).y S(i).y], 'k-', 'linewidth', .25);
        set(gca, pref);
    end
    set(h(1), 'ytick', yt, 'yticklabel', ytl, 'fontname', 'arial', 'fontsize', 14);
    axes(h(1));
    hy = ylabel('Cumulative Percentage of Sequences', 'fontname', 'arial', 'fontsize', 18);
    set(hy, 'position', [-.8 50 1]);

    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');
    
return
