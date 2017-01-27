function igseq_isoclonalityplots(db, iso, poplist, xmin, satcol)

% IGSEQ_ISOCLONALITYPLOTS show many clonality plots highlighting isotypes.
%
%    IGSEQ_ISOCLONALITYPLOTS(dbname, isotype, pops, xmin, color)
%    generates mulitple clonality plots for populations in the specified
%    database and colorizes lineages by percent composition of the
%    specified isotype.  The first argument is a sqlite database name
%    (including the full path) as a string.  The second argument is the
%    isotype to highlight (leaving out the "Ig" - should be "A", "E", "G",
%    "M", etc.  The third argument is a list of population names to plot -
%    they are plotted in the order that they appear in this list, from top
%    to bottom.  The fourth argument is the minimum value (as a percent) of
%    the x-axis (e.g., specifying "20" means that you are showing the
%    clonality plots from 20% to 100%).  The fifth argument is the color to
%    use for 100% composition of the specified isotype - it is used to
%    generate a gradient of the specified color (100%) to white (0%).
%
%    Example:
%        >> d = 'Users/test/Desktop/myfile.igseqdb';
%        >> pop = {'42_N', '42_SM', '42_PB'};
%        >> igseq_isoclonalityplots(d, 'A', pop, 20, [0 .6 .7]);
%
%    If you want to show all the populations in the database, you can use:
%        >> igseq_isoclonalityplots(d, 'A', igseq_upop(d), 20, [0 .6 .7];
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_UPOP, IGSEQ_CLONALITYPLOT, IGSEQ_ISOTYPES

% Chris Fucile
% Alex Rosenberg
% 1 May 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % make sure color specified ok
    if length(satcol) ~= 3, error('Color must be three-element vector'); end
    if max(satcol) > 1 || min(satcol) < 0, error('Elements of color vector must be between 0 and 1'); end
    
    % flip poplist so plots start from top, not bottom
    poplist = flipud(poplist(:))';
    
    % abort if populations not in db
    if ~isempty(setdiff(poplist, igseq_upop(db)))
        error('One or more specified population not in database');
    end
        
    % number of populations
    n = length(poplist);

    % colormap - 0 to 100
    
    if satcol(1) < 1, r = (1:-(1 - satcol(1))/100:satcol(1))'; else r = ones(101, 1); end
    if satcol(2) < 1, g = (1:-(1 - satcol(2))/100:satcol(2))'; else g = ones(101, 1); end
    if satcol(3) < 1, b = (1:-(1 - satcol(3))/100:satcol(3))'; else b = ones(101, 1); end
    cmap = [r g b];
        
    % populate histograms
    H = struct;
    conn = igseq_db(db);
    for i = 1:n
        disp(poplist{i});
        q = ['SELECT count(lineageID) AS C, sum(CASE WHEN isotype = "' iso '" THEN 1 ELSE 0 END) FROM Sequence WHERE population = "' poplist{i} '" GROUP BY lineageID ORDER BY C'];        
        a = cell2mat(get(fetch(exec(conn, q)), 'data'));
        H(i).lins  = a(:, 1);
        H(i).linsi = a(:, 2);
        H(i).color = cmap(1 + round((a(:, 2) ./ a(:, 1)) * 100), :);
    end
    close(conn);
    
    % threshold for drawing isotype colors (norm lin size - percent)
    thr = .05;
    
    % figure constants
    left = 12;
    pw = 1200;
    ph = 40;
    bot = 100;
    sp = 15; 
    top = 2;
    right = 140;
    cbarw = 250;
    cbarh = 30;
    coff = 50;
    cbot = 5;
    fw = left + pw + right;
    fh = bot + (n * ph) + ((n - 1) * sp) + top;
    
    % subplot preferences
    pref = struct(...
        'ticklength', [0 0],...
        'yticklabel', {{}},...
        'xlim', [xmin 100],...
        'ylim', [0 1],...
        'box', 'on',...
        'linewidth', 1.5,...
        'fontname', 'arial',...
        'fontsize', 12,...
        'layer', 'top');
    preft = struct(...
        'ticklength', [0 0],...
        'xticklabel', {{}},...
        'yticklabel', {{}},...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'box', 'off',...
        'xlim', [0 1],...
        'ylim', [0 1]);
    prefc = struct(...
        'ticklength', [0 0],...
        'xticklabel', {{}},...
        'yticklabel', {{}},...
        'box', 'on',...
        'linewidth', .5,...
        'xlim', [.5 101.5],...
        'ylim', [0 1],...
        'layer', 'top');
    
    % figure window
    figure('position', [20 20 fw fh]);
    
    % make each clonality bar
    for i = 1:n
        yoff = bot + ((i - 1) * (ph + sp));
        subplot('position', [(left + pw) / fw, yoff / fh, right / fw, ph / fh]);
            text(.05, .5, poplist{i}, 'fontname', 'arial', 'fontsize', 16, 'interpreter', 'none');
            set(gca, preft);
        subplot('position', [left / fw, yoff / fh, pw / fw, ph / fh], 'nextplot', 'add');
            normlin = 100 * H(i).lins ./ sum(H(i).lins);
            c = cumsum(normlin);
            ibig = find(normlin > thr);
            for q = 1:length(ibig)
                lup = c(ibig(q));
                llo = c(ibig(q) - 1);
                fill([llo llo lup lup llo], [0 1 1 0 0], H(i).color(ibig(q), :), 'edgecolor', 'none');
            end                
            plot([c'; c'], [0 1], '-', 'color', [1 1 1] * .2, 'linewidth', .25);
            set(gca, pref);
            if i > 1, set(gca, 'xticklabel', {}); end
            if i == 1, xlabel({' '; 'Cumulative Percentage of Sequences'}, 'fontname', 'arial', 'fontsize', 16); end
    end
    
    % make colorbar
    subplot('position', [(fw - (cbarw + coff)) / fw, cbot / fh, cbarw / fw, cbarh / fh], 'nextplot', 'add');
    for i = 1:101
        x = [-.5 .5 .5 -.5 -.5] + i;
        y = [0 0 1 1 0];
        fill(x, y, cmap(i, :), 'edgecolor', cmap(i, :));
    end
    set(gca, prefc);
    title(['Percent Ig' iso], 'fontname', 'arial', 'fontsize', 16);
    
return    

