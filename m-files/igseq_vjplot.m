function maxp = igseq_vjplot(db, type, varargin)

% IGSEQ_VJPLOT make VJ-usage heatmaps for each population in analysis unit.
%
%    c = IGSEQ_VJPLOT(dbname, type) creates a VJ-usage heat map for each
%    population in the specified database.  Rows are v-genes and columns
%    are j-genes. The first argument is a sqlite database name (including
%    the full path) as a string.  The second input argument is either
%    'sequence' or 'lineage', depending on whether you want to count the
%    sequences or lineages with a particular v- and j-gene usage.  The
%    frequency of a particular v-gene and j-gene combination is represented
%    by color using the parula color scale (black = 0, blue = low
%    frequency, yellow = high frequency).  The maximum of the color scale
%    is set to the 99.5th percentile of all frequencies in the data set.
%    The return argument is the maximum of the color scale.
%
%    c = IGSEQ_VJPLOT(dbname, type, cmax) creates plots as above, but
%    allows the user to specify the maximum frequency (0-100) for the color
%    scale for all plots.
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_VJCOUNT

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % get v- and j-gene usage counts
    a = igseq_vjcount(db, type);
    
    % matrix dimensions
    nv = length(a.vgene);
    nj = length(a.jgene);
    np = length(a.population);

    % compute frequencies
    for i = 1:np
        dat = squeeze(a.count(i, :, :));
        a.frq(i, :, :) = 100 * dat / sum(dat(:));        
    end
    
    % determine max color scale - either entered, or 99.5th percentile
    nargs = nargin;
    if nargs > 2
        maxp = varargin{1};
    else
        maxp = prctile(a.frq(:), 99.5);
    end
        
    % figure constants
    pix = 12;
    pw = nj * pix;
    ph = nv * pix;
    bottom = 60;
    left = 2;
    right = 90;
    top = 60;
    fw = left + pw + right;
    fh = bottom + ph + top;
    
    % colormap (set lowest color to black)
    cmap = parula(1024);
    cmap(1, :, :) = [0 0 0];
    
    % fontprefs
    fpref = struct('fontname', 'arial', 'interpreter', 'none');
    
    % create one figure per population
    for f = 1:np
        
        % make figure
        figure('position', [60 + ((f - 1) * fw), 60, fw, fh], 'color', 'w');
        
        % figure label
        subplot('position', [0, (bottom + ph) / fh, 1, top / fh]);
        text(.5, .8, a.population{f}, fpref, 'fontsize', 16, 'horizontalalignment', 'center');
        text(.5, .3, ['V-J usage by ' a.method], fpref, 'fontsize', 12, 'horizontalalignment', 'center');
        set(gca, 'visible', 'off', 'xlim', [0 1], 'ylim', [0 1]);
        
        % v-gene labels
        subplot('position', [(left + pw) / fw, bottom / fh, right / fw, ph / fh]);
        for i = 1:nv, text(.05, i, a.vgene{i}, fpref, 'fontsize', 10); end
        set(gca, 'visible', 'off', 'xlim', [0 1], 'ylim', [.5 nv + .5], 'ydir', 'reverse');
        
        % j-gene labels
        subplot('position', [left / fw, 0, pw / fw, bottom / fh]);
        for i = 1:nj, text(i, .95, a.jgene{i},fpref, 'fontsize', 10, 'rotation', 270); end
        set(gca, 'visible', 'off', 'xlim', [.5 nj + .5], 'ylim', [0 1]);
        
        % heat map
        subplot('position', [left / fw, bottom / fh, pw / fw, ph / fh]);
        im = imresize(squeeze(a.frq(f, :, :)), pix, 'nearest');
        imagesc(im, [0 maxp]);
        colormap(cmap);
        set(gca, 'visible', 'off', 'ydir', 'reverse');
        axis off;
        
    end
        
        

