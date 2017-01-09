function [nseq, nlin] = igseq_clonalityplot(db, pop, label, col, miny, maxx)

% IGSEQCLONALITYPLOT makes a clonality plot for a sample/population
%
% This function requires connection to the sequencing database.
%
% INPUT ARGUMENTS:
%
%    1: path to SQLite db
%    2: abbreviation of cell population (as appears in database)
%    3: title of plot
%    4: plot color
%    5: minimum y-axis percentage value (0 - 100)
%    6: maximum x value to use ([] for autoscale)
%
% RETURN VALUES:
%
%    1: number of sequences
%    2: number of lineages
% 
% USAGE EXAMPLE:
%
%    >> igseqclonalityplot('/test.db', 'PB', '185-PB', 'b', 30, []);
%
% by Alex Rosenberg
% University of Rochester Medical Center
% Division of Allergy Immunology and Rheumatology
% 24-FEB-2014
%
% see also IGSEQHIST, IGSEQDX

    % get histogram from database
    H = igseq_hist(db, pop);

    % figure constants (in pixels)
    sw = 180;
    sh = 550;
    spleft = 42;
    spright = 15;
    spbot = 32;
    sptop = 10;    
    titleh = 50;
    
    % font sizes
    fonttic = 24;
    fonttitle = 32;
 
    % calculate figure dimensions
    fw = spleft + sw + spright;
    fh = spbot + sh + sptop + titleh;
    
    % plot positions
    posMain  = [...
        (spleft / fw),...
        (spbot / fh),...
        (sw / fw),...
        (sh / fh)];
    posTitle = [...
        (spleft / fw),...
        ((spbot + sh + sptop) / fh),...
        (sw / fw),...
        (titleh / fh)];
        
    % get data
    rawlinsize = rldecode(H(:, 1)', H(:, 2)');
    cumseq = cumsum(rawlinsize);
    nseq = max(cumseq);
    nlin = sum(H(:, 2));
    cpercent = 100 * cumseq / nseq;
    normlinsize = 100 * rawlinsize / nseq;
        
    % if necessary, autoscale x-axis
    if isempty(maxx), maxx = max(normlinsize) * 1.002; end
       
    % create figure
    figure('position', [60 60 fw fh]);

    % plot to hold figure title
    subplot('position', posTitle);    
    text(1, 1.4, label,...
        'fontname', 'arial',...
        'fontsize', fonttitle,...
        'fontweight', 'bold',...
        'horizontalalignment', 'center',...
        'interpreter', 'none');
    set(gca,...
        'xlim', [0 2],...
        'ylim', [0 2],...
        'ticklength', [0 0],...
        'xticklabel', {},...
        'yticklabel', {},...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'box', 'off');
    
    % main plot
    subplot('position', posMain, 'nextplot', 'add');
    x = normlinsize';          
    y = cpercent';  
    xline = [zeros(1, length(x)); x'];
    yline = [1; 1] * y';        
    plot(xline,          yline,   '-',  'linewidth', .75, 'color', col);
    plot(x,              y,       '-',  'linewidth', 3.5, 'color', col);        
    plot([1 1] * x(end), [0 100], '--', 'linewidth', 1.5, 'color', col);
    set(gca,...
        'ticklength', [0 0],...
        'box', 'off',...
        'fontname', 'arial',...
        'fontsize', fonttic,...
        'xlim', [0 maxx],...
        'ylim', [miny 100 + ((100 - miny) * .001)]);   
        
return
                      
    % function for run-length decoding (based on code by Peter J. Acklam)
    function out = rldecode(a, b)
        a = a(:)';
        b = b(:)';
        i = cumsum([1 b]);
        j = zeros(1, i(end)-1);
        j(i(1:end-1)) = 1;
        out = a(cumsum(j));
    return
    
    
    
    