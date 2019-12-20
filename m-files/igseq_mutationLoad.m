function igseq_mutationLoad(db,t,pops)

% IGSEQ_MUTATIONLOAD mutation frequency histograms for specific populations.
%
%    h = IGSEQ_MUTATIONLOAD(dbname, type, population) returns a figure of
%    binned mutation frequency histograms for the given populations
%    within an analysis unit (database).  The first argument
%    is a sqlite database name (including the full path) as a string.  The
%    second argument is one of: "total", "silent" or "nonsilent" to specify
%    the type of mutations to consider.  The third argument is a cell array
%    specifying the populations of interest. 
%
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_ISOTYPES, IGSEQ_MUTHISTPLOT

% Chris Fucile
% Alex Rosenberg
% 04 Oct 2017
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    switch t
        case      'total',     f = 'V_Mutations';
        case      'silent',    f = 'V_SilentMutations';
        case      'nonsilent', f = 'V_NonsilentMutations';
    otherwise
        error('type must be one of: ''total'', ''silent'', ''nonsilent''');
    end

    pop_query = join(pops,''', ''');
    pop_query = char(strcat('(''',pop_query,''')'));

    q = ['select lineageID, population, V_Nucleotides, ' f ' from sequence where population in ' pop_query];
    conn = igseq_db(db);
    a = get(fetch(exec(conn, q)), 'data');
    close(conn);


    tmp = strtrim(cellstr(num2str(cell2mat(a(:,1)))));
    tmp = strcat(tmp(:,1),':',a(:,2));

    a(:,5) = tmp;

    % compute frequencies
    frq  = 100 * cell2mat(a(:, 4)) ./ cell2mat(a(:, 3));
    upop = unique(a(:,5));
    npop = length(upop);
    for i = 1:npop
        k = find(strcmp(a(:, 5), upop{i}));
        m(i,1) = upop(i);
        m(i,2) = num2cell(mean(frq(k)));
    end

    p = tabulate(a(:,5));
    p2 = tabulate(a(:,2));

    for i = 1:length(p)
        k = find(strcmp(m(:,1),p(i,1)));
        pop = strsplit(p{i,1},':');
        j = find(strcmp(p2(:,1),pop(1,2)));
        m(k,3) = num2cell(p{i,2} / p2{j,2} * 100);
    end

    for i = 1:length(m)
        pop = strsplit(m{i,1},':');
        m(i,4) = pop(2);
    end

    colorOrder = get(gca, 'ColorOrder');


    handles = [];
    for i = 1:length(pops)
        data = find(strcmp(m(:,4),pops{i}));
        col = colorOrder(mod(length(get(gca, 'Children')), size(colorOrder, 1))+1, :);
        handles(i) = plot(log10(cell2mat(m(data,3))),cell2mat(m(data,2)),'o ','Color',col);
        fig = gcf;
        hold on;

        binCenter = [];
        data_mean = [];
        x = log10(cell2mat(m(data, 3)));
        y = cell2mat(m(data, 2));
        binNum = round((x - .25) * 2) + 7;
        bins = unique(binNum)';
        for i = 1:length(bins)
            binCenter(i) = (bins(i) / 2) - 3.25;
            data_mean(i) = median(y(binNum==bins(i)));
        end
        data_mean_plot = plot(binCenter, data_mean,'Linewidth',2,'Color', col);
        hold on;
    end
    legend(handles, pops,'Interpreter', 'none');

    r = refline([0 5]);
    r.Color = [.5 .5 .5];
    r2 = refline([0 10]);
    r2.Color = [.5 .5 .5];

    y1 = get(gca,'ylim');
    plot([0 0], y1,'Color',[.5 .5 .5]);

end