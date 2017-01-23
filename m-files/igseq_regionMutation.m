function m = igseq_regionMutation(db, t)

% IGSEQ_regionMutation mutation frequency histograms for an analysis unit.
%
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
        case      'total',     f = '_Mutations';
        case      'silent',    f = '_SilentMutations';
        case      'nonsilent', f = '_NonsilentMutations';
        otherwise
            error('type must be one of: ''total'', ''silent'', ''nonsilent''');
    end
    
    l = {'FR1','CDR1','FR2','CDR2','FR3'};
    r_mean = cell(6,length(l));
    sd = r_mean;

    for j = 1:length(l)
        % open database and run query to get all data
        q = strcat('select population,', ' ', l{j}, '_Nucleotides,', ' ', l{j}, f, ' from sequence');
        conn = igseq_db(db);
        a = get(fetch(exec(conn, q)), 'data');
        close(conn);

        % compute frequencies
        frq  = 100 * cell2mat(a(:, 3)) ./ cell2mat(a(:, 2));
        frq(isnan(frq)) = 0;
        upop = unique(a(:,1));
        npop = length(upop);
        
        for i = 1:npop
            k = find(strcmp(a(:, 1), upop{i}));
            upop{i};
            r_mean{i,j} = mean(frq(k));
            sd{i,j} = std(frq(k));
        end
    end
    r_mean
    sd