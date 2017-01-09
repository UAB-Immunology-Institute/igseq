function ev = igseq_pielou(db)

% IGSEQ_PIELOU Species evenness for all population/isotype combinations
%
%    ev = IGSEQ_PIELOU(dbname) returns the species evenness value for each
%    population and isotype combination in the specified analysis unit. 
%    The first argument is a sqlite database name (including the full path)
%    as a string.  The return value is a nx4 cell array.  The first column
%    is the population.  The second column is the isotype - no value
%    indicates all isotypes were considered.  The third column is the
%    species richness, or in this context, the number of lineages (when
%    filtering on a specific isotype, lineages were still computed by the
%    pipeline using all isotypes, but here we are only considering
%    sequences of a particular isotype to extract lineages to compute
%    evenness).  The fourth column is species evenness: 1 is uniform
%    (completely polyclonal), and values close to 0 are highly oligoclonal.
%        
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_DX, IGSEQ_HILLPLOT, HILLDIV

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % open database and get populations, lineages IDs and isotypes
    conn = igseq_db(db);
    a = get(fetch(exec(conn, 'select population, lineageID, isotype from sequence')), 'data');
    close(conn);

    % get columns from query result
    pop = a(:, 1);
    lin = cell2mat(a(:, 2));
    iso = a(:, 3);
    
    % unique populations and isotypes
    upop = unique(pop);
    uiso = unique(iso);
    
    % eliminate "U" from isotypes
    uiso = setdiff(uiso, 'U');
    
    % numbers of populations and isotypes
    npop = length(upop);
    niso = length(uiso);
    
    % compute number of rows output
    n = npop * (1 + niso);
    
    % initialize output cell array
    ev = cell(n, 4);
   
    % column 1 are populations, column 2 are isotypes
    ev(:, 1) = reshape(repmat(upop', (1 + niso), 1), (npop * (1 + niso)), 1);
    ev(:, 2) = repmat([{''}; uiso], npop, 1);
    
    % compute evenness
    for i = 1:n
        if ~isempty(ev{i, 2})
            k = intersect(find(strcmp(ev{i, 2}, iso)), find(strcmp(ev{i, 1}, pop)));
        else
            k = find(strcmp(ev{i, 1}, pop));
        end
        if k > 0
            c = lin(k);
            t = tabulate(c);
            t = sort(t(t(:, 2) > 0, 2), 'descend');        
            p = t ./ sum(t);      
            H = -1 * sum(p .* log(p));
            S = length(p);
            Hmax = log(S);
            J = H / Hmax;
            ev{i, 3} = S;
            ev{i, 4} = J;
        else
            ev{i, 3} = 0;
            ev{i, 4} = NaN;
        end
    end
    
    return