function dx = igseq_dx(db, x)

% IGSEQ_DX Number of lineages in top x% of lineage-size-ranked sequences
%
%    dx = IGSEQ_DX(dbname, x) returns the Dx value for each population and
%    isotype combination in the specified analysis unit.  The first
%    argument is a sqlite database name (including the full path) as a
%    string.  The second argument is a list of "x" values to use.  The
%    return value is a nxm cell array, where m = 3 + the number of "x"
%    values used.  The first column is the population.  The second column
%    is the isotype - no value indicates all isotypes were considered.  The
%    third column is the species richness, or in this context, the number 
%    of lineages (when filtering on a specific isotype, lineages were still
%    computed by the pipeline using all isotypes, but here we are only 
%    considering sequences of a particular isotype to extract lineages to
%    compute Dx values).  The remaining columns are the Dx values
%    corresponding to the list of "x" parameters specified.
%
%    For example, the D20 is the number of lineages, that, when ordering
%    the lineages by size, include the 20% of sampled sequences starting
%    from the largest lineage.  If a lineage straddles the 20% point, it is
%    included in the count. 
%        
%    Requires Matlab Database Toolbox
%
%    See also IGSEQ_PIELOU, IGSEQ_HILLPLOT, HILLDIV

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % open database and get populations, lineages IDs and isotypes
    conn = igseq_db(db);
    a = table2cell(get(fetch(exec(conn, 'select population, lineageID, isotype from sequence where vgene != "IGHV1-2"')), 'data'));
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
    
    % number of dx values
    nx = length(x);
    
    % initialize output cell array
    dx = cell(n, 3 + nx);
   
    % column 1 are populations
    dx(:, 1) = reshape(repmat(upop', (1 + niso), 1), (npop * (1 + niso)), 1);
    dx(:, 2) = repmat([{''}; uiso], npop, 1);
    
    % compute dx
    x = 100 - x;
    for i = 1:n
        if ~isempty(dx{i, 2})
            k = intersect(find(strcmp(dx{i, 2}, iso)), find(strcmp(dx{i, 1}, pop)));
        else
            k = find(strcmp(dx{i, 1}, pop));
        end
        if k > 0
            c = lin(k);
            t = tabulate(c);
            t = sort(t(t(:, 2) > 0, 2));
            pct = 100 * cumsum(t) / sum(t);
            dx{i, 3} = length(pct);
            for j = 1:nx
                dx{i, j + 3} = length(pct) - find(pct <= x(j), 1, 'last');
            end
        else
            dx{i, 3} = 0;
            for j = 1:nx
                dx{i, j + 3} = NaN;
            end
        end            
    end
    
    return
            
            
            
            
            
            
            
            
            
            