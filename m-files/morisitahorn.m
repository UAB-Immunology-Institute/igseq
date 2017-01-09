function o = morisitahorn(r)

% MORISITAHORN Morisita-Horn overlap index.
%
%    o = MORISITAHORN(reps) returns a square matrix containing the
%    Moritsita-Horn overlap indices for all possible pairs of repertoire
%    comparisons.  Each column of the input matrix represents a repertoire
%    and each row represents a species (e.g. a grouped lineage).  The
%    values of the matrix represent the number of individuals (e.g.
%    sequences) for the corresponding species and repertoire.  For a
%    pairwise repertoire comparison, the index varies from 0 (no overlap)
%    to 1 (identical species and corresponding frequencies).

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % number of species and number of repertoires
    [~, nr] = size(r);

    % initialize output
    o = nan(nr);

    % loop through all pairs of repertoires
    for j = 1:nr
        for i = 1:j
            numerator = 2 * sum(prod(r(:, [i j]), 2));
            X = sum(r(:, i));
            Y = sum(r(:, j));
            sx = sum(r(:, i) .^ 2) / (X ^ 2);
            sy = sum(r(:, j) .^ 2) / (Y ^ 2);
            denominator = (sx + sy) * X * Y;
            o(i, j) = numerator / denominator;
        end
    end

    % make symmetric
    o(isnan(o)) = 0;
    o = triu(o, 1)' + o;
        
    return


