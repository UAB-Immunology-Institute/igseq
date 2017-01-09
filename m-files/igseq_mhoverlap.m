function t = igseq_mhoverlap(db)

% IGSEQ_MHOVERLAP compute Morisita-Horn overlaps in an analysis unit
%
%    t = IGSEQ_MHOVERLAP(dbname) returns a square matrix containing the
%    Morisita-Horn overlap indices for all possible pairs of repertoires in
%    the specified database (corresponding to an analysis unit.  The output
%    matrix is a cell array with the first row and column indicating the
%    population names.  For a pairwise repertoire comparison, the index
%    varies from 0 (no overlap) to 1 (identical species and corresponding
%    frequencies).  The input argument is a sqlite database name (including
%    the full path) as a string.
%
%    Requires Matlab Database Toolbox
%
%    See also MORISITAHORN

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % get grouped lineage table
    g = igseq_groupedlin(db);

    % compute morisita-horn indices for each repertoire pair
    o = morisitahorn(cell2mat(g(2:end, 4:end)));

    % extract population names
    pop = g(1, 4:end);

    % prepare output cell array
    t = cell(1 + size(o));
    t(2:end, 2:end) = num2cell(o);
    t(1, 1) = {''};
    t(1, 2:end) = pop;
    t(2:end, 1) = pop';

return