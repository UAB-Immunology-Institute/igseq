function [d1, q1] = hilldiv(composition, orders)

% HILLDIV Hill diversity index and several orders.
%
%    [d, q] = HILLDIV(comp, orders) returns the Hill diversity indices for
%    the supplied composition (first argument) at all of the specified
%    orders (second argument).  The first return value is a vector of
%    diversity values associated with the orders (second return value).

% Chris Fucile
% Alex Rosenberg
% 19 Feb 2016
% University of Rochester, Div. of Allergy, Immunology and Rheumatology
% Copyright (c) 2016. All rights reserved.
% This software is offered with no guarantees of any kind.

    % make sure composition truly adds up to 100%
    p = composition(:) / sum(composition);

    % deal with non-special cases first
    q = setdiff(orders, [0 1 Inf]);

    % compute Hill diversity at all specified orders
    qr = repmat(q(:)', length(p), 1);
    p = repmat(p, 1, length(q));
    s = sum(p .^ qr);
    d = s .^ (1 ./ (1 - q));
    q = q(:);
    d = d(:);

    % order zero (number of items in composition)
    if ~isempty(find(orders == 0, 1))
        d = [d; length(composition)];
        q = [q; 0];
    end

    % special case for order 1 
    if ~isempty(find(orders == 1, 1))
        d = [d; 1 / prod(p(:, 1) .^ p(:, 1))];
        q = [q; 1];
    end

    % special case for order infinity
    if ~isempty(find(isinf(orders), 1))
        d = [d; 1 / max(p(:, 1))];
        q = [q; Inf];
    end

    % re-sort so that special cases in right order
    [~, isort] = sort(q);
    q1 = q(isort);
    d1 = d(isort);

    return