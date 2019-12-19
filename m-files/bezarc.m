function [x, y] = bezarc(x1, y1a, y1b, x2, y2a, y2b, offset)


% offset is a number between 0 and .5


%This function constructs a Bezier points using n input control points. 
%
% - Sagar Aiya (f2008149@bits-goa.ac.in)

    PA(1, :) = [x1, x1 + ((x2 - x1) * offset), x2 - ((x2 - x1) * offset), x2];
    PB(1, :) = PA(1, :);

    PA(2, :) = [y1a y1a y2a y2a];
    PB(2, :) = [y1b y1b y2b y2b];

    PA(3, :) = 0;
    PB(3, :) = 0;



    cbot = getpoly(PA);
    ctop = getpoly(PB);

    x = [cbot(1, :) fliplr(ctop(1, :)) cbot(1, 1)];
    y = [cbot(2, :) fliplr(ctop(2, :)) cbot(2, 1)];

return


    function A = getpoly(P)

        n = 4;
        count = 1;

%         div = 50; % # of segments of curve; increase for smoother curve
        div = 200; % # of segments of curve; increase for smoother curve

        for u = 0:(1/div):1
            sum = [0 0 0]';
            for i = 1:n
                B = nchoosek(n,i-1)*(u^(i-1))*((1-u)^(n-i+1)); %B is the Bernstein polynomial value
                sum = sum + B*P(:,i);
            end
            B = nchoosek(n,n)*(u^(n));
            sum = sum + B*P(:,n);
            A(:,count) = sum; %the matrix containing the points of curve as column vectors. 
            count = count+1;  % count is the index of the points on the curve.
        end
        
    return

    
   