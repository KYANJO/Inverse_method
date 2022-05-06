
% ssq=objfun(u0);
%
% Computes the sum of squares of the differences between u(x,T) and
% the data points.

function ssq=objfun(u0)

    % Global variables.

    global m;
    global xpoints;
    global d;
    global Fvec;
    global alpha;
    global L;

    % sole the forward model.

    ufinal = forward(m,u0);

    % Compute the differences

    Fvec = ufinal(xpoints)-d;

    % Compute the sum of squares and add the regularization term.

    ssq = norm(Fvec)^2 + alpha*norm(L*u0)^2;

end