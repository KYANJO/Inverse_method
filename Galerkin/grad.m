
function g = grad(u0)

    % Global variables.

    global m;
    global n;
    global xpoints;
    global Fvec;
    global alpha;
    global L;
    
    % Now, use the discrete adjoint equation to find the gradient.
    
    v = zeros(n,1);
    v(xpoints) = 2*Fvec;
    g = adjoint(m,v);
    
    % Add the gradient of the regularization term.
    
    g = g+2*alpha*L'*(L*u0);
    
end