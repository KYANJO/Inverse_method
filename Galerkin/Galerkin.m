function [A,B,ntime,x,dt] = Galerkin(D, Tfinal,ax,bx,dt)

    addpath('/Users/mathuser/Documents/PHD/Spring_2022/Inverse Theory/Project/Galerkin/Galerkin_functions')
    
    %Input Data
    nelem = 333; %Number of Elements == 22
    nop = 3;    %Interpolation Order == 5

    integration_type=1; %=1 is inexact and =2 is exact

    %Store Constants
    ngl = nop + 1;
    
    npoin = nop*nelem + 1;
    
    if(nargin < 5)
        cfl = 160;
        dx = (bx-ax)/(nelem);
        dt_est = cfl*dx;
        ntime = floor(Tfinal/dt_est);
        dt = Tfinal/ntime;
        
    else
    
        ntime = floor(Tfinal/dt);
        
    end
    
    %Compute Interpolation and Integration Points
    [xgl,wgl]=legendre_gauss_lobatto(ngl);

    if (integration_type == 1)
        noq=nop;
    elseif (integration_type == 2)
        noq=nop+1;
    end
    nq=noq + 1;
    [xnq,wnq]=legendre_gauss_lobatto(nq);


    %Compute Lagrange Polynomial and derivatives
    [psi,dpsi] = lagrange_basis3(ngl,nq,xgl,xnq);

    %Create Grid
    [coord,intma]=create_grid(ngl,nelem,xgl,-1,1);

    %Form Global Mass and Differentiation Matrices

    Mmatrix = create_Mmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi);

    Lmatrix = create_Lmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,dpsi);
    
    A = Mmatrix;
    B = Mmatrix  + D*dt*Lmatrix;
        
    A = sparse(A);
    B = sparse(B);
    %Generate Grid Points
    x = zeros(npoin,1);
    for e=1:nelem
        for i=1:ngl
            ip=intma(i,e);
            x(ip) = coord(i,e);
           
        end
    end

end