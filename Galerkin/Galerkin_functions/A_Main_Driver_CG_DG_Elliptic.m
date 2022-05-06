%---------------------------------------------------------------------%
%This code computes the 1D Diffusion Equation using the 
%CG method with 3rd Order RK.
%This version constructs the Global Matrices

%---------------------------------------------------------------------%
clear all; 
close all;
clc
tic

%Input Data
nelem=32; %Number of Elements == 22
nop=4;    %Interpolation Order == 5

iplot_solution=1; %Switch to Plot or Not.

integration_type=1; %=1 is inexact and =2 is exact

icase=1; %case number: 1=Homogeneous Dirichlet with Sine;
                      %2=Homogeneous Dirichlet with Half-Cosine;
                      %3=Non-Homogeneous Dirichlet with Cosine
                      
D = 1;
ax = 0;
bx = 2*pi;
Tfinal = 1.0;
%Store Constants
ngl=nop + 1;

npoin=nop*nelem + 1;
cfl = 0.01;
dx = (bx-ax)/npoin;
dt_est = cfl*dx^2;
ntime = floor(Tfinal/dt_est);
dt = Tfinal/ntime;


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
[coord,intma_cg]=create_grid(ngl,nelem,xgl,ax,bx);

%Form Global Matrix Pointer

intma=intma_cg;

%Create Local/Element Mass and Differentiation Matrices
Mmatrix = create_Mmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi);

Lmatrix = create_Lmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,dpsi);

%Apply Dirichlet BC

Lmatrix(1,:)=0;
Lmatrix(1,1)=1;
Lmatrix(npoin,:)=0;
Lmatrix(npoin,npoin)=1;

Mmatrix(1,:)=0;
Mmatrix(1,1)=1;
Mmatrix(npoin,:)=0;
Mmatrix(npoin,npoin)=1;

% Initial solution
time = 0;
qe = exact_solution(coord,nelem,ngl,intma,npoin,icase,time);

Rmatrix = Mmatrix  + D*dt*Lmatrix;
q0 = qe;

for n = 1:ntime
    
    time = time + dt;
    qp = Mmatrix\(Rmatrix*q0);
    
    qe = exact_solution(coord,nelem,ngl,intma,npoin,icase,time);
    
    qp(1) = qe(1);
    qp(end) = qp(end);
    q0 = qp;
end


%Compute a gridpoint solution
x_sol=zeros(npoin,1);
for ie=1:nelem
   for i=1:ngl
      ip=intma(i,ie);
      x_sol(ip)=coord(i,ie);
   end 
end

%Plot Solution
if (iplot_solution == 1)
    h=figure;
    figure(h);
    plot_handle=plot(x_sol,qp,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(x_sol,qe,'b--');
    set(plot_handle,'LineWidth',2);

    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);

    title_text=[' Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    
end
