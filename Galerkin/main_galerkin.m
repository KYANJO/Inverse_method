
clear all
clc

% Global variables.

global m;
global n;
global deltax;
global x;
global deltat;
global xpoints;
global d;
global D;
global Fvec;
global alpha;
global L;

global A;
global B;

% load data

load('data.mat')

% Setup model discretization

ax=-1;
bx=1;
deltax=0.002;
xx = (deltax:deltax:(1000*deltax))'-1;
xpoints1 = (25:25:975)';
xd1 = xx(xpoints1);

% Set the thermal diffusion coefficient (m^2/s).

D=1e-6;

% Set the time period T s.

T=1e4;

% Compute matrix A and B
[A,B,m,x,deltat] = Galerkin(D,T,ax,bx);

% Initial condition (true solution, -1<=x<=1).

n = length(x);
mdlx = floor(n/2);
p60 = floor(n*60/100);
p20 = floor(n*20/100);
p25 = floor(n*25/100);
p77 = floor(n*77.5/100);
p80 = floor(n*80/100);

u0true = zeros(size(x));


u0true(mdlx:p60) = 10;
u0true(p20:p25) = 5;
u0true(p77:p80) = 3;

% Compute the solution.

uT = forward(m,u0true);

% Get the data points.

xpoints = zeros(39,1);
for i = 1:39
    j = find(x<=xd1(i));
    xpoints(i) = j(end);
end

% Plot the solution at time T and the noisy data points. 

figure(1);
clf
plot(x,uT,'k');
hold on
plot(x(xpoints),d,'ko','LineWidth',1);
xlabel('x (m)')
ylabel('\Delta T (^oK)','LineWidth',1);
title('Galerkin method')

% Now, setup and solve the inverse problem. 
% First, setup the regularization

% %zeroth-order Tikhonov matrix for the objective function
L = speye(n);

%search over a range of regularization parameters
alphas = logspace(-5,-2.5,20);

[resid, mnorm] = deal(zeros(length(alphas),1));
u0save = zeros(n,1);

for i = 1:length(alphas)
    
    alpha = alphas(i);

    % Setup a first guess for u0(a uniform temperature perturbation)
    
    u0guess = 10*ones(n,1);
    
    % Call the conjugate gradient method to find the best solution.
    
    [u0min,fu0min] = conjg(@objfun,@grad,u0guess,1e-10);
    
    
    uTmin = forward(m,u0min);
    
    resid(i) = norm(d-uTmin(xpoints));
    disp(['Norm of Residual vector for alpha',num2str(i),' =',num2str(alphas(i)),' norm = ',num2str(resid(i))])
    u0save(:,i) = u0min;
    mnorm(i) = norm(u0min);

end

save('uoData1.mat', 'u0save');
save('mnormData1.mat', 'mnorm');
save('residData1.mat', 'mnorm');

noiselevel = 0.1;
discrep = sqrt(length(xpoints))*noiselevel;
disp(['Discrepancy Principle Residual Criterion: ',num2str(discrep)])

figure(2)
idx = find(resid<discrep);
id = idx(end);
alpha_optimal = alphas(id);
ming = resid(id);

loglog(alphas,resid,'LineWidth',2)

ax=axis; hold on
H=loglog(alpha_optimal,ming,'.r',[alpha_optimal,alpha_optimal],[ming/length(alphas),ming],'--r');
%text(alpha_optimal,0.1,num2str(alphas(id)),'fontsize',20);
set(H,'markersize',15);
hold off
axis(ax)
xlabel('\alpha');
ylabel('Residual Norm ||u(u_0) - d||_2')
title('Galerkin method')
 
% Plot the closest model to meeting the discrepancy principle

figure(3);
clf
plot(x,u0true,'k','LineWidth',1);
hold on
plot(x,u0save(:,id),'k--','LineWidth',1);
legend('u0 true',['recovered model for \alpha = ',num2str(alphas(id))]);
title('Galerkin method')
xlabel('x (m)')
ylabel('\Delta T (^oK)');
ylim([-2 13])


