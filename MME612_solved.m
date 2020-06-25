%% Introduction
% This in-class activity was developed by Alex Mazursky and Paul Goetze. It
% intends to introduce Dr. Amit Shukla's MME 612 class to the adaptive
% gridding scheme we applied to our research project.

% TO BE FILLED IN: 
% SEE LINES 103-115.

%% Problem Initialization
clear all
close all
clc
% Building spaces
x_l = 12; % Length of Plate in x
dx = 1.0; % Step size of grid in x
y_l = 12; % Length of Plate in x
dy = 1.0; % Step size of grid in x
x_space = 0:dx:x_l; % assemble space in x
y_space = 0:dy:y_l; % assemble space in y
Lx = length(x_space); % Counts number of values in space
Ly = length(y_space);

S = zeros(Lx,Ly); % preallocate S, the coarse solution.

% Build Boundary Conditions for entire problem
% where m is in x and n is in y
for m = 1:Lx
    for n = 1:Ly
        if m == 1 % Left side boundary
          S(n,m) = 600-50*n;
        elseif m == Lx % Right side boundary
            S(n,m) = 100;
        elseif n == 1 % Bottom boundary
            S(n,m) = 600-50*m;
        elseif n == Ly % Top boundary
            S(n,m) = 100;
        end
    end
end
tmax = max(max(S)); % Extract maximum temperature
% used for scaling the plot's colors.

%% Gauss-Seidel Solver
pmax = 10^7; % maximum iterations
tol = 10^-1; % desired tolerance to exit iterative solver

sdivideX = floor(Lx/2);
sdivideY = floor(Ly/2);

cThub = .5;
cPcb = .8;

for p = 1:pmax
    % solution loop
    for n = 2:Ly-1
        for m = 2:Lx-1
            coeff = 1;
        % Governing heat conduction equation
               if (n < sdivideY && m < sdivideX)
                   coeff = cPcb;  
               else
                  coeff = cThub;
              end
            S(m,n) = coeff*0.25*(S(m+1,n)+S(m-1,n)+S(m,n+1)+S(m,n-1));
            fprintf('S: %d\n',S(m,n));
        end
    end
    
    % residual loop
    for n = 2:Ly-1
        for m = 2:Lx-1
        % calculate the residual for all nodes
        r(m,n) = abs(0.25*(S(m+1,n)+S(m-1,n)+S(m,n+1)+S(m,n-1))-S(m,n));
        fprintf('r: %d\n',r(m,n));
        end
    end
    % Normalize the residual for the entire plate
    Res = sum(sum(r(2:Lx-1,2:Ly-1)))/sum(sum(S(2:Lx-1,2:Ly-1)));
    
    fprintf('Res = %e Iter = %8.0d - Res = %4.2e \n',Res, p, Res);
    if (Res < tol)
    break % Exit iteration loop because of convergence
    end
    if (p == pmax)
    disp('Warning: code did not converge')
    % Hopefully this doesn't happen!
    end
end

%% Plot the coarse solution
figure('Name','Coarse Solution','NumberTitle','off');
surface(x_space,y_space,S)
% shading interp
colormap('Jet')
c = colorbar;
c.Limits = [0 tmax];

%% Refinement of the Solution
% We see the concentrated temperature in the bottom left corner. 
% This is our point of interest for refining the solution.
% We use variables ending in '1' to denote that they are particular to the
% refinement.

% Rebuild the space about a smaller region.
x_l1 = 4;
dx1 = dx/2;
y_l1 = 4;
dy1 = dy/2;
x_space1 = 0:dx1:x_l1;
y_space1 = 0:dy1:y_l1;
Lx1 = length(x_space1);
Ly1 = length(y_space1);

S1 = zeros(Lx1,Ly1);
% Build Boundary Conditions
% Pulling from Solution S...
% Make the boundary conditions match.
%-----------------------------------------------------------
for m = 1:Lx1
    for n = 1:Ly1
        if m == 1  
            S1(n,m) = S(ceil(n*dy1),ceil(m*dx1));
        elseif m == Lx1 && n >= 2
            S1(n,m) = S(ceil(n*dy1)+1,ceil(m*dx1));
        elseif n == 1
            S1(n,m) = S(ceil(n*dy1),ceil(m*dx1));
        elseif n == Ly1 && n>=2
            S1(n,m) = S(ceil(n*dy1),ceil(m*dx1)+1);
        end
    end
end 
%-----------------------------------------------------------
% Having established the new, finer space and boundary conditions, use
% another Gauss-Seidel formulation to solve.
for p = 1:pmax % gauss-seidel
    % solution loop
    for n = 2:Ly1-1
        for m = 2:Lx1-1
        S1(m,n) = 0.25*(S1(m+1,n)+S1(m-1,n)+S1(m,n+1)+S1(m,n-1));
        end
    end
    % residual loop
    for n = 2:Ly1-1
        for m = 2:Lx1-1
        r1(m,n) = abs(0.25*(S1(m+1,n)+S1(m-1,n)+S1(m,n+1)+S1(m,n-1))-S1(m,n));
        end
    end
    % Normalized Residual
    Res1 = sum(sum(r1(2:Lx1-1,2:Ly1-1)))/sum(sum(S1(2:Lx1-1,2:Ly1-1)));
    fprintf('Iter = %8.0d - Res1 = %4.2e \n', p, Res1);
    if (Res1 < tol)
    break % Exit iteration loop because of convergence
    end
    if (p == pmax)
    disp('Warning: code did not converge')
    end
end

%% Plot the fine solution and overlay the two.
figure('Name','Refined Solution','NumberTitle','off');
surface(x_space1,y_space1,S1)

% shading interp
colormap('Jet')
c = colorbar;
c.Limits = [0 tmax];

hold on

figure('Name','Full Solution','NumberTitle','off');

surface(x_space(x_l1+1/dx:length(x_space)),y_space(x_l1+1/dx:length(y_space)),S(x_l1+1/dx:length(x_space),x_l1+1/dx:length(y_space)))
surface(x_space(1:x_l1+1),y_space(x_l1+1:length(y_space)),S(x_l1+1:length(y_space),1:x_l1+1))
surface(x_space(x_l1+1:length(x_space)),y_space(1:x_l1+1),S(1:x_l1+1,x_l1+1:length(x_space)))
surface(x_space1(1:length(x_space1)),y_space1(1:length(y_space1)),S1(1:length(y_space1),1:length(x_space1)))

% shading interp
colormap('Jet')
c = colorbar;
c.Limits = [0 tmax];

% Postprocessing and Beautification

title({'T $(^{\circ}$C)'},...
'FontUnits','points',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,...
'FontName','Times')

ylabel({'y-position'},...
'FontUnits','points',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,...
'FontName','Times')

xlabel('x-position',...
'FontUnits','points',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,...
'FontName','Times')

set(gca, ...
  'FontName','Times',...
  'FontSize'    , 16, ...
  'Box'         , 'off'     , ...
  'Ytick', [0 0.5 1 1.5 2],...
  'LineWidth'   , 1         );