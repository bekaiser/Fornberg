% Derivative tests for Fornberg (1998) Lagrange polynomial derivatives
% Bryan Kaiser
% 12/20/15

close all
clear all
clc

% =============================================================================
% domain parameters

Lx = 3000.0; % km, domain size
Lxcenter = 0.0; % x value @ the center of the grid
N = 2^16; % series length (must be at least even)
dx = Lx/N; % km, uniform longitudinal grid spacing
x = [0.5*dx:dx:dx*N]-(Lx/2.0-Lxcenter); % km, centered uniform grid

% =============================================================================
% signal

k = 2.0*pi/(Lx/5.0);
s = sin(x.*k);
dsdx = cos(x.*k).*k;
d2sdx2 = -sin(x.*k).*k^2;

figure(1)
plot(x,s,'b',x,dsdx,'r',x,d2sdx2,'g')
xlabel('x')
title('signal')


% =============================================================================
% compute the uniform grid weights

c = weights(x(3),x(1:5),4,2)'; % 5 point stencil weights
dsF = zeros(1,N); % derivative of the input function
m = 1; % order of the derivative

dsF = zeros(1,N);
for j = 3:N-2      
        jm2 = j-2;
        jp2 = j+2;
        dsF(1,j) = c(m+1,:)*s(jm2:jp2)';
end
error = abs(dsdx-dsF);

figure(2)
semilogy(x(3:N-2),error(3:N-2),'k')
xlabel('x')
title('Lagrange polynomial 1st derivative error')

