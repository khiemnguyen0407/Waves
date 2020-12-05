%BENJAMIN  Solve the Benjamin equation
%   
%
% BENJAMIN implements a numerical scheme for solving the Benjamin equation
%
%   u_t + u_x + 2*u*u_x - alpha*Hu_{xx} - beta * u_{xxx} = 0.   (1)
%
% As the name suggests, the above equation is derived by Benjamin in [1].
% The dispersive term of the Benjamin equation contains the dispersion of
% the Benjamin-Ono type and the KdV type. Particularly, the first
% dispersive term, namely Hu_{xx}, corresponds to the BO equation and the
% second, namely u_{xxx}, to the KdV equation. 
%
% The script is a complete solver of the Benjamin equation. That is, it
% starts with setting up a uniform mesh on which the numerical solution
% will be computed, discretize the solution in space and solve the resulted
% system of ODEs in time. Unlike the solver for KdV equation (see kdv.m),
% this equation cannot be solved by other methods such as the finite
% difference and finite element methods due to the appearance of the
% Hilbert transform.  The Hilbert transform of function f = f(x) is defined
% as the convolution of a singular kernel with the function itself. Thus,
% although it is possible to approximate a Hilbert transform in a numerical
% basis, it is not trivial to integrate such approximations into a
% numerical scheme for the wave solution. For more explanation on the
% numerical strategy, type "help kdv".
%
% Equation (1) is a prototype of non-integrable wave equation and thus
% there is no exact solution of this equation. For this reason, we use the
% jump initial condition that is typically used to generate the dispersive
% shock waves. The jump condition is smoothed out by using the
% regularized tanh function instead of the straight constant lines.

%% Space-time mesh
a = -150;       % left boundary
b = 150;        % right boundary
L = b - a; 
N = 2^10;     % number of modes
h = L/N;
scale = L/(2*pi);
if mod(N, 2) == 0   % define collocation & wave numbers
    x = (a : h : b-h)';
    k = [0:N/2-1, 0, -N/2+1 : -1]'/scale;
else
    x = (a+h/2 : h : b-h/2)';
    k = [0:(N-1)/2, -(N-1)/2: -1]'/scale;
end
tmax = 20;

alpha = +1; % Parameters for Benjamin equation.
beta  = +1;

% Initial condition: We use the jump condition for generating DSWs.
jumpLocation = [a, 0];
uLimits = [1, 0];
delta = 1.5*ones(1,2);
iu = 0.5*(uLimits(1) - uLimits(2))*(tanh(delta(1)*(x-jumpLocation(1))) ...
        - tanh(delta(2)*(x-jumpLocation(2)))) + uLimits(2);

%% Time integration in the Fourier space
tstart = tic;

ik = 1i*k;
q = -(1i*k + alpha*1i*k.*abs(k) + 1i*beta*k.^3);
odefunc = @(t, u_fourier) q.*u_fourier - ik.*fft(ifft(u_fourier).^2);
tspan = linspace(0, tmax, 51);
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

[t, sol_fourier] = ode45(odefunc, tspan, fft(iu));
u = ifft(sol_fourier, [], 2);

telapsed = toc(tstart); fprintf('Time processed: %f \n', telapsed);
%% Plot the solution.
plot(x, u(1,:), 'b-', x, u(end,:), 'k-', 'LineWidth', 1.5);
