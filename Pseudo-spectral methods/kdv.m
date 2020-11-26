%KDV  Solve the Korteweg-de Vries equation
%   u_t + alpha * u u_x + u_{xxx} = 0
%
% This script implements a numerical scheme for solving the Korteweg-de
% Vries equation. The script is a complete solver of the KdV equation in
% that it starts with setting up a uniform mesh on which the numerical
% solution will be computed. The periodic spectral method, which is
% efficiently implemented via the Fast Fourier Transform, will be used for
% the spatial discretization. After such discretization, we obtain a system
% of ordinary differential equations (ODEs) for the solution at each node
% (a degree of freedom, or a DoF). This system can be written in the
% (discrete) Fourier space or in the physical space and can be solved by
% various time integration schemes. Among them, the 4th-order Runge-Kutta
% method is an ideal choice due to several reasons. (i) It is an explicit
% scheme with high accuracy, the computational time is acceptable for a
% large simulation time. (ii) The scheme is available as a built-in
% function in MATLAB. In fact, the algorithm is highly optimized that a
% naive implementation of a user cannot compete with this package. That
% said, one can also choose other time integration schemes that are also
% avaialable as built-in functions in MATLAB.
%
% The solution will be propagated in the Fourier space. That is, we solve
% the discretized equation in the Fourier space because it produces a
% better numerically stable procedure. The built-in ode45 takes care of the
% time step under the hood to achieve the stable solution or throws the
% warning and stops if the solution is not numerical stable. Thus, a stable
% numerical procedure also means less computational resources.
%
% To improve numerical stability, we adopt the strategy outlined in [1].
% Specifically, we eliminate the dispersive term u_{xxx} in the Fourier
% space by multiplying the entire equation (in Fourier space) by the
% integrating factor exp(-i*k^3*t) and applying the transformation 
%   fft(v) = exp(-i *k^3 * t) fft(u)
% The new unknown would be fft(v) = \hat{v} and the equation must be
% slightly modified (see [1], Chapter 10, page 111). The solution is
% finally obtained by the inverse transformation
%   u = ifft( exp(i * k^3 * t) * fft(v) ).
%
% [1] <a href="http://www.audentia-gestion.fr/Matlab/10.1.1.473.7647.pdf"
% > Spectral Methods in MATLAB </a>
% [2] <a href="https://www.math.upenn.edu/~jucurry/papers/soliton.pdf"
% > Exact one-soliton solution of the KdV equation </a>
% 
% Remark: For high efficiency, the number of grid points should be an even
% number of the power of 2, that is 2^N. The KdV equation will be solved on
% the interlval [a, b].
%
% See also BO

%% Space-time mesh
a = -20;
b = 120;
N = 2^10;
L = b - a;
scale = 0.5*L/pi;
h = L/N;
if mod(N, 2) == 0
    x = (a : h : b-h)';
    k = [0:N/2-1, 0, -N/2+1 : -1]'/scale;
else
    x = (a+h/2 : h : b-h/2)';
    k = [0:(N-1)/2, -(N-1)/2: -1]'/scale;
end
tmax = 50;
alpha = 6;      % typical form of the KdV equation.

% Initial condition: We use one-soliton solution to generate the initial
% condition. One can use any arbitrary (sufficiently regular) functions for
% the initial condition.
c = 1;          % velocity of solitary wave
iu = 0.5*c * sech(sqrt(c)/2 * (x - 0)).^2;
%% Time integration using pseudo-spectral method
tstart = tic;
%--------------------------------------
ik = 1i*k;
q = -1i * k.^3;
odefunc = @(t, vFourier) ...
    -0.5*alpha*ik.*exp(q*t) .* fft( ifft(exp(-q*t).*vFourier).^2  );
tspan = linspace(0, tmax, 101);
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-5);
[t, vFourier] = ode45(odefunc, tspan, fft(iu), options);
u = ifft(exp(-t * (q.')) .* vFourier, [], 2);
%--------------------------------------
telapsed = toc(tstart); fprintf('Time processed: %f \n', telapsed);

%% Plot the solution
close all
[xx, tt] = meshgrid(x, t);
mesh(xx, tt, u); view(3)

%% Time integration using pseudo-spectral method
% The following scheme is much slower than and less numerically stable than
% the above numerical scheme. This numerical scheme is based on a naive
% implementation of the 4th-RK method. The solution is sought as the direct
% Fourier transform of u. Thus, the master function @func(t,u) for the
% built-in ode45 is constructed in terms of fft(u). The condition for
% numerical stability is much stricter than the above numerical schem.
% Since the package ode45 takes care of the time step to achieve stable
% solution (or giving warning and stopping), this method takes much more
% time to execute. The corresponding code is given here for the interest of
% the educators and the students. One should use the above instead

% Uncomment to try the code.
%--------------------------------------

% tstart = tic;
% ik = 1i*k;
% ik3 = 1i * k.^3;
% odefunc = @(t, UFourier) ik3.*UFourier - 0.5*alpha*ik.*fft( ifft(UFourier).^2 );
% tspan = linspace(0, tmax, 101);    % time at which the solution is stored.
% options = odeset('RelTol', 1e-3, 'AbsTol', 1e-5);
% % Solve the wave equation in Fourier space using ode45 in MATLAB
% [t, uFourier] = ode45(odefunc, tspan, fft(iu), options);
% telapsed = toc(tstart); fprintf('Time processed: %f', telapsed);
