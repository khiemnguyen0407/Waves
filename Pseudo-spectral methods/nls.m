%NLS  Solve nonlinear Schroedinger equation
%
% The nonlinear Schroedinger equation in the dimensionaless form reads
% 
%   1i * \partial_{t} u = -\frac{1}{2} \partial_{x}^2 u + |u|^2 u,
% 
% where 1i stands for the imaginary unit, that is \sqrt{-1} = 1i. The
% solution u is complex-valued.
%
% In theoretical physics, the (one-dimensional) nonlinear Schrödinger
% equation (NLSE) is a nonlinear variation of the Schrödinger equation. It
% is a classical field equation whose principal applications are to the
% propagation of light in nonlinear optical fibers and planar waveguides
% and to Bose–Einstein condensates confined to highly anisotropic
% cigar-shaped traps, in the mean-field regime. This equation is integrable
% and was solved by the inverse scattering transform introduced in the
% celebrated paper by Zakharov and Shabat (1972) (see [2]).
%
% This script implements a numerical scheme for solving the Korteweg-de
% Vries equation. The script is a complete solver of the KdV equation. That
% is, it starts with setting up a uniform mesh on which the numerical
% solution will be computed, discretize the solution in space and solve the
% resulted system of ODEs in time. For more explanation on the
% numerical strategy, type "help kdv".
%
% To improve numerical stability, we adopt the strategy outlined in [1].
% Specifically, we eliminate the dispersive term u_{xx} in the Fourier
% space by multiplying the entire equation (in Fourier space) by the
% integrating factor exp(q(t)), where 
%       q(t) = 0.5*1i*k^2*t,
% and applying the transformation
%       fft(v) = exp(q(t)) fft(u).
% The new unknown would be fft(v) = \hat{v} and the equation must be
% slightly modified (see [1], Chapter 10, page 111). The solution is
% finally obtained by the inverse transformation
%       u = ifft( exp(-q(t)) * fft(v) ).
%
% References
% [1] <a href="https://en.wikipedia.org/wiki/Nonlinear_Schr%C3%B6dinger_equation"
% > Nonlinear Schroedinger equation on Wikipedia </a>
% [2] <a
% href="http://zakharov75.itp.ac.ru/static/local/zve75/zakharov/1972/1972-05-e_034_01_0062.pdf"
% > Exact Theory of Two-dimensional Self-focusing and One-dimensional
% Self-modulation of Waves in Nonlinear Media </a>

%% Problem domain and initial condition
[x, k] = generateMesh1D(-pi, pi, 2^6);
tmax = 1;           % maximum simulation time

% Initial condition: We use one-soliton solution to generate the initial
% condition. One can use any arbitrary (sufficiently regular) functions for
% the initial condition.
A = 2;  kk = 1;  omega = A^2 + kk^2/2; 
t0 = 0;
theta0 = kk*x - omega*t0;
iu = A * exp(1i * theta0);
theta = kk*x - omega*tmax;
u_exact = A * exp(1i * theta);

%% Time integration in the Fourier space
tstart = tic;
%--------------------------------------
q = 0.5 * 1i*k.^2;
myfunc = @(t, vFourier) odefunc(t, vFourier, q);
tspan = linspace(0, tmax, 101);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);
[t, vFourier] = ode45(myfunc, tspan, fft(iu), options);
u = ifft(exp(-t * (q.')) .* vFourier, [], 2);
%--------------------------------------
telapsed = toc(tstart); fprintf('Time processed: %f \n', telapsed);

%% Plot the solution
% close all
figure
[xx, tt] = meshgrid(x, t);
subplot(2,2,1)
mesh(xx, tt, real(u)); view(3)
subplot(2,2,2)
mesh(xx, tt, imag(u)); view(3);
subplot(2,2,3)
u_diff = u(end,:) - transpose(u_exact);
plot(x, real(u_diff), x, imag(u_diff), 'LineWidth', 1);

%% Helping function 
function y = odefunc(t, vFourier, q)
    uFourier = exp(-q*t).*vFourier;
    y = -1i*exp(q*t) .* fft( ifft(uFourier) .* abs(ifft(uFourier)).^2  );
end
