%COMPLEX2BO   Solve the 2-Benjamin-Ono equation in the complex form
%
% COMPLEX2BO implements a numerical scheme for solving the bidirectional
% Benjamin-Ono (2BO) equation given in the complex form. The equation
% represents the complex form of the Calogero–Sutherland (CS) model (see
% [1] and see [2]). The 2BO equation reads
% 
% \frac{1i}{g} u_t = -\frac{1}{2} u_{xx} + \frac{\pi^2}{2} u | u |^4 
%                    + \pi \HilbertTransform[\partial_{x} |u|^2 ] u].
%
% This equation can be written in the hydrodynamic form using the density
% variable \rho and the velocity v. Such formulation is rather complicated
% and can be given by system (6)-(7) in the paper [1]. In the complex form,
% the 2BO equation is similar to the nonlinear Schroedinger (NLS) equation.
% Thus, the solution is also complex-valued like the solution of the NLS
% equation.
%
% The script is a complete solver of the 2BO equation. That is, it starts
% with setting up a uniform mesh on which the numerical solution will be
% computed, discretize the solution in space and solve the resulted system
% of ODEs in time. Unlike the solver for KdV equation (see kdv.m), this
% equation cannot be solved by other methods such as the finite difference
% and finite element methods due to the appearance of the Hilbert
% transform.  The Hilbert transform of function f = f(x) is defined as the
% convolution of a singular kernel with the function itself. Thus, although
% it is possible to approximate a Hilbert transform in a numerical basis,
% it is not trivial to integrate such approximations into a numerical
% scheme for the wave solution. For more explanation on the numerical
% strategy, type "help kdv".
%
% The Calogero–Sutherland (CS) model and the associated 2BO equation are
% studied extensively in [1]. The 2BO equation is integrable and has
% multi-soliton solutions that can be found also in [1]. The dispersive
% shock waves (DSWs) governed by the hydrodynamic form of the CS model are
% studied in [2].
%
% References
% [1] <a
% href="https://iopscience.iop.org/article/10.1088/1751-8113/42/13/135201"
%> Integrable hydrodynamics of CS model: 2BO equation</a>
% [2] <a href="https://arxiv.org/abs/1705.09996"
% > DSWs in systems with nonlocal dispersion of Benjamin-Ono type </a>

%% Problem domain and initial condition
[x, k] = generateMesh1D(-pi/(sqrt(2)/2), pi/(sqrt(2)/2), 2^10);
HilbertSign = 1;  g = 1/2;        % coefficients of the equation
tmax = 0.5;
% Initial condition: We use the exact one-phase periodic solution.
A = 3; 
kk = sqrt(2)/2; 
omega = g*(kk^2 + A^4*pi^2)/2;
t0 = 0; theta0 = 1i*(kk*x - omega*t0);
iu = A * exp(theta0);

%% Time integration in the Fourier space
tstart = tic;

ig = 1i*g; 
half_kSquare = 0.5*k.^2;
absk = abs(k);
half_PiSquare = 0.5*pi*pi;
myfunc = @(t, uFourier) ...
    odefunc(t, uFourier, HilbertSign, ig, half_kSquare, absk, half_PiSquare);
tspan = linspace(0, tmax, 51);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

[t, uFourier] = ode45(myfunc, tspan, fft(iu), options);
u = ifft(uFourier, [], 2);

telapsed = toc(tstart); fprintf('Time processed: %f \n', telapsed);
%% Plot the solution
close all
[xx, tt] = meshgrid(x, t);
subplot(2,2,1)
mesh(xx, tt, real(u)); view(3)
subplot(2,2,2)
mesh(xx, tt, imag(u)); view(3);

theta = 1i*(kk*xx - omega*tmax);
u_exact = A * exp(theta);
u_diff = u(end,:) - u_exact;
subplot(2,2,3)
plot(x, real(u_diff), x, imag(u_diff));

%% Helper function
function y = odefunc(~, uFourier, HilbertSign, ig, half_kSquare, absk, half_PiSquare)
    u = ifft(uFourier);
    a = half_kSquare .* uFourier;
    b = half_PiSquare .* fft(u .* abs(u).^4);
    c = -HilbertSign * pi * fft(u.*ifft( absk .* fft(abs(u).^2) ));
    y = -ig*(a + b + c);
end
