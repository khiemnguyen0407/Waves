%NLS  Solve nonlinear Schroedinger equation

%% Space-time mesh
a = -pi;
b = pi;
N = 2^6+1;
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
tmax = 1;

% Initial condition: We use one-soliton solution to generate the initial
% condition. One can use any arbitrary (sufficiently regular) functions for
% the initial condition.
A = 2;  kk = 1;  omega = A^2 + kk^2/2; 
t0 = 0;
theta0 = kk*x - omega*t0;
iu = A * exp(1i * theta0);
theta = kk*x - omega*tmax;
u_exact = A * exp(1i * theta);

%% Time integration using pseudo-spectral method
tstart = tic;
q = 0.5 * 1i*k.^2;
myfunc = @(t, vFourier) odefunc(t, vFourier, q);
tspan = linspace(0, tmax, 101);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);
[t, vFourier] = ode45(myfunc, tspan, fft(iu), options);
u = ifft(exp(-t * (q.')) .* vFourier, [], 2);
telapsed = toc(tstart);
fprintf('Time processed: %f \n', telapsed);

%% Plot the solution
close all
[xx, tt] = meshgrid(x, t);
subplot(2,2,1)
mesh(xx, tt, real(u)); view(3)
subplot(2,2,2)
mesh(xx, tt, imag(u)); view(3);
subplot(2,2,3)
u_diff = u(end,:) - transpose(u_exact);
plot(x, real(u_diff), x, imag(u_diff), 'LineWidth', 2);

function y = odefunc(t, vFourier, q)
    uFourier = exp(-q*t).*vFourier;
    y = -1i*exp(q*t) .* fft( ifft(uFourier) .* abs(ifft(uFourier)).^2  );
end