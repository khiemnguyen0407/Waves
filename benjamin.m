%BENJAMIN  Solve the Benjamin-Ono equation
%   u_t + u_x + 2*u*u_x - alpha*Hu_{xx} - beta * u_{xxx} = 0

% Space-time mesh of the wave problem
% The equation is solved on the interval [-L,L].
L = 2000*pi;  	
N = 2^14+1;
disp(['L = ', num2str(L/pi), '*pi']);
disp(['N = 2^', num2str(length(factor(N-1))), '+1']);
% Scale factor transforming the Fourier Transform on [-pi,pi] --> [-L,L] 
scale = L/pi;
x = scale*(2*pi/N)* (-(N-1)/2 : (N-1)/2);
k = [0: (N-1)/2,  -(N-1)/2: -1]/scale;
T = 500;

%% Parameters for Benjamin equation.
alpha = -1; % alpha = +1; 
beta = -1;  % beta = +1;

%% Initial condition.
% For dispersive shock waves.
uLimit = [1, 0];
delta = 2.5*ones(1,2);

jumpLocation = [-L/2, L/2];
iu = 0.5*(uLimit(1) - uLimit(2))*(tanh(delta(1)*(x-jumpLocation(1))) ...
        - tanh(delta(2)*(x-jumpLocation(2)))) + uLimit(2);
iu = iu(:); % Colum vector.


% Numerical integration for wave solution using pseudo-spectral method.
ik = 1i*k;
q = -(1i*k + alpha*1i*k.*abs(k) + 1i*beta*k.^3);
% Put the algorithmic variables into the column vector.
ik = ik(:);
q = q(:);

% Setup the right-hand side function f(t,y) in the equation y'(t) = f(t,y)
% in the Fourier space.
ODEFUNC = @(t, u_fourier) q.*u_fourier - ik.*fft(ifft(u_fourier).^2);
tspan = linspace(0, T, 101);
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-4);
% Solve the wave equation in Fourier space using ODE45 by MATLAB.
[t, sol_fourier] = ode45(ODEFUNC, tspan, fft(iu), options);
% Retrieve the solution in real space.
u = ifft(sol_fourier, [], 2);
%% Plot the solution.
close all
tIndex = 26:25:101;
figure('Position', [100, 100, 1600, 800])
for i = 1:length(tIndex)
    subplot(2,2,i)
    tt = t(tIndex);
    xx = linspace(min(x), max(x), 16*(N-1)+1);
    uu = interp1(x, u(tIndex(i),:), xx, 'spline');
    plot(xx, uu, 'k-', x, u(1,:), '--', 'LineWidth', 1), grid on
    text(0, 2.25, ['$t = ', num2str(t(tIndex(i))), '$'], ...
        'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'center');
end