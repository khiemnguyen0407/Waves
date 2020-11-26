%bboRK4  Spectral method for the Boussinesq Benjamin-Ono equation in the
%        potential form.

%% Algorithmic parameters.
% Coefficients for the equation.
p = 10.0;     
disp(['p = ', num2str(p)]);
HilbertSign = -1;	disp(['Sign of Hilbert Transform = ', num2str(HilbertSign)]);

% Space-time mesh for the 1D wave equation.
% The equation is solved on the interval [-L, L].
% L = 1000*pi;
disp(['L = ', num2str(L/pi), '*pi']);
% N = 2^14;      % number of grid points.
h = 2*L/N;      % mesh spacing
% Scale factor transforming the Fourier transform on [-pi,pi] --> [-L, L]
scale = L/pi;   % To be precise, scale = 2*L / 2*pi
% T = 1500;         % total time for simulation. The last time is also T.
% x = scale*(2*pi/N)*(-N/2:N/2-1);
x = scale*(2*pi/N)* (-(N-1)/2 : (N-1)/2);

% Wavenumbers
ko = [0: N/2-1,  0,  -N/2+1: -1]/scale;     % for odd-order derivative.
ke = [0: N/2,         -N/2+1:-1]/scale;     % for even-order derivative.
% k = [0: (N-1)/2,  -(N-1)/2: -1]/scale;
% Initial condition for exact periodic solution.
initialCondition = 'dsw';
if strcmpi(initialCondition, 'periodic')
    kk = 1;  A = 1;  Q = 3;  beta = 0;
    P = 2*kk*sqrt(Q*Q - A*A);
    velocity = sqrt( (kk*Q)/sqrt(Q*Q - A*A) + p + beta);
    omega = kk * velocity;
    iu = P./(Q + A*cos(kk*x)) + beta;
    iv = -velocity*iu;
    u_exact = P./(Q + A*cos(kk*x - omega*T)) + beta;
    v_exact = -velocity*u_exact;
elseif strcmpi(initialCondition, 'soliton')
    beta = 0;
    velocity = 1.25;
    P = 4*(velocity^2 - p - beta);
    Q = 1;
    A = (velocity^2 - p - beta)^2;
    iu = P./(Q + A*x.^2) + beta;
    iv = -velocity*iu;
    u_exact = P./(Q + A.*(x - velocity*T).^2) + beta;
    v_exact = -velocity * u_exact;
elseif strcmpi(initialCondition, 'dsw')
    % Initial condition for generating Dispersive shock waves.
    vLimit = [0,  0];
    % uLimit = [2.00,  -0.25];
    
    % For the one fast DSW.
%     vLimit = [0, 10];
%     uLimit = [-10 + 5^(2/3)*(49 + 12*sqrt(10))^(1/3), 0];
    
    disp(['vLimit = ', num2str(vLimit)]);
    disp(['uLimit = ', num2str(uLimit)]);
    delta = 1.5*ones(1, 2);
    jumpLocation = [-L, 0];
    iv = 0.5*(vLimit(1) - vLimit(2))*(tanh(delta(1)*(x-jumpLocation(1))) ...
        - tanh(delta(2)*(x-jumpLocation(2)))) + vLimit(2);
    iu = 0.5*(uLimit(1) - uLimit(2))*(tanh(delta(1)*(x-jumpLocation(1))) ...
        - tanh(delta(2)*(x-jumpLocation(2)))) + uLimit(2);
    
    % Intermediate level
    uIntermediate = (0.75*(vLimit(2) - vLimit(1)) + 0.5*(p + uLimit(1))^(1.5)...
        + 0.5*(p+uLimit(2))^(1.5))^(2/3) - p;
    vIntermediate = 0.5*(vLimit(1)+vLimit(2)) + 1/3 * ((p + uLimit(2))^(3/2) - (p + uLimit(1))^(3/2));
    disp(['uIntermediate = ', num2str(uIntermediate)]);
    disp(['vIntermedidte = ', num2str(vIntermediate)]);
end
iv = iv(:);     % put initial solution into column vector
iu = iu(:);     % column vector

%% 4th-order Runge-Kutta method for solution over time.
q = 1i*p*ko + HilbertSign*1i*ke.* ke .* sign(ko);
ikHalf = 0.5*1i*ko;
ik = 1i*ko;

% q = 1i*p*k + HilbertSign*1i *k .* abs(k);
% ikHalf = 0.5*1i*k;
% ik = 1i*k;

% Put the algorithmic variables into the column vector.
q = q(:);
ikHalf = ikHalf(:);
ik = ik(:);

% Setup the right-hand side function f(t,y) in the equation y'(t) = f(t,y)
% in the Fourier space.
v_index = 1:N;
u_index = (1:N) + v_index(end);
ODEFUNC = @(t, u_fourier) ...
    vertcat(q.*u_fourier(u_index) + ikHalf.*fft( ifft(u_fourier(u_index)).^2 ), ...
    ik.*u_fourier(v_index));

% tspan = [linspace(0, T-1.1, 95), linspace(T-1, T, 11)];u
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-4);
[t, solution_fourier_matrix] = ode45(ODEFUNC, tspan, vertcat(fft(iv), fft(iu))); % default options
% [t, solution_fourier_matrix] = ode45(ODEFUNC, tspan, vertcat(fft(iv), fft(iu)), options);
% [t, solution_fourier_matrix] = ode15s(ODEFUNC, tspan, vertcat(fft(iv), fft(iu)), options);
% [t, solution_fourier_matrix] = ode23s(ODEFUNC, tspan, vertcat(fft(iv), fft(iu)), options);
% [t, solution_fourier_matrix] = ode23t(ODEFUNC, tspan, vertcat(fft(iv), fft(iu)), options);
% [t, solution_fourier_matrix] = ode23tb(ODEFUNC, tspan, vertcat(fft(iv), fft(iu)), options);

v = ifft(solution_fourier_matrix(:,v_index), [], 2);
u = ifft(solution_fourier_matrix(:,u_index), [], 2);

%% Plot the solution
% close all
fs = 20;    % font size
lw = 0.5;   % line width

rFunc = @(u,v) v - 2/3 * (p+ u)^(1.5);
lFunc = @(u,v) v + 2/3 * (p + u)^(1.5);
Vr = @(u, v) sqrt(p + u);
Vl = @(u, v) -sqrt(p + u);

%-----------------------------------------
figure('Position', [100, 100, 1400, 900]);
subplot(2,1,1)
plot(x, u(end,:), '-', x, u(1,:), '--', 'LineWidth', lw), grid on, hold on
if exist('u_exact', 'var'), plot(x, u_exact, 'ko-'); end
myTitle(['$u(x, t = ', num2str(t(end)), '),\quad p = ', num2str(p), '$'], fs)
generate2DLabels('$x$', '$u$', fs);
%-----------------------------------------
subplot(2,1,2)
plot(x, v(end,:), '-', x, v(1,:), '--', 'LineWidth', lw), grid on, hold on
if exist('v_exact', 'var'), plot(x, v_exact, 'ko-'); end

myTitle(['$v(x, t = ', num2str(t(end)), '),\quad p = ', num2str(p), '$'], fs) 
generate2DLabels('$x$', '$v$', fs);

%% 
function generate2DLabels(xl, yl, fs)
    xlabel(xl, 'Interpreter', 'latex', 'FontSize', fs)
    ylabel(yl, 'Interpreter', 'latex', 'FontSize', fs);
end

function myTitle(mytitle, fs)
    title(mytitle, 'Interpreter', 'latex', 'FontSize', fs); 
end
