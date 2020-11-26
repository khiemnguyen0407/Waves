%bboRK4  Spectral method for the Boussinesq Benjamin-Ono equation in the
%        potential form.

%% Algorithmic parameters.

% Coefficients for the equation.
p = 10;
HilbertSign = 1;

% Space-time mesh for the 1D wave equation.
% The equation is solved on the interval [-L, L].
L = 2*pi;     % problem domain is [-L, L]
N = 2^6;       % number of grid points.
h = 2*L/N;      % mesh spacing
% Scale factor transforming the Fourier transform on [-pi,pi] --> [-L, L]
scale = L/pi;   % To be precise, scale = 2*L / 2*pi
dt = 1e-6;      % time step
T = 2;         % total time for simulation. The last time is also T.
nstep = round(T/dt);    % number of time steps
nslide = 100;   % number of solution slices.
tslide = 1:round(nstep/nslide) : nstep +1;  % vector of the time slices
nplot = length(tslide);

x = scale*(2*pi/N)*(-N/2:N/2-1);

% Wavenumbers
ko = [0: N/2-1,  0,  -N/2+1: -1]/scale;     % for odd-order derivative.
ke = [0: N/2,         -N/2+1:-1]/scale;     % for even-order derivative.

% Solution vectors (to be precise, two matrices)
v = zeros(nplot, N);
u = zeros(nplot, N);

% Initial condition for exact periodic solution.
initialOption = 'exact';
if strcmpi(initialOption, 'exact')
    k = 0.5;
    A = 0.5;
    Q = 1;
    beta = 0;
    P = 2*k*sqrt(Q*Q - A*A);
    velocity = sqrt(k * Q /sqrt(Q*Q - A*A) + p + beta);
    omega = k * velocity;
    u(1,:) = P./(Q + A*cos(k*x)) + beta;
    v(1,:) = omega/k*u(1,:);
elseif strcmpi(initialOption, 'dsw')
    % Initial condition for generating Dispersive shock waves.
    bL = [3, 0];
    delta = [5, 5];
    jumpLocation = 0.5 * [-L, L];
    v(1,:) = zeros(1, N);
    u(1,:) = 0.5*(bL(1) - bL(2))*(tanh(delta(1)*(x-jumpLocation(1))) ...
        - tanh(delta(2)*(x-jumpLocation(2)))) + bL(2);
end

%% 4th-order Runge-Kutta method for solution over time.
% Index for indicating layer of storage.
j = 2;

% Time recorder.
t = zeros(1, nplot);
% Fourier transforms of a & b.
v_fourier = fft( v(1,:) );
u_fourier = fft( u(1,:) );
q = 1i*p*ko + HilbertSign*1i*ke.* ke .* sign(ko);
ikHalf = 0.5*1i*ko;
ik = 1i*ko;
dtHalf = 0.5*dt;
for n = 1:nstep
    g1a = dt*(q.*u_fourier + ikHalf.*fft( ifft(u_fourier).^2 ));
    g1b = dt*ik.*v_fourier;
    
    lambda2 = (u_fourier + 0.5*g1a);
    g2a = dt*(q.*lambda2 + ikHalf.*fft( ifft(lambda2).^2 ));
    g2b = dt*ik.*( v_fourier + 0.5*g1b );
    
    lambda3 = u_fourier + 0.5*g2a;
    g3a = dt*(q.*lambda3 + ikHalf.*fft( ifft(lambda3).^2 ));
    g3b = dt*ik.*( v_fourier + 0.5*g2b );
    
    lambda4 = u_fourier + g3a;
    g4a = dt*(q.*lambda4 + ikHalf.*fft( ifft(lambda4).^2 ));
    g4b = dt*ik.*( v_fourier + g3b );
    
    % Update the solution in the Fourier space.
    v_fourier = v_fourier + 1/6 * (g1a + 2*(g2a + g3a) + g4a);
    u_fourier = u_fourier + 1/6 * (g1b + 2*(g2b + g3b) + g4b);
    
    % Check whether the numerical solution is stable.
    if ~isempty( find(isnan(v_fourier), 1) )
        disp(['The numerical solution is unstable upto time T = ', num2str(n*dt)]);
        return
    end
    
    if ~isempty( find(n+1 == tslide, 1) )
        v(j, :) = real( ifft(v_fourier) );
        u(j, :) = real( ifft(u_fourier) );
        t(j) = n * dt;
        j = j + 1;
        disp(['Running process: ', num2str(round(n/nstep * 100)), '%']);
    end
end

%%
close all
if strcmpi(initialOption, 'exact')
    figure('Position', [100, 75, 1200, 500])
    subplot(1,2,1)
    plot(x, u(1,:),'k-', 'LineWidth', 1.0), grid on
    subplot(1,2,2)
    plot(x, u(end,:), x, P./(Q + A*cos(k*x - omega*T)), 'LineWidth', 1.0);
    grid on
end