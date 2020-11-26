%BENJAMINRK4  Solve the Benjamin-Ono equation
% u_t + u_x + 2*u*u_x - alpha*Hu_{xx} - beta * u_{xxx} = 0

% Space-time mesh of the wave problem
% The equation is solved on the interval [-L,L].
L = 800*pi;  	disp(['L = ', num2str(L/pi), '*pi']);
N = 2^13;       disp(['N = 2^', num2str(length(factor(N)))]);
h = L/(N/2);    disp(['h = ', num2str(h)]);
% Scale factor transforming the Fourier Transform on [-pi,pi] --> [-L,L] 
scale = L/pi;
x = scale*(2*pi/N)*(-N/2:N/2-1);

dt = 5e-4;
disp(['dt = ', num2str(dt)]);

T = 50;
disp(['T = ', num2str(T)]); disp(['T/dt = ' num2str(T/dt)]);
nStep = round(T/dt);

nSlide = 100; tSlide = 1:round(nStep/nSlide):nStep+1;
nplot = length(tSlide); t = zeros(nplot,1);

%% Parameters for Benjamin equation.
alpha = -1;
beta = -1;
% alpha = +1;
% beta = +1;

% Solution vector.
u = zeros(nplot,N);

%% Initial condition.
% For dispersive shock waves.
uLimit = [1, 0];
delta = 2.5*ones(1,2);

jumpLocation = [-L, 0];
iu = 0.5*(uLimit(1) - uLimit(2))*(tanh(delta(1)*(x-jumpLocation(1))) ...
        - tanh(delta(2)*(x-jumpLocation(2)))) + uLimit(2);
u(1,:) = iu;

% Numerical integration for wave solution using pseudo-spectral method.
% Wave number for odd-order derivative.
ko = [0:N/2-1  0  -N/2+1:-1]/scale;
% Wave number of even-order derivative.
ke = [0:N/2  -N/2+1:-1]/scale;
% Index for indicating layer of storage.
j = 2;
%LOOP FOR WAVE SOLUTION.
g = -1i * dt * ko;
h = 1i*ko + alpha*1i*ke.*ke.*sign(ko) + 1i*beta*ko.^3;
E = exp(-0.5*dt*h);
E2 = E.^2;
uhat = fft(u(1,:));
for n = 1:nStep
    b1 = real( ifft(uhat) );
    a1 = g.*fft( b1.^2 );
    b2 = real( ifft(E.*(uhat + 0.5*a1)) );
    a2 = g.*fft( b2.^2 );
    b3 = real( ifft(E.*uhat + 0.5*a2) );
    a3 = g.*fft( b3.^2 );
    b4 = real( ifft(E2.*uhat + E.*a3) );
    a4 = g.*fft( b4.^2 );
    % Update solution in Fourier space at next time step.
    uhat = E2.*uhat + (E2.*a1 + 2*E.*(a2 + a3) + a4)/6;
    if ~isempty(find(isnan(uhat),1))
        disp(['The solution is unstable up to time T = ', num2str(n*dt)]);
        return
    end
    if ~isempty(find(n+1==tSlide,1))
        u(j,:) = real( ifft(uhat) );
        nAs = max(u(j,x > L2)) - uL(2);
        t(j) = n*dt;
        j = j+1;
        disp(['Running process:  '  num2str(round(n/nStep*100)) '%']);
    end
end
% If the solution is stable up the last time instant, display "success".
if j > size(u,1)
    disp(['The solution is stable up to this step T = ', num2str(t(end)) '.']);
end