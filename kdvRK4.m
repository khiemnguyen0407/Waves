%kdvRK4 Solve the KdV equation with a spectral method
% 
% This is a script implementing using the 4th-order Runge-Kutta integration
%scheme in time and spectral method in space.
% Set up grid and initial data:
L = 250*pi;
N = 2^12;
h = L/(N/2);
% Scale factor transforming the Fourier Transform on [-pi,pi] --> [-L,L] 
scale = L/pi;
x = scale*(2*pi/N)*(-N/2:N/2-1);
% disp(['Suggested time step dt = ' num2str(h^3*0.1)]);
dt = 0.001; T = 100; nStep = round(T/dt);
nSlides = 100; tSlides = 1:round(nStep/nSlides):nStep+1;
nplot = length(tSlides); t = zeros(nplot,1);

u = zeros(nplot,N);     % Solution vector

%% INITIAL CONDITION.
% For regularized step-like condition (top hat initial condition).
uL = [2   1.0];
% Generate initial condition at time t = 0.
delta = 20*h;
hL = 2*L/3;
negIdxL = find(x <= -hL-delta/2);
negIdxR = find(x >= -hL+delta/2 & x <= -delta/2);
posIdx = find(x >= +delta/2);
midIdxL = find(-hL-delta/2 < x & x < -hL+delta/2);
midIdxR = find(-delta/2 < x & x < delta/2);
iu = zeros(size(x));
iu(negIdxL) = uL(2)*ones(size(negIdxL));
iu(negIdxR) = uL(1)*ones(size(negIdxR));
iu(posIdx) = uL(2)*ones(size(posIdx));
iu(midIdxL) = (uL(1)-uL(2))/delta*(x(midIdxL)+hL) + (uL(1)+uL(2))/2;
iu(midIdxR) = (uL(2)-uL(1))/delta*x(midIdxR) + (uL(1)+uL(2))/2;
% Initial condition is stored at the first vector of solution.
u(1,:) = iu;
% Wavenumber.
k = [0:N/2-1 0 -N/2+1:-1]/scale;

%% LOOP FOR SOLUTION.
j = 2;
g = -0.5i*dt*k;  % u_t + u u_x + u_{xxx} = 0.
ik3 = 1i*k.^3;
E = exp(dt*ik3/2);
E2 = E.^2;
uhat = fft(u(1,:));
hwait = waitbar(0, 'please wait ...');
for n = 1:nStep
    b1 = real( ifft(uhat) );
    a1 = g.*fft( b1.^2 );
    b2 = real( ifft(E.*(uhat+0.5*a1)) );
    a2 = g.*fft( b2.^2);
    b3 = real( ifft(E.*uhat + 0.5*a2) );
    a3 = g.*fft( b3.^2);
    b4 = real( ifft(E2.*uhat + E.*a3) ); 
    a4 = g.*fft( b4.^2);
    % Update solution in Fourier space at next time step.
    uhat = E2.*uhat + (E2.*a1 + 2*E.*(a2 + a3) + a4)/6;
    if ~isempty(find(isnan(uhat),1))
        fprintf('The solution is unstable up to time T = %f', n*dt);
        delete(hwait); 
        return
    end
    if ~isempty(find(n+1==tSlides, 1, 'first'))
        u(j,:) = real( ifft(uhat) );
        t(j) = n*dt;
        j = j+1;
        waitbar(n/nStep,hwait, [num2str(round(n/nStep*100)) '%']);
    end
end
% If the solution is stable up the last time instant, display "success".
if j > size(u,1)
    disp(['The solution is stable up to this step T = ', num2str(t(end)) '.']);
    delete(hwait);
end