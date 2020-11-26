% DNLS  Spectral method for the defocusing Nonlinear Schroedinger equation.

% Space-time mesh for the 1D wave equation.
% The equation is solved on the interval [-L,L].
L = 300*pi;     disp(['L = ', num2str(L/pi), '*pi']);
N = 2^17;       disp(['N = 2^', num2str(length(factor(N)))]);
h = L/(N/2);    disp(['h = ', num2str(h)]);
% Scale factor transforming the Fourier Transform on [-pi,pi] --> [-L,L]
gamma = L/pi;
x = gamma*(2*pi/N)*(-N/2:N/2-1);
% disp(['Suggested time step dt = ' num2str(h^3*0.1)]);
dt = 1.0e-4;      disp(['dt = ', num2str(dt)]);
T = 40;         disp(['T = ', num2str(T)]); disp(['T/dt = ' num2str(T/dt)]);
disp(['dt/h^2 = ' num2str(dt/h^2)]);
nstep = round(T/dt);
nslide = 500; tSlide = 1:round(nstep/nslide):nstep+1;
nplot = length(tSlide); t = zeros(nplot,1);
% Preallocate the vectors of solution.
u = zeros(nplot,N);

% INITIAL CONDITION.
%====================================
% % I. For exact periodic solution.
% A = 2; k = sqrt(2); omega = A^2 + k^2/2;
% t0 = 0; theta0 = k*x - omega*t0;
% u(1,:) = A*exp(1i*theta0);
% theta = k*x - omega*T;
% uExact = A*exp(1i*theta);

% II. For Dispersive Shock Waves (Using tanh-function).
rL = [1   0.5];
vL = [-2   2];
rm = rL(1); rp = rL(2);
vm = vL(1); vp = vL(2);
disp(num2str(vertcat(rL,vL)));
% Construct initial condiution using tanh regularization.
zeta = 2;   disp(['zeta = ' num2str(zeta)]);
stretch = 17.5;  disp(['stretch = ', num2str(stretch)]);
delta = round(stretch/h)*h;

L1 = -L/2;
L2 = +L/2;
idxL = find(x <= L1-delta/2);
idxM = find(x >= L1+delta/2 & x <= L2-delta/2);
idxR = find(x >= L2+delta/2);
midIdxL = find(L1-delta/2 < x & x < L1+delta/2);
midIdxR = find(L2-delta/2 < x & x < L2+delta/2);
% Initial condition for \rho.
ir = zeros(size(x));
ir(idxL) = rp*ones(size(idxL));
ir(idxM) = rm*ones(size(idxM));
ir(idxR) = rp*ones(size(idxR));
ir(midIdxL) = 0.5*(rm + rp) + 0.5*(rm - rp)*tanh(zeta*(x(midIdxL)-L1));
ir(midIdxR) = 0.5*(rm + rp) + 0.5*(rp - rm)*tanh(zeta*(x(midIdxR)-L2));

% Initial condition for \phi which is defined by \phi_x = g v.
iphi = zeros(size(x));
iphi(idxL) = ( vp*x(idxL) + L1*(vm - vp) );
iphi(idxM) = ( vm*x(idxM) );
iphi(idxR) = ( vp*x(idxR) + L2*(vm - vp) );
AL = 0.5*(vm + vp);
BL = (0.25*(vm - vp)/zeta)*((2*L1+delta)*zeta - 2*log(cosh(0.5*delta*zeta)));
AR = 0.5*(vm + vp);
BR = (0.25*(vm - vp)/zeta)*((2*L2-delta)*zeta + 2*log(cosh(0.5*delta*zeta)));
iphi(midIdxL) = ( AL*x(midIdxL) ...
    + (0.5*(vm-vp)/zeta)*log(cosh(zeta*(x(midIdxL) - L1))) + BL );
iphi(midIdxR) = ( AR*x(midIdxR) ...
    + (0.5*(vp-vm)/zeta)*log(cosh(zeta*(x(midIdxR) - L2))) + BR );
iv = horzcat([vp, diff(iphi)/h]);
% Consider initial condition as the solution at the current time step.
u(1,:) = sqrt(ir).*exp(1i*iphi);

% Wave number for odd-order derivative.
ko = [0:N/2-1  0  -N/2+1:-1]/gamma;
% Wave number for even-order derivative.
ke = [0:N/2   -N/2+1:-1]/gamma;

% LOOP FOR WAVE SOLUTION.
% Index for indicating layer of storage.
j = 2;
% Intermediate constants for the loop.
g = -1i*dt;
ik2 = 1i*ke.^2;
E = exp(-0.25*dt*ik2);
E2 = E.^2;
uhat = fft(u(1,:));
for n = 1:nstep
    b1 = ifft( uhat );
    a1 = g.*fft( abs(b1).^2 .* b1 );
    b2 = ifft( E.*(uhat + 0.5*a1) );
    a2 = g.*fft( abs(b2).^2 .* b2 );
    b3 = ifft( E.*uhat + 0.5*a2 );
    a3 = g.*fft( abs(b3).^2 .* b3 );
    b4 = ifft( E2.*uhat + E.*a3 );
    a4 = g.*fft( abs(b4).^2 .* b4 );
    uhat = E2.*uhat + (E2.*a1 + 2*E.*(a2 + a3) + a4)/6;
    
    % Double-check if the numerical solution is unstable to stop the loop.
    if ~isempty(find(isnan(uhat),1))
        disp(['Solution is unstable up to time T = ', num2str(n*dt)]);
        return
    end
    
    % If the time step coincides the SLIDE that needs to be stored, we keep
    % it in the solution vector.
    if ~isempty(find(n+1 == tSlide,1))
        u(j,:) = ifft(uhat);
        t(j) = n*dt;
        j = j+1;
        % Print out intermediate results.
        disp(['Running process:  ' num2str(round(n/nstep*100)), '%']);
        fprintf('\n');
    end
end
if j > size(u,1)
    disp(['The solution is stable up to time T = ', num2str(t(end)), '.']);
end