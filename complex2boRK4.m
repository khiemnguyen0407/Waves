%complex2boRK4  Spectral method for the 2BO equation in complex form.

%% Preamble.
% Coefficients for the equaiton.
g = 1/2;
HilbertSign = +1;
% Space-time mesh for the 1D wave equation.
% The equation is solved on the interval [-L,L].
L = 800*pi;      disp(['L = ', num2str(L/pi), '*pi']);
N = 2^19;       disp(['N = 2^', num2str(length(factor(N)))]);
h = L/(N/2);    disp(['h = ', num2str(h)]);
% Scale factor transforming the Fourier Transform on [-pi,pi] --> [-L,L]
gamma = L/pi;
x = gamma*(2*pi/N)*(-N/2:N/2-1);
% disp(['Suggested time step dt = ' num2str(h^3*0.1)]);
dt = 0.4e-4;      disp(['dt = ', num2str(dt)]);
T = 50;         disp(['T = ', num2str(T)]); disp(['T/dt = ' num2str(T/dt)]);
disp(['dt/h^2 = ' num2str(dt/h^2)]);
nstep = round(T/dt);
nslide = 100;
tslide = 1:round(nstep/nslide):nstep+1;
nplot = length(tslide); t = zeros(nplot,1);
% Preallocate the vectors of solution.
u = zeros(nplot,N);

% INITIAL CONDITION.
% % I. For periodic solution.
% A = 3; k = sqrt(2)/2; omega = g*(k^2 + A^4*pi^2)/2;
% t0 = 0; theta0 = 1i*(k*x - omega*t0);
% uNow = A*exp(theta0);
% theta = 1i*(k*x - omega*T);
% uExact = A*exp(theta);
% % Consider initial condition as the solution at the current time step.
% u(1,:) = uNow;

% II. For dispersive shock waves with nearly discontinuous jump.
rL = [+1.5    +0.5];  
vL = [-3.0    +3.0];  
rm = rL(1);  rp = rL(2);  
vm = vL(1);  vp = vL(2);
disp(vertcat(rL,vL));
% Compute intermediate step according to the theory.
ri = 1/2*((rm + rp) + (vm - vp)/(pi*g));
vi = pi*g/2*((vm - vp) + (vm + vp)/(pi*g));
% Construct initial condiution using tanh regularization.
zeta = 1; disp(['zeta = ', num2str(zeta)]);
stretch = 40.0; disp(['stretch = ', num2str(stretch)]);

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
ir2 = (rm - rp)/2*( tanh(zeta*(x-L1)) - tanh(zeta*(x-L2)) ) + rp;

% Initial condition for \phi which is defined by \phi_x = g v.
iphi = zeros(size(x));
iphi(idxL) = ( vp*x(idxL) + L1*(vm - vp) )/g;
iphi(idxM) = ( vm*x(idxM) )/g;
iphi(idxR) = ( vp*x(idxR) + L2*(vm - vp) )/g;
AL = 0.5*(vm + vp);
BL = (0.25*(vm - vp)/zeta)*((2*L1+delta)*zeta - 2*log(cosh(0.5*delta*zeta)));
AR = 0.5*(vm + vp);
BR = (0.25*(vm - vp)/zeta)*((2*L2-delta)*zeta + 2*log(cosh(0.5*delta*zeta)));
iphi(midIdxL) = ( AL*x(midIdxL) ...
    + (0.5*(vm-vp)/zeta)*log(cosh(zeta*(x(midIdxL) - L1))) + BL ) /g;
iphi(midIdxR) = ( AR*x(midIdxR) ...
    + (0.5*(vp-vm)/zeta)*log(cosh(zeta*(x(midIdxR) - L2))) + BR ) /g;
iv = horzcat([vp, g*diff(iphi)/h]);
% Consider initial condition as the solution at the current time step.
u(1,:) = sqrt(ir).*exp(1i*iphi);

% Wave number for odd-order derivative.
ko = [0:N/2-1  0  -N/2+1:-1]/gamma;
% Wave number for even-order derivative.
ke = [0:N/2   -N/2+1:-1]/gamma;

%% Iteration for soluion over time.
% Index for indicating layer of storage.
j = 2;
piSH = pi*pi/2;
absk = abs(ko);
kappa = -1i*g*dt;
ik2 = 1i*ke.^2;
E = exp(-0.25*dt*ik2*g);
E2 = E.^2;
uhat = fft(u(1,:));
for n = 1:nstep
    b1 = ifft( uhat );
    a1 = kappa.*( piSH*fft( abs(b1).^4.*b1 ) ...
        - HilbertSign*pi*fft( b1.*ifft(absk.*fft(abs(b1).^2)) ) );
    b2 = ifft( E.*(uhat + 0.5*a1) );
    a2 = kappa.*( piSH*fft( abs(b2).^4.*b2 ) ...
        - HilbertSign*pi*fft( b2.*ifft(absk.*fft(abs(b2).^2)) ) );
    b3 = ifft( E.*uhat + 0.5*a2 );
    a3 = kappa.*( piSH*fft( abs(b3).^4.*b3 ) ...
        - HilbertSign*pi*fft( b3.*ifft(absk.*fft(abs(b3).^2)) ) );
    b4 = ifft( E2.*uhat + E.*a3 );
    a4 = kappa.*( piSH*fft(abs(b4).^4.*b4 ) ...
        - HilbertSign*pi*fft( b4.*ifft(absk.*fft(abs(b4).^2)) ) );
    % Update solution.
    uhat = E2.*uhat + (E2.*a1 + 2*E.*(a2 + a3) + a4)/6;
    
    % Double-check if the numerical solution is unstable to stop the loop.
    if ~isempty(find(isnan(uhat),1))
        disp(['Solution is stable up to time T = ', num2str(n*dt)]);
        return
    end
    
    % If the time step coincides the SLIDE that needs to be stored, we keep
    % it in the solution vector.
    if ~isempty(find(n+1 == tslide,1))
        u(j,:) = ifft(uhat);
        t(j) = n*dt;        
        j = j+1;
        % Print out process of simulation.
        disp(['Running process:  ' num2str(n/nstep*100,3), '%']);
        % Print out intermediate numerical results.
    end
end
disp(['The solution is stable up to time T = ', num2str(t(end)), '.']);