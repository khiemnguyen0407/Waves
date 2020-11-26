%bo  Numerical scheme for solving the Benjamin-Ono equation
% u_t + alpha*uu_x + delta*Hu_xx = 0 

% This sign dictates the way we define the Hilbert transform.
HilbertSign = 1;

% Space-time mesh of the wave problem
% The equation is solved on the interval [-L,L].
L = 800*pi;  	disp(['L = ', num2str(L/pi), '*pi']);
N = 2^12;       disp(['N = 2^', num2str(length(factor(N)))]);
h = L/(N/2);    disp(['h = ', num2str(h)]);
% Scale factor transforming the Fourier Transform on [-pi,pi] --> [-L,L] 
scale = L/pi; 
x = scale*(2*pi/N)*(-N/2:N/2-1);

dt = 1e-4;
disp(['dt = ', num2str(dt)]);

T = 100;  
disp(['T = ', num2str(T)]); disp(['T/dt = ' num2str(T/dt)]);

nStep = round(T/dt);
nSlide = 100; tSlide = 1:round(nStep/nSlide):nStep+1;
nplot = length(tSlide); t = zeros(nplot,1);

% Solution vector.
u = zeros(nplot,N);

% INITIAL CONDITION.
% % For exact 1-phase periodic solution (u_t + 2 u u_x + H u_{xx} = 0).
% c = -1; a = 1; b = 2;
% t0 = 0; theta0 = (a - b)*x - (a^2 - b^2)*t0;
% u(1,:) = c + a - b + 2*(b-a)*(b-c - sqrt((b-c)*(a-c))*cos(theta0))...
%     ./ (a+b-2*c - 2*sqrt((b-c)*(a-c))*cos(theta0));
% [xx,tt] = meshgrid(x,0:(T/nSlide):T);
% theta = (a-b)*xx - (a^2 - b^2)*tt;
% uExact = c + a - b + 2*(b-a)*(b-c - sqrt((b-c)*(a-c))*cos(theta))...
%     ./ (a+b-2*c - 2*sqrt((b-c)*(a-c))*cos(theta));

% For dispersive shock waves (Using strict lines).
uL = [2, 1];
disp(uL);
Vs = 2*uL(1);
Cg = 2*uL(2);
As = 4*(uL(1) - uL(2));
% % Generate the initial-condition at time t = 0.
jumpLocation = [-L, 0];

% Initial condition is stored at the first vector of solution.
% Construct initial condition using tanh regularization.
delta = 5*ones(1,2);
iu = (uL(1)-uL(2))/2*( tanh(delta(1)*(x-jumpLocation(1))) - ...
    tanh(delta(2)*(x-jumpLocation(2))) ) + uL(2);
u(1,:) = iu;

% Numerical integration for wave solution using pseudo-spectral method.
% Wave number for odd-order derivative.
ko = [0:N/2-1  0  -N/2+1:-1]/scale;
% Wave number of even-order derivative.
ke = [0:N/2  -N/2+1:-1]/scale;
% Index for indicating layer of storage.
j = 2;
%LOOP FOR WAVE SOLUTION.
g = -1i*dt*ko; % u_{t} + 2*u u_x + Hu_{xx} = 0.
ik2 = 1i*ke.*ke.*sign(ko);
E = exp(HilbertSign * 0.5*dt*ik2);
E2 = E.^2;
uhat = fft(u(1,:));
tstart = tic;
for n = 1:nStep
    b1 = real( ifft(uhat) );
    a1 = g.*fft( b1.^2 );
    b2 = real( ifft(E.*(uhat+ 0.5*a1)) );
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
        % nAs = max(u(j,x > L2)) - uL(2);
        t(j) = n*dt;
        j = j+1;
        % Print out intermediate results.
        disp(['Running process:  '  num2str(round(n/nStep*100)) '%']);
    end
end
% If the solution is stable up the last time instant, display "success".
if j > size(u,1)
    disp(['The solution is stable up to this step T = ', num2str(t(end)) '.']);
end
telapsed = toc(tstart);

%% Plot the solution
% close all
% % Plot direct numerical solution.
% tIndex = round(nSlide/4) * (1:4)+1;
% figure('Position', [100, 75, 1200, 800])
% for i = 1:4
%     subplot(2,2,i)
%     plot(x, u(tIndex(i), :)), hold on
%     grid on
% end

% % Plot the modulation solution on top of direct numerical solution.
% figure
% tIndex = length(t);
% tt = t(tIndex);
% xx = linspace(0, tt, 1e4+1);
% aa = (0.5/tt)*xx;
% bb = 0.5*uL(1);
% cc = 0;
% theta = (bb - aa).*xx - (bb.*bb - aa.*aa)*tt;
% uu = 2*(bb-aa).^2 ./ (aa + bb - 2*cc - 2*sqrt((aa-cc).*(bb-cc)).*cos(theta)) + 2*cc;
% plot(x, u(tIndex,:), xx, uu, 'LineWidth', 1.0), grid on
% legend('  Direct numerical solution', '  Modulation solution', ...
%     'FontSize', 16, 'Location', 'best')
