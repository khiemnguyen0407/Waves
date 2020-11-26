%BO  Solve the Benjamin-Ono equation
%   u_t + alpha*u u_x + Hu_{xx} = 0
%
% BO implements a numerical scheme for solving the Benjamin-Ono (BO)
% equation. The script is a complete solver of the BO equation in that it
% starts with setting up a uniform mesh on which the numerical solution
% will be computed. Unlike the solver for KdV equation (see KDV), this
% equation cannot be solved by other methods such as the finite difference
% and finite element methods due to the appearance of the Hilbert
% transform. The Hilbert transform of function f = f(x) is defined as the
% convolution of a singular kernel with the function itself. Thus, although
% it is possible to approximate a Hilbert transform in a numerical basis,
% it is not trivial to integrate such approximations into a numerical
% scheme for the wave solution. For more explanation on the numerical
% strategy, type "help kdv".
%
% To improve numerical stability, we adopt the strategy outlined in [1].
% Specifically, we eliminate the dispersive term Hu_{xx} in the Fourier
% space by multiplying the entire equation (in Fourier space) by the
% integrating factor exp(-i*k*abs(k)*t) and applying the transformation
%       fft(v) = exp(-i*k*abs(k)*t) fft(u)
% The new unknown would be fft(v) = \hat{v} and the equation must be
% slightly modified (see [1], Chapter 10, page 111). The solution is
% finally obtained by the inverse transformation
%       u = ifft( exp(i*k*abs(k)*t) * fft(v) ).
%
% References
% [1] <a href="http://www.audentia-gestion.fr/Matlab/10.1.1.473.7647.pdf"
% > Spectral Methods in MATLAB </a>
% [2] <a href="https://link.springer.com/article/10.1007/BF02510262"
% > A numerical method for the Benjamin-On equation </a>
%
% See also KDV, NLS, DNLS, COMPLEX2BO, BENJAMIN

%%
% This sign dictates the way we define the Hilbert transform.
HilbertSign = +1;

% Space-time mesh of the wave problem
% The equation is solved on the interval [-L,L].
L = 1*pi;
N = 2^9;
% Scale factor transforming the Fourier Transform on [-pi,pi] --> [-L,L]
scale = L/pi;
x = scale*(2*pi/N)*(-N/2:N/2-1);
nSlide = 101;
% x = scale * (2*pi/N) * (-(N-1)/2 : (N-1)/2);
T = 5;

% Initial condition.
initialCondition = 'exact';
switch initialCondition
    case 'exact'
        % For exact 1-phase periodic solution (u_t + 2 u u_x + H u_{xx} = 0).
        c = -1; a = 1; b = 2;
        t0 = 0; theta0 = (a - b)*x - (a^2 - b^2)*t0;
        iu = c + a - b + 2*(b-a)*(b-c - sqrt((b-c)*(a-c))*cos(theta0))...
            ./ (a+b-2*c - 2*sqrt((b-c)*(a-c))*cos(theta0));
        [xx,tt] = meshgrid(x,0:(T/nSlide):T);
        theta = (a-b)*xx - (a^2 - b^2)*tt;
        uExact = c + a - b + 2*(b-a)*(b-c - sqrt((b-c)*(a-c))*cos(theta))...
            ./ (a+b-2*c - 2*sqrt((b-c)*(a-c))*cos(theta));
    case 'dsw'
        % For dispersive shock waves (Using strict lines).
        uL = [2, 1];
        disp(uL);
        Vs = uL(1);
        Cg = uL(2);
        As = 4*(uL(1) - uL(2));
        % % Generate the initial-condition at time t = 0.
        jumpLocation = [-L, 0];
        % Initial condition is stored at the first vector of solution.
        % Construct initial condition using tanh regularization.
        delta = 10* ones(1,2);
        iu = (uL(1)-uL(2))/2*( tanh(delta(1)*(x-jumpLocation(1))) ...
            - tanh(delta(2)*(x-jumpLocation(2))) ) + uL(2);
end

%% Solve the wave equation.
tstart = tic;
% Numerical integration for wave solution using pseudo-spectral method.
ko = [0:N/2-1  0  -N/2+1:-1]/scale;     % for odd-order derivative.
ke = [0:N/2  -N/2+1:-1]/scale;          % for even-order derivative.
% k = [0: (N-1)/2,  -(N-1)/2 : -1]/scale;
ik2 = transpose(1i*ke.*ke.*sign(ko));
ik = transpose(1i*ko);

% ik2 = transpose(1i * k.*k.*sign(k));
% ik = transpose(1i * k);
q = -HilbertSign*ik2;
odefunc = @(t, u_fourier) ...
    -ik.*fft( ifft(u_fourier).^2 ) + ik2 .* u_fourier;    % u_t + 2u u_x + Hu_{xx}
% odefunc = @(t, U_fourier) ...
%     -0.5*ik.*exp(q*t).*fft( ifft( exp(-q*t).*U_fourier ).^2  ); % u_t + 2u u_x + Hu_{xx}
    
% [t, uFourierMatrix] = ode45(ODEFUNC_1, 0:0.5:T, transpose(fft(iu)) );
% u = ifft(uFourierMatrix, [], 2);  % Compute  solution in the real space.
dt = 0.20; 
tspan = linspace(0, T, T/dt+1);
[t, UFourierMatrix] = ode45(odefunc, tspan, transpose(fft(iu)) );
% u = zeros(length(t), N);
% for i = 1:length(t)
%     u(i,:) = ifft( exp(-transpose(q)*t(i)).*UFourierMatrix(i,:) );
% end
u = ifft(UFourierMatrix, [], 2);

telapsed = toc(tstart)

%% Plot the solution


% Plot the modulation solution on top of direct numerical solution.
% close all
% figure('Position', [100, 75, 1000, 900])
% fs = 14;
% tt = t(end);
% % Generate the modulation solutio (given by Whitham theory).
% xx = linspace(0, tt, 1e4+1);
% aa = (0.5/tt)*xx;
% bb = 0.5*uL(1);
% cc = 0;
% theta = (bb - aa).*xx - (bb.*bb - aa.*aa)*tt;
% uu = 2*(bb-aa).^2 ./ (aa + bb - 2*cc - 2*sqrt((aa-cc).*(bb-cc)).*cos(theta)) + 2*cc;
% v_soliton = uL(1);
% v_linear = uL(2);
% 
% tplot = [250, 350, 450, 500];
% 
% for i = 1:length(tplot)
%     subplot(2,2,i)
%     plot(x, u(t == tplot(i), :), '-')
%     grid on, hold on
%     plot(v_soliton*ones(1,2)*tplot(i), [uL(2), uL(1)],  'k-', ...
%         v_linear*ones(1,2)*tplot(i), [uL(2), uL(1)], 'k-', 'LineWidth', 2)
%     title(['$a(x, t = ', num2str(tplot(i)), ')$'], ...
%         'Interpreter', 'latex', 'FontSize', fs)
% end
% legend('  Direct numerical solution', '  Modulation solution', ...
%     'FontSize', 16, 'Location', 'best')