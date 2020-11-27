%BO  Solve the Benjamin-Ono equation
%   u_t + alpha*u u_x + Hu_{xx} = 0
%
% BO implements a numerical scheme for solving the Benjamin-Ono (BO)
% equation. The script is a complete solver of the BO equation in that it
% starts with setting up a uniform mesh on which the numerical solution
% will be computed, discretize the solution in space and solve the resulted
% system of ODEs in time. Unlike the solver for KdV equation (see KDV),
% this equation cannot be solved by other methods such as the finite
% difference and finite element methods due to the appearance of the
% Hilbert transform. The Hilbert transform of function f = f(x) is defined
% as the convolution of a singular kernel with the function itself. Thus,
% although it is possible to approximate a Hilbert transform in a numerical
% basis, it is not trivial to integrate such approximations into a
% numerical scheme for the wave solution. For more explanation on the
% numerical strategy, type "help kdv".
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
% [3] <a
% href="https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.566.6982&rep=rep1&type=pdf"
% > Modulation solutions for the Benjamin–Ono equation </a>
%
% See also KDV, NLS, DNLS, COMPLEX2BO, BENJAMIN

%% Space-time mesh
HilbertSign = 1; % This variable corresponds how we define the Hilbert transform.
a = -pi; b = pi; 
L = b - a; 
N = 2^7+1; h = L/N;
scale = L/(2*pi);
if mod(N, 2) == 0
    x = (a : h : b - h)';
    k = [0:N/2-1, 0, -N/2+1 : -1]'/scale;
else
    x = (a + h/2 : h : b -h/2)';
    k = [0:(N-1)/2, -(N-1)/2: -1]'/scale;
end

alpha = 2;
tmax = 2;
% Initial condition: We use the exact 1-phase periodic solution or the
% algebraic 1-soliton solution. These solutions can be found in [3].
p = -1; r = 1; s = 2;
t0 = 0;  theta0 = (r - s)*x - (r^2 - s^2)*t0;
iu = p + r - s + 2*(s-r)*(s-p - sqrt((s-p)*(r-p))*cos(theta0))...
    ./ (r+s-2*p - 2*sqrt((s-p)*(r-p))*cos(theta0));

%% Time integration using pseudo-spectral method
tstart = tic;
%--------------------------------------
ik = 1i*k;
q = -HilbertSign*1i*k.*k.*sign(k);
tspan = linspace(0, tmax, 101);
odefunc = @(t, vFourier) ...
    -0.5*alpha*ik.*exp(q*t) .* fft( ifft( exp(-q*t).*vFourier ).^2  );
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-10, 'Vectorized', 'on');
[t, vFourier] = ode45(odefunc, tspan, fft(iu), options); 
u = ifft(exp(-t*(q.')).*vFourier, [], 2);  
%--------------------------------------
telapsed = toc(tstart); fprintf('Time processed: %f \n', telapsed);

%% Plot the solution
close all
[xx, tt] = meshgrid(x, t);
theta = (r-s)*xx - (r^2 - s^2)*tt;
uExact = p + r - s + 2*(s-r)*(s-p - sqrt((s-p)*(r-p))*cos(theta))...
    ./ (r+s-2*p - 2*sqrt((s-p)*(r-p))*cos(theta));
subplot(1,2,1), mesh(xx, tt, u); view(3)
xlabel('x'); ylabel('t'); zlabel('u');
subplot(1,2,2)
mesh(xx, tt, u - uExact), view(3);
xlabel('x'); ylabel('t'); zlabel('$u - u_{\mathrm{exact}}$', 'Interpreter', 'latex');

%% Different initial conditions
% Dispersive shock waves
%----------------------------------------------
% uLimits = [1.0,  0.0];      % Jump between the left limit and right limit
% jumpLocations = [j1, j2];
% delta = 10* ones(1,2);
% iu = (uLimits(1)-uLimits(2))/2*( tanh(delta(1)*(x-jumpLocations(1))) ...
%     - tanh(delta(2)*(x-jumpLocations(2))) ) + uLimits(2);

% 1-soliton solution
% r = 1; iu = 4*r./(1 + 4*r*r*(x - 2*r*t0).^2);
