function [initialSolution, initialVelocity, initialAcceleration] = ...
    getInitialCondition(initialCondition)

switch initialCondition
    case 'kink'
        t0 = 0;
        x0 = 0;
        c = 0.5*sqrt(3);
        c2 = c*c;
        m = sqrt(1 - c2);
        m2 = 1 - c2;
        initialSolution = @(x) 4 * atan( exp( (x - x0 - c*t0)./m) );
        initialVelocity = @(x) -2*c/m*sech( (x - x0 - c*t0)./m );
        initialAcceleration = @(x) -2*c2/m2...
            *sech((x - x0 - c*t0)./m).*tanh((x-x0-c*t0)./m);
                
    case 'antikink'
        t0 = 0;
        x0 = 0;
        c = 0.5*sqrt(3);
        c2 = c*c;
        m = sqrt(1 - c2);
        m2 = 1 - c2;
        initialSolution = @(x) 4 * atan( exp( -(x - x0 - c*t0)./m) );
        initialVelocity = @(x) 2*c/m*sech( (x - x0 - c*t0)./m );
        initialAcceleration = @(x) 2*c2/m2...
            *sech((x - x0 - c*t0)./m).*tanh((x - x0 - c*t0)./m);
        
    case 'kink-kink'
        t0 = 0;
        t1 = 10;
        x0 = 0;
        c = 0.5*sqrt(3);
        m = sqrt(1 - c*c);
        initialSolution = @(x) 4.*atan(c.*sech(c.*m.^(-1).*(t0+(-1).*t1)).*sinh(m.^(-1).*(x+(-1) ...
            .*x0)));
        initialVelocity = @(x) (-4).*c.^2.*sech(c.*m.^(-1).*(t0+(-1).*t1)).*sinh(m.^(-1).*(x+(-1) ...
            .*x0)).*(m+c.^2.*m.*sech(c.*m.^(-1).*(t0+(-1).*t1)).^2.*sinh(m.^( ...
            -1).*(x+(-1).*x0)).^2).^(-1).*tanh(c.*m.^(-1).*(t0+(-1).*t1));
        initialAcceleration = @(x) (-2).*c.^3.*((-1)+c.^2).^(-1).*((-3)+c.^2+cosh(2.*c.*m.^(-1).*(t0+ ...
            (-1).*t1))+(-1).*c.^2.*cosh(2.*m.^(-1).*(x+(-1).*x0))).*sech(c.* ...
            m.^(-1).*(t0+(-1).*t1)).^3.*sinh(m.^(-1).*(x+(-1).*x0)).*(1+c.^2.* ...
            sech(c.*m.^(-1).*(t0+(-1).*t1)).^2.*sinh(m.^(-1).*(x+(-1).*x0)) ...
            .^2).^(-2);
        
    case 'kink-antikink'
        t0 = 0;
        t1 = 10;
        x0 = 0;
        c = 0.5*sqrt(3);
        m = sqrt(1 - c*c);
        initialSolution = @(x) 4.*atan(c.^(-1).*sech(m.^(-1).*(x+(-1).*x0)).*sinh(c.*m.^(-1).*( ...
            t0+(-1).*t1)));
        initialVelocity = @(x) 4.*cosh(c.*m.^(-1).*(t0+(-1).*t1)).*sech(m.^(-1).*(x+(-1).*x0)).*( ...
            m+c.^(-2).*m.*sech(m.^(-1).*(x+(-1).*x0)).^2.*sinh(c.*m.^(-1).*( ...
            t0+(-1).*t1)).^2).^(-1);
        initialAcceleration = @(x) (-2).*c.^3.*((-1)+c.^2).^(-1).*((-3)+c.^2+(-1).*cosh(2.*c.*m.^(-1) ...
            .*(t0+(-1).*t1))+c.^2.*cosh(2.*m.^(-1).*(x+(-1).*x0))).*sech(m.^( ...
            -1).*(x+(-1).*x0)).^3.*sinh(c.*m.^(-1).*(t0+(-1).*t1)).*(c.^2+ ...
            sech(m.^(-1).*(x+(-1).*x0)).^2.*sinh(c.*m.^(-1).*(t0+(-1).*t1)) ...
            .^2).^(-2);
        
    case 'standing-breather'
        t0 = 0;
        t1 = 0;
        x0 = 0;
        omega = 2*sqrt(2)/3;
        m = sqrt(1 - omega*omega);
        m2 = m*m;
        m3 = m2*m;
        omega2 = omega*omega;
        omega3 = omega2*omega;
        initialSolution = @(x) 4 * atan( m/omega * sin(omega*(t0 - t1))...
            ./cosh(m*(x - x0) ) );
        initialVelocity = @(x) 4*m*cos(omega*(t0-t1)).*sech(m*(x-x0))...
            ./( 1 + m2/omega2*sech(m*(x-x0)).^2 .* sin(omega*(t0-t1)).^2 );
        initialAcceleration = @(x) -4*omega3*sin(omega*(t0-t1)).*...
            (m*omega2*sech(m*(x-x0))+(2*m3*cos(omega*(t0-t1))^2 ...
            + m3*sin(omega*(t0-t1))^2)*sech(m*(x-x0)).^3)./(( omega2 + m2*sin(omega*(t0-t1))^2 ...
            .*sech(m*(x-x0)).^2 ).^2);
        
    case 'moving breather'
        t0 = 0;
        x0 = -5;
        omega = 2*sqrt(2)/3;
        m = sqrt(1 - omega*omega);
        c = 0.5*sqrt(3);
        initialSolution = @(x) 4.*atan(m.*omega.^(-1).*sech((1+(-1).*c.^2).^(-1/2).*m.*(c.*t0+( ...
            -1).*x+x0)).*sin((1+(-1).*c.^2).^(-1/2).*omega.*(t0+c.*((-1).*x+ ...
            x0))));
        initialVelocity = @(x) 8.*(1+(-1).*c.^2).^(-1/2).*m.*omega.*cosh((1+(-1).*c.^2).^(-1/2).* ...
            m.*(c.*t0+(-1).*x+x0)).*((-1)+m.^2.*cos(2.*(1+(-1).*c.^2).^(-1/2) ...
            .*omega.*(t0+(-1).*c.*x+c.*x0))+((-1)+m.^2).*cosh(2.*(1+(-1).* ...
            c.^2).^(-1/2).*m.*(c.*t0+(-1).*x+x0))).^(-1).*((-1).*omega.*cos(( ...
            1+(-1).*c.^2).^(-1/2).*omega.*(t0+c.*((-1).*x+x0)))+c.*m.*sin((1+( ...
            -1).*c.^2).^(-1/2).*omega.*(t0+c.*((-1).*x+x0))).*tanh((1+(-1).* ...
            c.^2).^(-1/2).*m.*(c.*t0+(-1).*x+x0)));
        initialAcceleration = @(x) 4.*((-1)+c.^2).^(-1).*omega.*(omega.^2+m.^2.*sech((1+(-1).*c.^2) ...
            .^(-1/2).*m.*(c.*t0+(-1).*x+x0)).^2.*sin((1+(-1).*c.^2).^(-1/2).* ...
            omega.*(t0+c.*((-1).*x+x0))).^2).^(-2).*(c.^2.*m.^5.*sech((1+(-1) ...
            .*c.^2).^(-1/2).*m.*(c.*t0+(-1).*x+x0)).^5.*sin((1+(-1).*c.^2).^( ...
            -1/2).*omega.*(t0+c.*((-1).*x+x0))).^3+m.*omega.^2.*sech((1+(-1).* ...
            c.^2).^(-1/2).*m.*(c.*t0+(-1).*x+x0)).*(omega.^2.*sin((1+(-1).* ...
            c.^2).^(-1/2).*omega.*(t0+c.*((-1).*x+x0)))+(-2).*c.*m.*omega.* ...
            cos((1+(-1).*c.^2).^(-1/2).*omega.*(t0+c.*((-1).*x+x0))).*tanh((1+ ...
            (-1).*c.^2).^(-1/2).*m.*((-1).*c.*t0+x+(-1).*x0))+(-1).*c.^2.* ...
            m.^2.*sin((1+(-1).*c.^2).^(-1/2).*omega.*(t0+c.*((-1).*x+x0))).* ...
            tanh((1+(-1).*c.^2).^(-1/2).*m.*((-1).*c.*t0+x+(-1).*x0)).^2)+ ...
            m.^3.*sech((1+(-1).*c.^2).^(-1/2).*m.*(c.*t0+(-1).*x+x0)).^3.*sin( ...
            (1+(-1).*c.^2).^(-1/2).*omega.*(t0+c.*((-1).*x+x0))).*((1/2).* ...
            omega.^2.*(3+2.*c.^2+cos(2.*(1+(-1).*c.^2).^(-1/2).*omega.*(t0+( ...
            -1).*c.*x+c.*x0)))+(-1).*c.*m.*omega.*sin(2.*(1+(-1).*c.^2).^( ...
            -1/2).*omega.*(t0+c.*((-1).*x+x0))).*tanh((1+(-1).*c.^2).^(-1/2).* ...
            m.*(c.*t0+(-1).*x+x0))+c.^2.*m.^2.*sin((1+(-1).*c.^2).^(-1/2).* ...
            omega.*(t0+c.*((-1).*x+x0))).^2.*tanh((1+(-1).*c.^2).^(-1/2).*m.*( ...
            c.*t0+(-1).*x+x0)).^2));
        
    case 'ring-shaped kink'
        t0 = 0;
        r0 = 10;
        c = sqrt(3)/2;
        c2 = c*c;
        m = sqrt(1 - c2);
        m2 = 1 - c2;
        initialSolution = @(r) 4 * atan( exp( (r - r0 - c*t0)./m) );
        initialVelocity = @(r) -2*c/m*sech( (r - r0 - c*t0)./m );
        initialAcceleration = @(r) -2*c2/m2...
            *sech((r - r0 - c*t0)./m).*tanh((r-r0-c*t0)./m) ...
            - 2*c/m * (sech((r - r0 - c*t0)./m)./r);
        
%     case 'ring-shaped antikink'
%         t0 = 0;
%         r0 = 60;
%         c = -sqrt(3)/2;
%         c2 = c*c;
%         m = sqrt(1 - c2);
%         m2 = 1 - c*c;
%         initialSolution = @(r) 4 * atan( exp( -(r - r0 - c*t0)./m ) );
%         initialVelocity = @(r) 2*c/m * sech( (r - r0 - c*t0)./m );
%         initialAcceleration = @(r) -2*c2 / m2 ...
%             *sech((
        
    case 'five ring-shaped kinks'
        c0 = 0.9 : -0.1 : 0.5;
        tReturn = 50;     % Returning time.
        r0 = tReturn * sqrt(1 - c0.*c0) ./ asin(c0);
        c2 = c0.*c0;
        m = sqrt(1 - c2);
        % m2 = 1 - c2;
        initialSolution = @(r) 4 * atan( exp( (r - r0(1))./m(1) ) ) ...
            + 4 * atan( exp( (r - r0(2))./m(2) ) ) ...
            + 4 * atan( exp( (r - r0(3))./m(3) ) ) ...
            + 4 * atan( exp( (r - r0(4))./m(4) ) ) ...
            + 4 * atan( exp( (r - r0(5))./m(5) ) );
        initialVelocity = @(r) -2*c0(1)/m(1) * sech( (r - r0(1))./m(1) ) ...
            -2*c0(2)/m(2) * sech( (r - r0(2))./m(2) ) ...
            -2*c0(3)/m(3) * sech( (r - r0(3))./m(3) ) ...
            -2*c0(4)/m(4) * sech( (r - r0(4))./m(4) ) ...
            -2*c0(5)/m(5) * sech( (r - r0(5))./m(5) );
        initialAcceleration = @(r) 0;
end