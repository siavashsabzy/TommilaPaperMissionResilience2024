classdef simplePropagator
    properties


    end
    methods (Static)
        function xdot = gravityJ2(  ~ , x )
            Re  = 6378.137;
            J2  = 1.08262668e-3;
            mu  = 398600.4415;
            r   = norm(x(1:3));
            xdot    = [ x(4);
                x(5);
                x(6);
                -mu*x(1)/r^3*(1+1.5*J2*(Re/r)^2*(1-5*(x(3)/r)^2));
                -mu*x(2)/r^3*(1+1.5*J2*(Re/r)^2*(1-5*(x(3)/r)^2));
                -mu*x(3)/r^3*(1+1.5*J2*(Re/r)^2*(3-5*(x(3)/r)^2))];
        end

        function [ t, y ] = rk4(  dydt, tspan, y0, n )
            y0 = transpose ( y0(:) );
            m = size ( y0, 2 );
            t = zeros ( n + 1, 1 );
            y = zeros ( n + 1, m );
            tfirst = tspan(1);
            tlast = tspan(2);
            dt = ( tlast - tfirst ) / n;
            t(1,1) = tspan(1);
            y(1,:) = y0(1,:);
            for i = 1 : n
                f1 = dydt ( t(i,1),            y(i,:) );
                f2 = dydt ( t(i,1) + dt / 2.0, y(i,:) + dt * transpose ( f1 ) / 2.0 );
                f3 = dydt ( t(i,1) + dt / 2.0, y(i,:) + dt * transpose ( f2 ) / 2.0 );
                f4 = dydt ( t(i,1) + dt,       y(i,:) + dt * transpose ( f3 ) );
                t(i+1,1) = t(i,1) + dt;
                y(i+1,:) = y(i,:) + dt * ( ...
                    transpose ( f1 ) ...
                    + 2.0 * transpose ( f2 ) ...
                    + 2.0 * transpose ( f3 ) ...
                    +       transpose ( f4 ) ) / 6.0;
            end
            return
        end

        function [positionOut, velocityOut] = propagate(secToGo, positionIn, velocityIn)
            positionOut = zeros(1 , 3 );
            velocityOut = zeros(1 , 3 );
            if secToGo <=0
                positionOut = positionIn;
                velocityOut = velocityIn;
                return;
            end
            step = ceil(secToGo / 6);
            ephvec_in = [positionIn', velocityIn'].*0.001;
            [ ~, y ] = simplePropagator.rk4( @simplePropagator.gravityJ2 , [0 secToGo], ephvec_in , step );
            positionOut(1,1:3) = y(end, 1:3).*1000.0;
            velocityOut(1,1:3) = y(end, 4:6).*1000.0;
        end

    end


end