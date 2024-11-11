function epheci = multiCoeToEciW(ceos)

mu  = 3.986005e14;      % [m^3/s^2]
nbSats = size(ceos,1);
nbCol = size(ceos,2);
epheci = zeros(nbSats,nbCol);

%% MAIN ALGORITHM
for ou = 1:nbSats
    a = ceos(ou,1);
    e = ceos(ou,2);
    i = ceos(ou,3) ;
    OMEGA = ceos(ou,4);
    w  = ceos(ou,5);
    f = ceos(ou,6);
    % Handle special cases:
    
    % Circular equitorial orbit.
    if e == 0 && i == 0
        w     = 0;
        OMEGA = 0;
        
        % Circular inclined orbit.
    elseif e == 0
        w     = 0;
        
        % Elliptical equitorial.
    elseif i == 0
        OMEGA = 0;
    end
    
    % Compute the semi-latus rectum.
    p = a*(1-e^2);
    
    
    % Define the position vector in perifocal PQW coordinates.
    X(1,1) = p*cos(f)/(1+e*cos(f));
    X(2,1) = p*sin(f)/(1+e*cos(f));
    X(3,1) = 0;
    
    % Define the velocity vector in perifocal PQW coordinates.
    V(1,1) = -sqrt(mu/p)*sin(f);
    V(2,1) =  sqrt(mu/p)*(e+cos(f));
    V(3,1) =  0;
    
    % Define Transformation Matrix To IJK.
    Trans(1,1) =  cos(OMEGA)*cos(w)-sin(OMEGA)*sin(w)*cos(i);
    Trans(1,2) = -cos(OMEGA)*sin(w)-sin(OMEGA)*cos(w)*cos(i);
    Trans(1,3) =  sin(OMEGA)*sin(i);
    
    Trans(2,1) =  sin(OMEGA)*cos(w)+cos(OMEGA)*sin(w)*cos(i);
    Trans(2,2) = -sin(OMEGA)*sin(w)+cos(OMEGA)*cos(w)*cos(i);
    Trans(2,3) = -cos(OMEGA)*sin(i);
    
    Trans(3,1) =  sin(w)*sin(i);
    Trans(3,2) =  cos(w)*sin(i);
    Trans(3,3) =  cos(i);
    % Transform to the ECI coordinate system.
    X = Trans*X;
    V = Trans*V;
    epheci(ou,1:6) = [X',V'];
end