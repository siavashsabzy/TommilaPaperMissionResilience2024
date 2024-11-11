function allRsoFixPositions = getRsoPostions()
activeFileId = fopen('active.txt');
tline = fgetl(activeFileId);
%% preallocation for faster search
% initializing with redundant zero rows.
meoCoe = zeros(50, 8); % including revNumber for later filtering.
geoCoe = zeros(1000,8); % including revNumber for later filtering.
% each row --> [sma, ecc, inc, raan, aop, meanAn, jdayEpoch, revNum]
%%
earthMu = 3.986005e14;      % [m^3/s^2]
deg2rad = pi / 180.0;         %  0.01745329251994330;  % [deg/rad]
xpdotp  = (2.0*pi) / 86400.0;   %  [rev/day] -> [rad/sec]
mm = 0;
ecc = 0;
lineNb = 0;
nbMeo = 0;
nbGeo = 0;
totalS = 0;
while(tline~=-1)
    lineNb = str2double(tline(1:2));
    if lineNb == 2
        totalS = totalS + 1;
        tline(26) = '.';
        for j = 27:33
            if (tline(j) == ' ')
                tline(j) = '0';
            end
        end
        mm = str2double( tline(52:63) );
        ecc = str2double( tline(26:33) );
        if ( 1440/mm <= 800.0 && 1440/mm >= 600.0 ) && (ecc < 0.25) ...
                && ( strcmp(rsoName(1:7), 'NAVSTAR') ) ||...
                ( ( mm <= 1.01 && mm >= 0.99 ) && (ecc < 0.01) )
            epochY = str2double( firstLine(19:20) ) + 2000;
            epochD = str2double( firstLine(21:32) );
            revNb = str2double( tline(64:68) );
            if mm <= 1.01
                nbGeo = nbGeo + 1;
                geoCoe(nbGeo, 1) = ( earthMu / ( ( mm * xpdotp ) ^ 2 ) ) ^ (1/3);
                geoCoe(nbGeo, 2) = ecc;
                geoCoe(nbGeo, 3) = str2double( tline(8:16) ) * deg2rad; % inc
                geoCoe(nbGeo, 4) = str2double( tline(17:25) ) * deg2rad; % raan
                geoCoe(nbGeo, 5) = str2double( tline(34:42) ) * deg2rad; % aop
                geoCoe(nbGeo, 6) = str2double( tline(43:51) ) * deg2rad; % ma
                geoCoe(nbGeo, 7) = days2jDay( epochY,epochD);
                geoCoe(nbGeo, 8) = revNb;
            else
                nbMeo = nbMeo + 1;
                meoCoe(nbMeo, 1) = ( earthMu / ( ( mm * xpdotp ) ^ 2 ) ) ^ (1/3);
                meoCoe(nbMeo, 2) = ecc;
                meoCoe(nbMeo, 3) = str2double( tline(8:16) ) * deg2rad; % inc
                meoCoe(nbMeo, 4) = str2double( tline(17:25) ) * deg2rad; % raan
                meoCoe(nbMeo, 5) = str2double( tline(34:42) ) * deg2rad; % aop
                meoCoe(nbMeo, 6) = str2double( tline(43:51) ) * deg2rad; % ma
                meoCoe(nbMeo, 7) = days2jDay( epochY,epochD);
                meoCoe(nbMeo, 8) = revNb;
            end

        end

    elseif lineNb == 1
        firstLine = tline;
    else
        rsoName = tline;
    end
    tline = fgetl(activeFileId);
end
%%
tempArr = zeros(nbGeo, 8); %#ok
newMeoCoe = zeros(nbMeo, 7); %#ok
newGeoCoe = zeros(nbMeo, 7); %#ok same rows as meo
allRso = zeros(nbMeo*2, 7); %#ok
geoFilteringStep = ceil(nbGeo / nbMeo);
tempArr = sortrows(geoCoe(1:nbGeo,:), 8);
newMeoCoe = meoCoe(1:nbMeo , 1:7); 
newGeoCoe = tempArr(1:geoFilteringStep:end,1:7);
allRso = [newMeoCoe; newGeoCoe];
clearvars -except allRso
%% propagation for a day with a minute step time
% find the most fresh coe
nbRso = size(allRso,1);
allRso = multiCoeToEci(allRso);
latestEpoch = max(allRso(:,7));
allRsoEci = zeros(nbRso, 6); % leaving out the last column
% propagate all RSOs to epoch
for j = 1:nbRso
    secToGo = ( latestEpoch - allRso(j, 7) ) * 86400.0;
    [positionOut, velocityOut] = simplePropagator.propagate(secToGo, ...
        allRso(j, 1:3), allRso(j, 4:6) );
    allRsoEci(j, 1:3) = positionOut;
    allRsoEci(j, 4:6) = velocityOut;
end
% propagation of all RSOs for a day
latestMjDay = latestEpoch - 2400000.5; % modify julian Day
step = 60 / 86400; % a minute
% no need for recording velocities
% allPositiosEixed = zeros(nbRso, 3, 1440); % an universal Array
preRsoEci = allRsoEci;
allRsoFixPositions = zeros(1440, nbRso, 4);
for tt = 1:1440
    % in this loop the fist sample will be excluded - for no reason :)
    dt = tt * step; 
    curretMjDay = latestMjDay + dt;
    curretJday = latestEpoch + dt;
    toEcef = hpFrameConversion( curretMjDay );
    for rr = 1:nbRso
        [positionOut, velocityOut] = simplePropagator.propagate(step * 86400, ...
            preRsoEci(rr, 1:3), preRsoEci(rr, 4:6) );
        preRsoEci(rr, 1:3) = positionOut;
        preRsoEci(rr, 4:6) = velocityOut;
        [fixedPos, ~] = toEcef.convertToFixed(positionOut, velocityOut);
        allRsoFixPositions(tt,rr,1:3) = fixedPos;
        allRsoFixPositions(tt,rr,4) = curretJday;
    end
end

end