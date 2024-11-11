function allScenarioStruc = loadScenariosData(allRsoFixPositions)
epochJday = allRsoFixPositions(1,1,4) - (60/86400); % one min back
nbRso = size(allRsoFixPositions, 2);
epochMjDay = epochJday - 2400000.5; % modify julian Day
% a script for loading observers' positions
% buiding 1100 walker constellation
nbSce = 0;
step = 60/86400;
inclinations = 0:10:100;
fovs = [90 - 0.5, 90 - 1.5, 90 - 2.5];
allScenarioStruc(3300) = struct('nbSbObs', [], 'nbPlane', [], ...
    'nbSatPp', [], 'phase', [], ...
    'inclination', [], 'fov', [], ...
    'sbPos', [], 'fullAccessArray', []);
for f = 1:3
    for i = 1:11
        for j = 1:10
            for k = 1:10
                nbSce = nbSce + 1;
                nbPlane = j;
                nbSatPerPlane = k;
                nbSb = nbPlane*nbSatPerPlane;
                phase = floor(nbPlane/2);
                allScenarioStruc(nbSce).nbSbObs = nbSb;
                allScenarioStruc(nbSce).nbPlane = nbPlane;
                allScenarioStruc(nbSce).nbSatPp = nbSatPerPlane;
                allScenarioStruc(nbSce).phase = phase;
                allScenarioStruc(nbSce).inclination = inclinations(i);
                allScenarioStruc(nbSce).fov = fovs(f);
                allScenarioStruc(nbSce).sbPos = zeros(1440, nbSb, 4);
                allScenarioStruc(nbSce).fullAccessArray = ...
                    uint8( zeros(1440, nbSb, nbRso) );
            end
        end
    end
end
delete(gcp('nocreate'));
parpool;
clc;
D = parallel.pool.DataQueue;
afterEach(D, @ParforProgressbar1100);
parfor i = 1 : 1100
    incl = allScenarioStruc(i).inclination;
    refCoe = [500000+6378137, 0.00001, incl, 0, 0, 0];
    nbPlane = allScenarioStruc(i).nbPlane;
    nbSatPerPlane = allScenarioStruc(i).nbSatPp;
    nbSb = allScenarioStruc(i).nbSbObs;
    phase = allScenarioStruc(i).phase;
    allSbCoe = zeros(nbSb,6); %#ok
    allSbCoe = makeWalker(nbPlane, nbSatPerPlane, phase, refCoe);
    allSbEci = zeros(nbSb,6); %#ok
    allSbEci = multiCoeToEciW(allSbCoe);
    allSbObsFixPositions = zeros(1440, nbSb, 4);
    % propagation for a day
    for tt = 1:1440
        % in this loop the fist sample will be excluded - there is
        % a reson: epoch is a minute sooner!
        dt = tt * step;
        curretMjDay = epochMjDay + dt;
        curretJday = epochJday + dt;
        toEcef = hpFrameConversion( curretMjDay );
        for rr = 1:nbSb
            [positionOut, velocityOut] = simplePropagator.propagate(step * 86400, ...
                allSbEci(rr, 1:3), allSbEci(rr, 4:6) );
            allSbEci(rr, 1:3) = positionOut;
            allSbEci(rr, 4:6) = velocityOut;
            [fixedPos, ~] = toEcef.convertToFixed(positionOut, velocityOut);
            allSbObsFixPositions(tt, rr,1:3) = fixedPos;
            allSbObsFixPositions(tt, rr,4) = curretJday;
        end
    end
    allScenarioStruc(i).sbPos = allSbObsFixPositions;
    send(D, i);
end
indx = 0;
for j = 1101:3300
    if j > 2200
        indx = j - 2200;
    elseif j > 1100
        indx = j - 1100;
    else
        pass;
    end
    allScenarioStruc(j).sbPos = allScenarioStruc(indx).sbPos;
end
clc;
end