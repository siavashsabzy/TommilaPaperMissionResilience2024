clear all; %#ok for persistant vars
clc
%%
allRsoFixPositions = getRsoPostions();
allScenarioStruc = loadScenariosData(allRsoFixPositions);
sunFixPosArray = getSunFixPositions ( allRsoFixPositions(:,1,4) );
%%
clc;
delete(gcp('nocreate'));
parpool;
clc;
D2 = parallel.pool.DataQueue;
afterEach(D2, @ParforProgressbar3300);
nbScenario = size(allScenarioStruc,2);
parfor s = 1:nbScenario
    sbPos = allScenarioStruc(s).sbPos;
    mask = allScenarioStruc(s).fov;
    allScenarioStruc(s).fullAccessArray = ...
        getFullAccArray(sbPos, allRsoFixPositions, mask, sunFixPosArray);
    send(D2, s);
end
save('allScenarioStruc.mat', 'allScenarioStruc');
delete(gcp('nocreate'));
%%
gbObsPos = getGbObsPos();
gbObsAccessArray = getGbFullAccArray(gbObsPos, allRsoFixPositions, 10, sunFixPosArray);
clc;
nbMeoRso = 38;
nbGeoRso = 38;
nbRso = nbMeoRso + nbGeoRso;
nbScenario = size(allScenarioStruc,2);
allTotalLatSb = zeros(nbScenario, 38);
allTotalLatGb = zeros(nbScenario, 38);
meanTotalLat = zeros(nbScenario, 1);
maxLatGbInc = zeros(nbScenario, 1);
maxLatSb = zeros(nbScenario, 1);
maxLatMeo = zeros(nbScenario, 1);
maxLatGeo = zeros(nbScenario, 1);
meanTotalLatMeo = zeros(nbScenario, 1);
meanTotalLatGeo = zeros(nbScenario, 1);
MIT = zeros(nbScenario, 3); % meanTotalLat, incl, totalSat
MITmeo = zeros(nbScenario, 3); % meanTotalLat, incl, totalSat
MITgeo = zeros(nbScenario, 3); % meanTotalLat, incl, totalSat
MNF = zeros(nbScenario, 3); % maxLat, nSat, fov
MNP = zeros(nbScenario, 3); % maxLat, nSat, nbPlane
MMT = zeros(nbScenario, 3); % maxLat, nSat, nbPlane
gbInclude = true;
for s = 1:nbScenario
    accArray = allScenarioStruc(s).fullAccessArray;
    fov = (90 - allScenarioStruc(s).fov) * 2;
    nbPlane = allScenarioStruc(s).nbPlane;
    nbObser = size(accArray, 2);
    incl = allScenarioStruc(s).inclination;
    curMaxLatGbInc = 0;
    curMaxLatSb = 0;
    curMaxLatMeo = 0;
    curMaxLatGeo = 0;
    for r = 1:76
        latCouter = 0;
        latCoutersbOnly = 0;
        totalLatSb = 0;
        totalLatGbIncl = 0;
        for t = 1:1440
            seen = false;

            for o = 1:nbObser
                if accArray(t,o,r) == 1 || seen
                    seen = true;
                    break;
                end
            end

            if seen
                latCoutersbOnly = 0;
            else
                latCoutersbOnly = latCoutersbOnly + 1;
                totalLatSb = totalLatSb + 1;
                if latCoutersbOnly > curMaxLatMeo && r < 39
                    curMaxLatMeo = latCoutersbOnly;
                end
                if latCoutersbOnly > curMaxLatGeo && r > 38
                    curMaxLatGeo = latCoutersbOnly;
                end
                if latCoutersbOnly > curMaxLatSb
                    curMaxLatSb = latCoutersbOnly;
                end
            end

            if (gbInclude && ( max( gbObsAccessArray(t,:,r) ) > 0 ) )
                seen = true;
            end

            if seen
                latCouter = 0;
            else
                latCouter = latCouter + 1;
                totalLatGbIncl = totalLatGbIncl + 1;

                if latCouter > curMaxLatGbInc
                    curMaxLatGbInc = latCouter;
                end

            end
        end
        allTotalLatSb(s,r) = totalLatSb;
        allTotalLatGb(s,r) = totalLatGbIncl;
    end
    meanTotalLatMeo(s) = mean(allTotalLatSb(s,1:38));
    meanTotalLatGeo(s) = mean(allTotalLatSb(s,39:76));
    maxLatGbInc(s) = curMaxLatGbInc;
    maxLatSb(s) = curMaxLatSb;
    maxLatMeo(s) = curMaxLatMeo;
    maxLatGeo(s) = curMaxLatGeo;
    MITmeo(s,:) = [meanTotalLatMeo(s), incl, nbObser];
    MITgeo(s,:) = [meanTotalLatGeo(s), incl, nbObser];
    MNF(s,:) = [maxLatSb(s), nbObser, fov];
    MNP(s,:) = [maxLatGbInc(s), nbObser, nbPlane];
    MMT(s,:) = [maxLatMeo(s), maxLatGeo(s), nbObser];
end

%%
nbResi1 = 0;
nbResi2 = 0;
meanTotalLatSb = zeros(nbScenario, 1);
meanTotalLatGb = zeros(nbScenario, 1);
for s = 1:nbScenario
    if maxLatSb(s) < (1440 / 2)
        nbResi1 = nbResi1 + 1;
    end
    meanTotalLatSb(s) = mean(allTotalLatSb(s,:));
    meanTotalLatGb(s) = mean(allTotalLatGb(s,:));
    if meanTotalLatSb(s) < 1441
        nbResi2 = nbResi2 + 1;
    end
end
ind = 0;
ind2 = 0;
NRI = zeros(nbResi1, 3);
resi = zeros(nbResi1, 1);
NRI2 = zeros(nbResi2, 3);
resi2 = zeros(nbResi2, 1);
for s = 1:nbScenario
    nbObser = allScenarioStruc(s).nbSbObs;
    incl = allScenarioStruc(s).inclination;
    if maxLatSb(s) < (1440 / 2)
        ind = ind + 1;
        um = 1 - (maxLatGbInc(s) / 1440);
        ump = 1 - (maxLatSb(s) / 1440);
        resi(ind) = 1 - (um - ump);
        NRI(ind,:) = [nbObser, resi(ind), incl];
    end

    if meanTotalLatSb(s) < 1441
        ind2 = ind2 + 1;
        um2 = 1 - meanTotalLatGb(s) / 1440;
        ump2 = 1 - meanTotalLatSb(s) / 1440;
        resi2(ind2) = 1 - (um2 - ump2);
        NRI2(ind2,:) = [nbObser, resi2(ind2), incl];
    end

end
%%
figure()
scatter(MITmeo(:,2), MITmeo(:,1), 30, MITmeo(:,3),"filled","Marker","o", "LineWidth",1)
hold on
colorbar;

figure()
scatter(MITgeo(:,2), MITgeo(:,1), 30, MITgeo(:,3),"filled","Marker","o", "LineWidth",1)
hold on
colorbar;

MNF = sortrows(MNF,3,"descend");
figure()
scatter(MNF(:,2), MNF(:,1), 30, MNF(:,3),"filled","Marker","o", "LineWidth",1)
hold on
colorbar;

MNP = sortrows(MNP,3,"descend");
figure()
scatter(MNP(:,2), MNP(:,1), 30, MNP(:,3),"filled","Marker","o", "LineWidth",1)
hold on
colorbar;

figure()
scatter(MMT(:,2), MMT(:,1), 30, MMT(:,3),"filled","Marker","o", "LineWidth",1)
hold on
colorbar;

figure()
scatter(NRI(:,2), NRI(:,1), 30, NRI(:,3),"filled","Marker","o", "LineWidth",1)
hold on
colorbar;

figure()
scatter(NRI2(:,2), NRI2(:,1), 30, NRI2(:,3),"filled","Marker","o", "LineWidth",1)
hold on
colorbar;