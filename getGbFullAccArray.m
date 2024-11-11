function fullAccArray = getGbFullAccArray(gbObsPos, svPosArr, mask, sunFixPosArray)
nbObs = size(gbObsPos, 1);
nbRso = size(svPosArr, 2);
fullAccArray = uint8( zeros(1440, nbObs, nbRso) );
for t = 1 : 1440
    for o = 1 : nbObs
        PosObs = gbObsPos(o,:);
        inSun = hasAccess_mex(PosObs, sunFixPosArray(t,:), -12);
        if inSun
            
        else
            for s = 1:nbRso
                PosSV = [svPosArr(t,s,1), svPosArr(t,s,2), svPosArr(t,s,3)];
                fullAccArray(t, o, s) = hasAccess_mex(PosObs, PosSV, mask);
            end
        end

    end
end
end