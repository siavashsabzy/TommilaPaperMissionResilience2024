function sats = makeWalker(nbPlane, nbSatPerPlane, phase, refCoe)
sats   = zeros(nbPlane*nbSatPerPlane,6);
a      = refCoe(1,1);
ecc    = refCoe(1,2);
inc    = refCoe(1,3);
raan_0 = refCoe(1,4);
aop    = refCoe(1,5);
nu_0   = refCoe(1,6);
nb = 0;
deg2rad = pi / 180.0;
for i = 1:nbPlane
    for j = 1:nbSatPerPlane
        raan = raan_0 + i * 360 / nbPlane;
        if raan == 360
            raan = 0;
        end
        nu = nu_0 + j * 360 / nbSatPerPlane + i * 360 * phase / (nbSatPerPlane * nbPlane);
        if nu >= 180
            nu = nu - 360;
        end
        nb = nb+1;
        sats(nb,:) = [a, ecc, inc * deg2rad, raan * deg2rad, aop * deg2rad, ...
            nu * deg2rad];
    end
end
end