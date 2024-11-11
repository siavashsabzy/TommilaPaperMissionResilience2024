function hasAcc = hasAccess(PosObs, PosSV, mask)
assert(isa(PosObs,'double'));
assert(all(size(PosObs) == [1 3]));
assert(isa(PosSV,'double'));
assert(all(size(PosSV) == [1 3]));
assert(isa(mask,'double'));
hasAcc = uint8(0);
R=PosSV-PosObs;
a = 6378137;
e = 8.1819190842622e-2;
% calculations:
b   = sqrt(a^2*(1-e^2));
ep  = sqrt((a^2-b^2)/b^2);
p   = sqrt(PosObs(1)^2+PosObs(2)^2);
th  = atan2(a*PosObs(3),b*p);
lon = atan2(PosObs(2),PosObs(1));
lat = atan((PosObs(3)+ep^2*b*(sin(th))^3)/(p-e^2*a*(cos(th))^3));
XYZ2ENU=[-sin(lon) cos(lon) 0;
    -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat);
    cos(lat)*cos(lon) cos(lat)*sin(lon)  sin(lat)];
ENU=XYZ2ENU*R';
Elevation=asin(ENU(3)/norm(ENU)) * 180 / pi;
if Elevation > mask
    hasAcc = uint8(1);
end
end