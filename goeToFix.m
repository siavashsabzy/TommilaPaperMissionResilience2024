function pos=goeToFix(lla)
lat = lla(1) * pi / 180;
lon = lla(2) * pi / 180;
alt = lla(3);
% WGS84 ellipsoid constants:
a = 6378137;
e = 8.1819190842622e-2;
% intermediate calculation
% (prime vertical radius of curvature)
N = a ./ sqrt(1 - e^2 .* sin(lat).^2);
% results:
x = (N+alt) .* cos(lat) .* cos(lon);
y = (N+alt) .* cos(lat) .* sin(lon);
z = ((1-e^2) .* N + alt) .* sin(lat);
pos = [x,y,z];
return