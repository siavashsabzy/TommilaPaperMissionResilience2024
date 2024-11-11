% ------------------------------------------------------------------------------
%
%                           function days2mdh + function jday.m
%
%  this function converts the day of the year, days, to the equivalent
% julian date
%
%  author        : david vallado                  719-573-2600   22 jun 2002
%
%  inputs          description                    range / units
%    year        - year                           1900 .. 2100
%    days        - julian day of the year         1.0  .. 366.0
%
%  outputs       :
%    jd          - julian date                    days from 4713 bc
% -----------------------------------------------------------------------------
function julianDay = days2jDay ( yr,days)
% --------------- set up array of days in month  --------------
for i= 1 : 12
    lmonth(i) = 31;
    if i == 2
        lmonth(i)= 28;
    end
    if i == 4 || i == 6 || i == 9 || i == 11
        lmonth(i)= 30;
    end
end
dayofyr= floor(days );
% ----------------- find month and day of month ---------------
if rem(yr-1900,4) == 0
    lmonth(2)= 29;
end
i= 1;
inttemp= 0;
while ( dayofyr > inttemp + lmonth(i) ) && ( i < 12 )
    inttemp= inttemp + lmonth(i);
    i= i+1;
end
mon= i;
day= dayofyr - inttemp;
% ----------------- find hours minutes and seconds ------------
temp = (days - dayofyr )*24.0;
hr  = fix( temp );
temp = (temp-hr) * 60.0;
min = fix( temp );
sec = (temp-min) * 60.0;
jd = 367.0 * yr  ...
    - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )   ...
    + floor( 275 * mon / 9.0 ) ...
    + day + 1721013.5;   % use - 678987.0 to go to mjd directly
jdfrac = (sec + min * 60.0 + hr *3600.0) / 86400.0;
% check jdfrac
if jdfrac > 1.0
    jd = jd + floor(jdfrac);
    jdfrac = jdfrac - floor(jdfrac);
end
julianDay = jd + jdfrac;
end
