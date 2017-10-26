function hs = Warren(time,element)

% Calculate snow thickness for the coordinates in 'element', for a given
% 'time', based on the Warren et al (1999) climatology
%
% hs = Warren('5-Mar-2008',element);
%

Jan = [28.01  0.1270 -1.1833 -0.1164 -0.0051  0.0243];
Feb = [30.28  0.1056 -0.5908 -0.0263 -0.0049  0.0044];
Mar = [33.89  0.5486 -0.1996  0.0280  0.0216 -0.0176];
Apr = [36.80  0.4046 -0.4005  0.0256  0.0024 -0.0641];
May = [36.93  0.0214 -1.1795 -0.1076 -0.0244 -0.0142];
Jun = [36.59  0.7021 -1.4819 -0.1195 -0.0009 -0.0603];
Jul = [11.02  0.3008 -1.2591 -0.0811 -0.0043 -0.0959];
Aug = [ 4.64  0.3100 -0.6350 -0.0655  0.0059 -0.0005];
Sep = [15.81  0.2119 -1.0292 -0.0868 -0.0177 -0.0723];
Oct = [22.66  0.3594 -1.3483 -0.1063  0.0051 -0.0577];
Nov = [25.57  0.1496 -1.4643 -0.1409 -0.0079 -0.0258];
Dec = [26.67 -0.1876 -1.4229 -0.1413 -0.0316 -0.0029];

% Select the current month
month = datestr(time,'mmm');
eval(['monmat = ', month,';'])
hs0 = monsel(monmat,element);

% We now assume that the climatology is correct for the 15th of each month and
% do a linear interpolation to the current date. This is discontinues on the
% scale of days.

% hs0 was for the 15th of 'month' - this is time0
timevec = datevec(time);
time0   = datenum([timevec(1:2) 15 zeros(1,3)]);

% Select the previous or next month for time1
if timevec(3) <= 15
    % Go back 17 days to be sure to get the previous month
    timevec = datevec(time0 - 17);
    time1   = datenum([timevec(1:2) 15 zeros(1,3)]);
    month   = datestr(time1, 'mmm');
else
    % Go forwards 18 days to be sure to get the next month
    timevec = datevec(time0 + 18);
    time1   = datenum([timevec(1:2) 15 zeros(1,3)]);
    month   = datestr(time1, 'mmm');
end
eval(['monmat = ', month,';'])
hs1 = monsel(monmat,element);

% Calculate the mean
dt    = abs(time0 - time1);
fac1  = abs(time - time0)/dt;
fac0  = abs(time - time1)/dt;

hs = (fac0*hs0 + fac1*hs1);

end

function hs = monsel(monmat,element)

H0 = monmat(1);
A  = monmat(2);
B  = monmat(3);
C  = monmat(4);
D  = monmat(5);
E  = monmat(6);

% Calculate the coordinates
m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
[lon, lat] = m_xy2ll(element.x/6378.273,element.y/6378.273);
r   = 90 - lat;
phi = lon;

x = r .* cosd(phi);
y = r .* sind(phi);

hs = 1e-2*max(0, H0 + A.*x + B.*y + C.*x.*y + D.*x.^2 + E.*y.^2);

end
