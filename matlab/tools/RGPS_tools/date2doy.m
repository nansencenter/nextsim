function doy = date2doy(dayV,monthV,yearV)
% The DATE2DOY function converts a Gregorian date to the day of year. 
%
% Syntax:
%
%   DOY = DATE2DOY(DAY,MONTH,YEAR)
%
% Inputs:
%
%   DAY - Day of month, enumerated in [1,31]
%   MONTH - Month, enumerated in [1,12]
%   YEAR - Year
%
% Outputs:
%   
%   DOY - Day of year, enumerated in [1,366]
%
% DATE2DOY was written by:
%
% Ethan Miller (esmiller@uiuc.edu)
% Remote Sensing and Space Sciences Group (http://rsss.csl.uiuc.edu/)
% Department of Electrical and Computer Engineering
% University of Illinois at Urbana-Champaign
% Current:  2009 February 3 / Version 1.1 (vectorized)
% Original:  before 2007 December 2 / Version 1.0

if ((length(dayV)~=length(monthV)) | (length(dayV)~=length(yearV)) | ...
        (length(monthV)~=length(yearV)))
    error('date2doy: inputs must have the same dimensions');
    return;
end

doy = zeros(size(dayV));

for i = 1:length(dayV)
    day = dayV(i);
    month = monthV(i);
    year = yearV(i);
    
    daysmonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    
    if (((mod(year,4)==0) & (mod(year,100)~=0)) | (mod(year,400)==0))
        daysmonth(2)=daysmonth(2)+1;
    end
    if (month > 1)
        doy(i) = sum(daysmonth(1:month-1)) + day;
    else
        doy(i) = day;
    end
end