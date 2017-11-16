
years=['2003_12_04_2004_06_01';
'1999_11_01_2000_05_14';
'2001_11_05_2002_06_01';
'1997_11_02_1998_06_01';
'2002_12_27_2003_06_01';
'1998_10_28_1999_05_17';
'2004_11_10_2005_06_01';
'1996_11_07_1997_06_01';
'2000_11_04_2001_06_01';
'2005_11_29_2006_06_01';
'2007_12_01_2008_06_01';
'2006_12_03_2007_06_01']

for year=1:12
%parfor year=1996:2007
    % We load the motion field from the RGPS Lagragian dataset
    %RGPS_kwok_vel_def(year,'Motion');
    
    startdate=years(year,1:10)
    enddate=years(year,12:21)
    %We compute and save the deformations for the currently processed period
    RGPS_drift_deformation(startdate,enddate,1);
end;


% parfor year=1996:2007
%     % We load the motion field from the RGPS Lagragian dataset
%     RGPS_kwok_vel_def(year,'Motion');
%     
%     %We compute and save the deformations for the currently processed period
%     RGPS_drift_deformation(startdate,enddate,1);
% end;
