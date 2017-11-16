
date_start=datenum('03-Dec-2006 00:00:00');
%date_start=datenum('06-Feb-2007 00:00:00');
date_start=datenum('05-Jan-2007 00:00:00');
date_start=datenum('22-Feb-2007 00:00:00');
%date_end=date_start+5;
%date_end=datenum('20-Jan-2007 00:00:00')+34;
date_end=datenum('15-Apr-2007 00:00:00');


period_length=10;
for date=date_start:date_end
    close all
    [saved_obs_name,saved_mobs_name,saved_mod_name]=main_script_C(datestr(date),datestr(date+period_length),'nextsim.log',0,1,0);
    %plots_def_map({saved_obs_name,saved_mobs_name,saved_mod_name},{['obs_' num2str(date)],['mobs_' num2str(date)],['mod_' num2str(date)]},'arctic',1,'png',date+0.5*period_length);
    plots_def_map({saved_mod_name},{['mod_' num2str(date)]},'arctic',0,'png',date+0.5*period_length);
    
    %plots_def_map({saved_obs_name},{['obs_' num2str(date)]},'arctic',1,'png',date+0.5*period_length);
end