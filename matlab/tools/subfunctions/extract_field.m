function [field_tmp, field_plotted]=extract_field(field,data_out,dirname,step,simul_in,no_error)

  if ~exist('no_error','var'); no_error = 0; end

  simul = simul_in.simul;

%---------------------------
% We extract the data fields
%---------------------------
  if strcmp(field,'Freezing_Temperature')
     fld = 'SST';
     [field_tmp]=get_and_check(fld,data_out,dirname,step,no_error);

     fld = 'SSS';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step,no_error);
     
     field_tmp = field_tmp+0.055*field_tmp2;
     field_plotted = 'Ocean degrees above the freezing point';
  elseif strcmp(field,'Skin_temperature')
     fld = 'SST';
     [sst]=get_and_check(fld,data_out,dirname,step,no_error);

     fld = 'Tice_0';
     [tsurf]=get_and_check(fld,data_out,dirname,step,no_error);
     
     fld = 'Tsurf_thin_ice';
     [tsurf_thin]=get_and_check(fld,data_out,dirname,step,no_error);
     
     fld = 'Concentration';
     [conc]=get_and_check(fld,data_out,dirname,step,no_error);
     
     fld = 'Concentration_thin_ice';
     [conc_thin]=get_and_check(fld,data_out,dirname,step,no_error);
     
     field_tmp = conc.*tsurf + conc_thin.*tsurf_thin + (1-conc-conc_thin).*sst;
     field_plotted = 'Ocean degrees above the freezing point';
%   elseif strcmp(field,'Qatm')
%      fld = 'Qsw';
%      [Qsw]=get_and_check(fld,data_out,dirname,step,no_error);
% 
%      fld = 'Qlw';
%      [Qlw]=get_and_check(fld,data_out,dirname,step,no_error);
%      
%      fld = 'Qsh';
%      [Qsh]=get_and_check(fld,data_out,dirname,step,no_error);
%      
%      fld = 'Qlh';
%      [Qlh]=get_and_check(fld,data_out,dirname,step,no_error);
%           
%      field_tmp = Qsw + Qlw + Qsh + Qlh;
%      field_plotted = 'Total atmospheric heat flux';
  elseif strcmp(field,'Ridged_volume_per_area')
     fld = 'Ridge_ratio';
     [field_tmp]=get_and_check(fld,data_out,dirname,step,no_error);

     fld = 'Thickness';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step,no_error);
     
     field_tmp = field_tmp.*field_tmp2;
     field_plotted = 'Ridged volume per area';
  elseif strcmp(field,'Total_thickness')
     fld = 'Thickness';
     [field_tmp]=get_and_check(fld,data_out,dirname,step,no_error);

     fld = 'Thin_ice';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step,no_error);
     
     field_tmp = field_tmp+field_tmp2;
     field_plotted='Thick and thin ice thickness';
  elseif strcmp(field,'Ice_thickness')
     fld = 'Thickness';
     [field_tmp1]=get_and_check(fld,data_out,dirname,step,no_error);

     fld = 'Concentration';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step,no_error);
     
     field_tmp = field_tmp1./field_tmp2;
     %field_tmp(field_tmp2<=0.01) = 0.;
     field_tmp(field_tmp1==0.0) = 0.;
     field_plotted='Thick and thin ice thickness';
  elseif strcmp(field,'Thin_ice_thickness')
     fld = 'Thin_ice';
     [field_tmp]=get_and_check(fld,data_out,dirname,step,no_error);

     fld = 'Concentration_thin_ice';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step,no_error);
     
     field_tmp = field_tmp./field_tmp2;
     field_tmp(field_tmp2<=0.01) = 0.;
     field_plotted='Thick and thin ice thickness';
  elseif strcmp(field,'Total_concentration')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step,no_error);

     fld = 'Concentration_thin_ice';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step,no_error);     
     field_tmp = field_tmp+field_tmp2;
     field_plotted='Thick and thin ice concentration';
  elseif strcmp(field,'Critical_external_stress')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step,no_error);

     fld = 'Thickness';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step,no_error);
     
     fld = 'Ridge_ratio';
     [field_tmp3]=get_and_check(fld,data_out,dirname,step,no_error);
     
     
     field_tmp = simul.cfix*field_tmp2*simul.scale_coef.*exp(simul.ridging_exponent*(1.-field_tmp));
     field_plotted='Critical_external_stress';   
  elseif strcmp(field,'New_critical_external_stress')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step,no_error);

     fld = 'Thickness';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step,no_error);
     
     fld = 'Ridge_ratio';
     [field_tmp3]=get_and_check(fld,data_out,dirname,step,no_error);
     
     
     field_tmp = (simul.cfix*(1-field_tmp3)+6*simul.cfix*field_tmp3).*field_tmp2*simul.scale_coef.*exp(simul.ridging_exponent*(1.-field_tmp));
     field_plotted='New_critical_external_stress';   
  else
    [field_tmp]=get_and_check(field,data_out,dirname,step,no_error);
    field_plotted=field;
  end
end
