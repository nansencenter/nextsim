function [field_tmp]=extract_field(field,data_out,dirname,step,simul_in)
%---------------------------
% We extract the data fields
%---------------------------
  if strcmp(field,'Freezing_Temperature')
     fld = 'SST';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'SSS';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp+0.055*field_tmp2;
     field_plotted='Freezing_Temperature';
  elseif strcmp(field,'Ridged_volume_per_area')
     fld = 'Ridge_ratio';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Thickness';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp.*field_tmp2;
     field_plotted='Ridged_volume_per_area';
  elseif strcmp(field,'Total_thickness')
     fld = 'Thickness';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Thin_ice';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp+field_tmp2;
     field_plotted='Thick and thin ice thickness';
  elseif strcmp(field,'Ice_thickness')
     fld = 'Thickness';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Concentration';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp./field_tmp2;
     field_tmp(field_tmp2<=0.01) = 0.;
     field_plotted='Thick and thin ice thickness';
  elseif strcmp(field,'Total_concentration')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Concentration_thin_ice';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp+field_tmp2;
     field_plotted='Thick and thin ice concentration';
  elseif strcmp(field,'Critical_external_stress')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Thickness';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     fld = 'Ridge_ratio';
     [field_tmp3]=get_and_check(fld,data_out,dirname,step);
     
     
     field_tmp = simul_in.cfix*field_tmp2*simul_in.scale_coef.*exp(simul_in.ridging_exponent*(1.-field_tmp));
     field_plotted='Critical_external_stress';   
  elseif strcmp(field,'New_critical_external_stress')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Thickness';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     fld = 'Ridge_ratio';
     [field_tmp3]=get_and_check(fld,data_out,dirname,step);
     
     
     field_tmp = (simul_in.cfix*(1-field_tmp3)+6*simul_in.cfix*field_tmp3).*field_tmp2*simul_in.scale_coef.*exp(simul_in.ridging_exponent*(1.-field_tmp));
     field_plotted='New_critical_external_stress';   
  else
    [field_tmp]=get_and_check(field,data_out,dirname,step);
  end
end
