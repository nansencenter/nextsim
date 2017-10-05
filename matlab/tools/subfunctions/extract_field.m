function [field_tmp, field_plotted]=extract_field(field,data_out,dirname,step,simul_in)

  simul = simul_in.simul;

%---------------------------
% We extract the data fields
%---------------------------
  if strcmp(field,'Freezing_Temperature')
     fld = 'SST';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'SSS';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp+0.055*field_tmp2;
     field_plotted = 'Ocean degrees above the freezing point';
  elseif strcmp(field,'Skin_temperature')
     fld = 'SST';
     [sst]=get_and_check(fld,data_out,dirname,step);

     fld = 'Tice_0';
     [tsurf]=get_and_check(fld,data_out,dirname,step);

     try % Try to read the thin ice 
        fld = 'Tsurf_thin_ice';
        [tsurf_thin]=get_and_check(fld,data_out,dirname,step,true);
     
        fld = 'Concentration_thin_ice';
        [conc_thin]=get_and_check(fld,data_out,dirname,step,true);
     catch ERR
         tsurf_thin = 0;
         conc_thin  = 0;
     end
     
     fld = 'Concentration';
     [conc]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = conc.*tsurf + conc_thin.*tsurf_thin + (1-conc-conc_thin).*sst;
     field_plotted = 'Ocean degrees above the freezing point';
%   elseif strcmp(field,'Qatm')
%      fld = 'Qsw';
%      [Qsw]=get_and_check(fld,data_out,dirname,step);
% 
%      fld = 'Qlw';
%      [Qlw]=get_and_check(fld,data_out,dirname,step);
%      
%      fld = 'Qsh';
%      [Qsh]=get_and_check(fld,data_out,dirname,step);
%      
%      fld = 'Qlh';
%      [Qlh]=get_and_check(fld,data_out,dirname,step);
%           
%      field_tmp = Qsw + Qlw + Qsh + Qlh;
%      field_plotted = 'Total atmospheric heat flux';
  elseif strcmp(field,'Ridged_volume_per_area')
     fld = 'Ridge_ratio';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Thickness';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp.*field_tmp2;
     field_plotted = 'Ridged volume per area';
  elseif strcmp(field,'Total_thickness')
     fld = 'Thickness';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);
     
     try
        fld = 'Thin_ice';
        [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     catch
         field_tmp2=0.;
     end
     field_tmp = field_tmp+field_tmp2; 
     field_plotted='Thick and thin ice thickness';
  elseif strcmp(field,'Ice_thickness')
     fld = 'Thickness';
     [field_tmp1]=get_and_check(fld,data_out,dirname,step);

     fld = 'Concentration';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp1./field_tmp2;
     %field_tmp(field_tmp2<=0.01) = 0.;
     field_tmp(field_tmp1==0.0) = 0.;
     field_plotted='Thick and thin ice thickness';
  elseif strcmp(field,'Thin_ice_thickness')
     fld = 'Thin_ice';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Concentration_thin_ice';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp./field_tmp2;
     field_tmp(field_tmp2<=0.01) = 0.;
     field_plotted='Thick and thin ice thickness';
  elseif strcmp(field,'Total_concentration')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     try
        fld = 'Concentration_thin_ice';
        [field_tmp2]=get_and_check(fld,data_out,dirname,step);     
     catch
         field_tmp2=0.;
     end
     field_tmp = field_tmp+field_tmp2;
     
     field_plotted='Thick and thin ice concentration';
  elseif strcmp(field,'Lead_fraction')
      fld = 'Concentration';
      [field_tmp]=get_and_check(fld,data_out,dirname,step);   
      
      field_tmp = (1.-field_tmp);
      field_plotted='Lead_fraction';
  elseif strcmp(field,'Critical_external_stress')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Thickness';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     fld = 'Ridge_ratio';
     [field_tmp3]=get_and_check(fld,data_out,dirname,step);
     
     
     field_tmp = simul.cfix*field_tmp2*simul.scale_coef.*exp(simul.ridging_exponent*(1.-field_tmp));
     field_plotted='Critical_external_stress';   
  elseif strcmp(field,'sat_concentration')
     a_param=1.2;
     b_param=19.5;

     fld = 'Thickness';
     [field_tmp]=get_and_check(fld,data_out,dirname,step)/10;

     
     fld = 'Concentration';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);

     field_tmp = field_tmp./field_tmp2;
     field_tmp(field_tmp2<=0.01) = 0.;
     
     SIC_SIT=1-abs(min(0.,(field_tmp-a_param)/a_param)).^b_param;
     
     field_tmp=field_tmp2.*SIC_SIT;
     
     fld = 'Tice_2';
     [field_tmp3]=get_and_check(fld,data_out,dirname,step);  
          
     field_tmp(field_tmp3>-2) = field_tmp(field_tmp3>-2) - 0.11 * (2 + field_tmp3(field_tmp3>-2));
     field_plotted='sat_concentration';        
     elseif strcmp(field,'melt_ponds')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Tice_2';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);  
          
     field_tmp = 0.11 * (2 + field_tmp2);
     field_tmp(field_tmp2<-2)=0;
     field_plotted='Melt pond fraction';        
     
  elseif strcmp(field,'New_critical_external_stress')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Thickness';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     fld = 'Ridge_ratio';
     [field_tmp3]=get_and_check(fld,data_out,dirname,step);
     
     
     field_tmp = (simul.cfix*(1-field_tmp3)+6*simul.cfix*field_tmp3).*field_tmp2*simul.scale_coef.*exp(simul.ridging_exponent*(1.-field_tmp));
     field_plotted='New_critical_external_stress';   
  else
    [field_tmp]=get_and_check(field,data_out,dirname,step);
    field_plotted=field;
  end
end
