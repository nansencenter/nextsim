function [field_tmp, field_plotted]=extract_field(field,data_out,dirname,step)
%---------------------------
% We extract the data fields
%---------------------------
  if strcmp(field,'Freezing_Temperature')
     fld = 'SST';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'SSS';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp+0.055*field_tmp2;
  elseif strcmp(field,'Ridged_volume_per_area')
     fld = 'Ridge_ratio';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Thickness';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp.*field_tmp2;
  elseif strcmp(field,'Total_thickness')
     fld = 'Thickness';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Thin_ice';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);
     
     field_tmp = field_tmp+field_tmp2;
     field_plotted='Thick and thin ice thickness';
  elseif strcmp(field,'Total_concentration')
     fld = 'Concentration';
     [field_tmp]=get_and_check(fld,data_out,dirname,step);

     fld = 'Concentration_thin_ice';
     [field_tmp2]=get_and_check(fld,data_out,dirname,step);

     field_tmp = field_tmp+field_tmp2;
     field_plotted='Thick and thin ice concentration';
  else
    [field_tmp]=get_and_check(field,data_out,dirname,step);
    field_plotted=field;
  end
end
