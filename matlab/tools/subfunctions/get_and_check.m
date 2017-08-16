function [field_tmp]=get_and_check(fld,data_out,dirname,step,no_error)

    if ~exist('no_error','var'); no_error = 0; end

    try
        field_tmp=data_out.(fld);
     catch ME
        disp([fld,' not present in ',dirname,'(step=',num2str(step),')']);
        disp(' ');
        disp('Available fields are:');
        disp(' ');

        flds   = fieldnames(data_out);
        for k=1:length(flds)
           disp(flds{k});
        end

        if no_error
           field_tmp = [];
        else
           error([fld,' not present in ',dirname,'(step=',num2str(step),')']);
        end
     end
end
