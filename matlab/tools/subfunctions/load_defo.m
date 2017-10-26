function [defo,ebox]=load_defo(defo_file)

check_entry=exist(defo_file)
switch(check_entry)
    case 2
        defo_test=load(defo_file);
        list_field=fieldnames(defo_test);
        defo=getfield(defo_test,list_field{1});
        ebox=getfield(defo_test,list_field{2});
    otherwise
        disp(check_entry)
        error('data type not treated')
end

