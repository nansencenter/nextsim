% Default path for the data
% This path has also to be in your matlab path
indir='/net/sverdrup-1/vol/sim/data/';
forecast_dir='/Home/phigri/forecasting/forecasts/'; % only works when loggued on Johansen

% If you have copy the data locally, you can define another path in a file
% that has to be named 'define_personal_data_path.m'
if(exist('define_personal_data_path.m'))
    define_personal_data_path();
end