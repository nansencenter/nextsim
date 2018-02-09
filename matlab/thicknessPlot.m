function thicknessPlot(files)

% function thicknessPlot(files)
% 
% files:    A file name pattern for the files we wish to plot see the Matlab function 'dir'. Example: 'longRun*.mat'

load CS2_SMOS.mat

fh = figure;

% Fit to display width
screenSize = get(0,'ScreenSize');
fh.Position = [1 925 screenSize(3) 420];

hold on

for year=2010:2015
    
    % Don't plot non-observed periods, and don't plot the first bit we know is wrong
    tmin = datenum(year,   09, 01);
    %tmin =     datenum(year,   11, 10);
    tmax = min(datenum(year+1, 05, 01), datenum(2015,12,31));
    
    tmask = tCS2_SMOS>=tmin & tCS2_SMOS<=tmax;
    
    shadedErrorBar(tCS2_SMOS(tmask), hCS2_SMOS_ICE(tmask), ...
        [ hmaxCS2_SMOS_ICE(tmask)-hCS2_SMOS_ICE(tmask); ...
          hCS2_SMOS_ICE(tmask)-hminCS2_SMOS_ICE(tmask) ] );
      
end

% Loop over input files
matFiles = dir(files);
legends = cell(1,length(matFiles));
ph = zeros(1,length(matFiles));

for i=1:length(matFiles)
    
    load(matFiles(i).name)
    if iscell(t) % Legacy support for Einar's old scripts
        ph(i) = plot(t{1}, hice_ICE{1});
    else
        ph(i) = plot(t, hice_ICE);
    end
    legends{i} = matFiles(i).name(1:end-4);
end

legend(ph, legends, 'location', 'northwest', 'Interpreter','none')

xlim([datenum(2010,11,15) datenum(2016,01,01)])
datetick('x','mmm ''yy')
grid
ylabel('Mean thickness [m]')

% Replace X with Time
hdt = datacursormode(fh);
set(hdt,'UpdateFcn', {@dtips})

end

function output_txt = dtips(~,event_obj)

% Replace X with Time

pos = get(event_obj,'Position');
% h = get(event_obj,'Target');

X = pos(1);
Y = pos(2);

output_txt = { ...
    ['Value: ', num2str(Y)] ...
    ['Time:  ', datestr(X)]};

end

