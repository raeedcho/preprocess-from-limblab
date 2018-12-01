function out = processCDSspikes(filename,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loads a CDS file and returns a field
spiking_chans  = 1:96;
exclude_units  = 255; % sort id of units to exclude
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_flag = false;

% load the CDS
if ~isempty(filename)
    load(filename);
else
    error_flag = true;
    disp(['ERROR: ' mfilename ': no filename provided']);
end

[data, wf] = deal(cell(1,length(cds.units)));
for unit = 1:length(cds.units)
    data{unit} = cds.units(unit).spikes.ts;
    wf{unit} = cds.units(unit).spikes.wave; % waveforms aren't supported right now in convertDataToTD
end

% assume right now that the blackrock sampling is 30kHz
t = (cds.meta.dataWindow(1):1/double(30000):cds.meta.dataWindow(2))';

labels = [vertcat(cds.units.chan) vertcat(cds.units.ID)];

% remove bogus units
bad_idx = ~ismember(labels(:,1),spiking_chans) | ismember(labels(:,2),exclude_units);
labels = labels(~bad_idx,:);
data = data(~bad_idx);
wf = wf(~bad_idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.meta   = [];
out.data   = data;
out.wf     = wf;
out.labels = labels;
out.t      = t;
out.error_flag = error_flag;

