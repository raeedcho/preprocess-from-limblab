function out = processCDSevents(filename,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loads a CDS file and returns a field

% parameters
event_names = {};
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

if isempty(event_names)
    % just pick out everything that has 'Time' in the name
    event_names = cds.trials.Properties.VariableNames(endsWith(cds.trials.Properties.VariableNames,'Time'));
end

% if events are time, it expects them in cells like spiking data
%   you can also give it already-binned events and it just passes them
%   along
data = cell(1,length(event_names));
for eventnum = 1:length(event_names)
    data{eventnum} = cds.trials.(event_names{eventnum});
end
t = (cds.meta.dataWindow(1):1/1000:cds.meta.dataWindow(2))';

% this part could be automated
labels = event_names;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.meta   = [];
out.data   = data;
out.labels = labels;
out.t      = t;
out.error_flag = error_flag;

