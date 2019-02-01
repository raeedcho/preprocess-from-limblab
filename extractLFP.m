%% extract lfp
lfp = cds.lfp(:,2:end);

% first add nans to start
new_t = (0:1/2000:cds.lfp.t(end))';
lfp_interp = interp1(cds.lfp.t,lfp.Variables,new_t);

%% now change elecs to chans to match units
lfp_labels = lfp.Properties.VariableNames;
chans = [cds.units.chan]';
elecs = sscanf([cds.units.label],'elec%d');
unit_guide = [elecs chans];
unit_guide = double(unique(unit_guide,'rows'));

label_elec = sscanf(horzcat(lfp_labels{:}),'LeftS1Area2elec%d');
label_chan = interp1(unit_guide(:,1),unit_guide(:,2),label_elec);
new_labels = strsplit(sprintf('chan%d ',label_chan),' ');
new_labels(end) = [];
lfp_table = array2table([new_t lfp_interp],'VariableNames',[{'t'} new_labels]);
