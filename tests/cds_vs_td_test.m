%% Testing signals + events
figure
ax(1) = subplot(2,1,1);
plot(cds_cell{1}.analog{3}.t,cds_cell{1}.analog{3}.Marker_1)
hold on
plot(repmat(cds_cell{1}.trials.targetStartTime',2,1),repmat(ylim',1,height(cds_cell{1}.trials)),'--k','linewidth',2)
plot(repmat(cds_cell{1}.trials.ctHoldTime',2,1),repmat(ylim',1,height(cds_cell{1}.trials)),'--r','linewidth',2)
ax(2) = subplot(2,1,2);
plot((1:length(trial_data.pos))*trial_data.bin_size,trial_data.markers(:,1:3))
hold on
plot(trial_data.bin_size*repmat(trial_data.idx_targetStartTime,2,1),repmat(ylim',1,length(trial_data.idx_targetStartTime)),'--k','linewidth',2)
plot(trial_data.bin_size*repmat(trial_data.idx_ctHoldTime,2,1),repmat(ylim',1,length(trial_data.idx_ctHoldTime)),'--r','linewidth',2)

linkaxes(ax,'xy')

%% Testing processCDS...
figure
ax(1) = subplot(2,1,1);
plot(cds_cell{1}.kin.t,cds_cell{1}.kin.x)
ax(2) = subplot(2,1,2);
plot(out.t,out.cont_data(:,1))

linkaxes(ax,'xy')

%% Testing spikes (line is faster than scatter)
figure
props = {'LineStyle','none','Marker','.','MarkerEdge','k','MarkerSize',6};

cds_units = cds.units;
cds_units = cds_units(vertcat(cds_units.ID)~=255); % exclude invalidated
ax(1) = subplot(2,1,1);
for unit_idx = 1:10
    ts = cds_units(unit_idx).spikes.ts;
    ts = ts(1:min(1000,length(ts)));
    line([ts ts],[unit_idx unit_idx],props{:})
    hold on
end
title 'CDS units'


td_time = trial_data.bin_size*(1:length(trial_data.pos))';
ax(2) = subplot(2,1,2);
for unit_idx = 1:10
    ts = td_time(trial_data.S1_spikes(:,unit_idx)>0);
    ts = ts(1:min(1000,length(ts)));
    line([ts ts],[unit_idx unit_idx],props{:})
    hold on
end
title 'TD units'

linkaxes(ax,'xy')