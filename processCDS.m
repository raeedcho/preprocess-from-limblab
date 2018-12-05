function out = processCDS(filename,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loads a CDS file and returns all continuous signals it can extract, along
% with all the event data (trial table entries ending in 'Time').

% params
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first the events...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% this part could be automated
event_labels = event_names;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now the continuous signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure out which signals to extract
signal_names = {};
if cds.meta.hasKinematics
    signal_names = [signal_names {'kin'}];
end
if cds.meta.hasForce
    signal_names = [signal_names {'force'}];
end
if cds.meta.hasEmg
    signal_names = [signal_names {'emg'}];
end
if cds.meta.hasLfp
    % signal_names = [signal_names {'lfp'}];
end
if cds.meta.hasAnalog
    % analog could be a lot of stuff...
    for analog_idx = 1:length(cds.analog)
        header = cds.analog{analog_idx}.Properties.VariableNames;
        % Figure out if we have motor control data
        if any(contains(header,'MotorControl'))
            signal_names = [signal_names {'motorcontrol'}];
            motorcontrol_idx = analog_idx;
            continue
        end
        if any(contains(header,'Frame')) || any(contains(header,'Marker'))
            signal_names = [signal_names {'markers'}];
            markers_idx = analog_idx;
            continue
        end
        if any(endsWith(header,'_ang'))
            signal_names = [signal_names {'joint_ang'}];
            joint_ang_idx = analog_idx;
        end
        if any(endsWith(header,'_vel'))
            signal_names = [signal_names {'joint_vel'}];
            joint_vel_idx = analog_idx;
        end
        if any(endsWith(header,'_len'))
            signal_names = [signal_names {'muscle_len'}];
            muscle_len_idx = analog_idx;
        end
        if any(endsWith(header,'_muscVel'))
            signal_names = [signal_names {'muscle_vel'}];
            muscle_vel_idx = analog_idx;
        end
        if any(contains(header,'_handPos'))
            signal_names = [signal_names {'hand_pos'}];
            hand_pos_idx = analog_idx;
        end
        if any(contains(header,'_handVel'))
            signal_names = [signal_names {'hand_vel'}];
            hand_vel_idx = analog_idx;
        end
        if any(contains(header,'_handAcc'))
            signal_names = [signal_names {'hand_acc'}];
            hand_acc_idx = analog_idx;
        end
        if any(contains(header,'_elbowPos'))
            signal_names = [signal_names {'elbow_pos'}];
            elbow_pos_idx = analog_idx;
        end
        if any(contains(header,'_elbowVel'))
            signal_names = [signal_names {'elbow_vel'}];
            elbow_vel_idx = analog_idx;
        end
        if any(contains(header,'_elbowAcc'))
            signal_names = [signal_names {'elbow_acc'}];
            elbow_acc_idx = analog_idx;
        end
    end
end

% extract signals
samp_rate = zeros(1,length(signal_names));
[timevec_cell,cont_data_cell,signal_labels] = deal(cell(1,length(signal_names)));
sig_guides = struct();
for signum = 1:length(signal_names)
    switch lower(signal_names{signum})
        case 'kin'
            % do kin stuff
            data_table = cds.kin;
            data_cols = startsWith(data_table.Properties.VariableNames,{'x','y','vx','vy','ax','ay'});
        case 'force'
            % do force stuff
            data_table = cds.force;
            data_cols = startsWith(data_table.Properties.VariableNames,{'fx','fy','fz','mx','my','mz'});
        case 'emg'
            % do emg stuff
            error('emg is not implemented')
        case 'motorcontrol'
            % do motor control stuff
            data_table = cds.analog{motorcontrol_idx};
            data_cols = startsWith(data_table.Properties.VariableNames,'MotorControl');
            sig_guides.motorcontrol_names = data_table.Properties.VariableNames(data_cols);
        case 'markers'
            % do marker stuff
            marker_table = cds.analog{marker_idx};
            marker_cols = ~strcmpi(marker_table.Properties.VariableNames,'Frame') & ~strcmpi(marker_table.Properties.VariableNames,'t');
            marker_names = cell(1,width(marker_table(:,marker_cols)));
            for marker_idx = 1:length(marker_cols)
                marker_names(((marker_idx-1)*3+1):(marker_idx*3)) = strcat(marker_cols(marker_idx),{'_y','_z','_x'});
            end
            
            marker_data = marker_table{:,marker_cols};
            
            % set in data_table
            data_table = table([t marker_data
            error('markers are not yet implemented')
        case 'joint_ang'
            % do opensim stuff
            data_table = cds.analog{joint_ang_idx};
            data_cols = endsWith(data_table.Properties.VariableNames,'_ang');
            sig_guides.joint_names = strrep(data_table.Properties.VariableNames(data_cols),'_ang','');
        case 'joint_vel'
            % do opensim stuff
            data_table = cds.analog{joint_vel_idx};
            data_cols = endsWith(data_table.Properties.VariableNames,'_vel');
            sig_guides.joint_names = strrep(data_table.Properties.VariableNames(data_cols),'_vel','');
        case 'muscle_len'
            % do opensim stuff
            data_table = cds.analog{muscle_len_idx};
            data_cols = endsWith(data_table.Properties.VariableNames,'_len');
            sig_guides.muscle_names = strrep(data_table.Properties.VariableNames(data_cols),'_len','');
        case 'muscle_vel'
            % do opensim stuff
            data_table = cds.analog{muscle_vel_idx};
            data_cols = endsWith(data_table.Properties.VariableNames,'_muscVel');
            sig_guides.muscle_names = strrep(data_table.Properties.VariableNames(data_cols),'_muscVel','');
        case 'hand_pos'
            % do opensim stuff
            data_table = cds.analog{hand_pos_idx};
            data_cols = find(contains(data_table.Properties.VariableNames,'_handPos'));
            % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
        case 'hand_vel'
            % do opensim stuff
            data_table = cds.analog{hand_vel_idx};
            data_cols = find(contains(data_table.Properties.VariableNames,'_handVel'));
            % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
        case 'hand_acc'
            % do opensim stuff
            data_table = cds.analog{hand_acc_idx};
            data_cols = find(contains(data_table.Properties.VariableNames,'_handAcc'));
            % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
        case 'elbow_pos'
            % do opensim stuff
            data_table = cds.analog{elbow_pos_idx};
            data_cols = find(contains(data_table.Properties.VariableNames,'_elbowPos'));
            % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
        case 'elbow_vel'
            % do opensim stuff
            data_table = cds.analog{elbow_vel_idx};
            data_cols = find(contains(data_table.Properties.VariableNames,'_elbowVel'));
            % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
        case 'elbow_acc'
            % do opensim stuff
            data_table = cds.analog{elbow_acc_idx};
            data_cols = find(contains(data_table.Properties.VariableNames,'_elbowAcc'));
            % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
        otherwise
            error('No idea what this signal is (%s)',signal_names{signum})
    end
    signal_labels{signum} = data_table.Properties.VariableNames(data_cols);
    cont_data_cell{signum} = data_table{:,data_cols};
    timevec_cell{signum} = data_table.t;
    samp_rate(signum) = 1/mode(diff(data_table.t));
end

% resample everything to the highest sampling rate
[final_rate,maxrate_idx] = max(samp_rate);
t_start = -inf;
t_end = inf;
for signum = 1:length(signal_names)
    [P,Q] = rat(final_rate/samp_rate(signum),1e-7);
    if P~=1 || Q~=1
        assert(signum~=maxrate_idx,'Something went wrong with the resample code...')
        
        % need to resample
        cont_data_cell{signum} = resample(cont_data_cell{signum},P,Q);
        
        % resample time vector
        assert(~any(timevec_cell{signum}<=0),'Why are there negative times in this CDS?')
        timevec = downsample(upsample(timevec_cell{signum},P),Q);
        timevec_cell{signum} = interp1(find(timevec>0),timevec(timevec>0),(1:length(timevec))');
        
        % get rid of NaNs on end from upsampling
        nanners = isnan(timevec_cell{signum});
        timevec_cell{signum}(nanners) = [];
        cont_data_cell{signum}(nanners,:) = [];
    end
    % collect largest min and smallest max time for trimming
    t_start = max(t_start,timevec_cell{signum}(1));
    t_end = min(t_end,timevec_cell{signum}(end));
end

% trim time vectors to be the same length, if they're not already
for signum = 1:length(signal_names)
    trim_pts = timevec_cell{signum}<t_start | timevec_cell{signum}>t_end;
    cont_data_cell{signum}(trim_pts,:) = [];
    timevec_cell{signum}(trim_pts) = [];
end

% try horizontally concatenating...If everything went well, things should
% be the right length...
cont_data = horzcat(cont_data_cell{:});
cont_labels = horzcat(signal_labels{:});
t = timevec_cell{maxrate_idx};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.meta   = sig_guides;
out.cont_data   = cont_data;
out.cont_labels = cont_labels;
out.event_data   = event_data;
out.event_labels = event_labels;
out.t      = t;
out.error_flag = error_flag;

