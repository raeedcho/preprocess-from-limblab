%% Pre-processing parameters
clear meta
meta.ranBy='Raeed';
meta.monkey='Han';
meta.date='20160315';
meta.task={'RW'}; % for the loading of cds
meta.taskAlias={'RW_001'}; % for the filename (cell array list for files to load and save)
meta.td_taskname = 'RW'; % for saving the TD
meta.EMGrecorded = false; % whether or not EMG was recorded
meta.motorcontrol = false; % whether or not motor control signals were recorded
meta.markered = true; % whether or not the colorTracking has been markered

%% Set up meta fields
meta.localdatafolder=fullfile('C:\Users\rhc307\data\'); % folder with data-td and working data folder
meta.cdslibrary=fullfile(meta.localdatafolder,'cds-library');
meta.tdlibrary=fullfile(meta.localdatafolder,'td-library');
meta.remotefolder=fullfile('Z:\limblab\User_folders\Raeed');

%% Make TD
spike_routine = @processCDSspikes;
cds_routine = @processCDS;

cont_signal_names = {...
    'pos',...
    'vel',...
    'acc',...
    'force',...
    };

cont_signal_labels = {...
    {'x','y'},...
    {'vx','vy'},...
    {'ax','ay'},...
    {'fx','fy','fz','mx','my','mz'},...
    };

if meta.motorcontrol
    cont_signal_names = [...
        cont_signal_names,...
        {'motor_control'},...
        ];
    
    cont_signal_labels = [...
        cont_signal_labels,...
        {{'MotorControlSho','MotorControlElb'}},...
        ];
end

if meta.markered
    cont_signal_names = [...
        cont_signal_names,...
        {...
            'markers',...
            'joint_ang',...
            'joint_vel',...
            'muscle_len',...
            'muscle_vel',...
            'opensim_hand_pos',...
            'opensim_hand_vel',...
            'opensim_hand_acc',...
            'opensim_elbow_pos',...
            'opensim_elbow_vel',...
            'opensim_elbow_acc',...
            }...
        ];

    cont_signal_labels = [...
        cont_signal_labels,...
        {...
            sort(getMarkerNames()),...
            strcat(getJointNames(),'_ang'),...
            strcat(getJointNames(),'_vel'),...
            strcat(getMuscleNames(),'_len'),...
            strcat(getMuscleNames(),'_muscVel'),...
            strcat({'X','Y','Z'},{'_handPos'}),...
            strcat({'X','Y','Z'},{'_handVel'}),...
            strcat({'X','Y','Z'},{'_handAcc'}),...
            strcat({'X','Y','Z'},{'_elbowPos'}),...
            strcat({'X','Y','Z'},{'_elbowVel'}),...
            strcat({'X','Y','Z'},{'_elbowAcc'}),...
            }...
        ];
end

if meta.EMGrecorded
    emg_signal_names = getEMGNames();
else
    emg_signal_names = {};
end

% parameters...
clear params
params = struct('bin_size',0.001);

%% load it in
td_cell = cell(1,length(meta.taskAlias));
for fileIdx = 1:length(meta.taskAlias)
    % get event names specific to task
    switch meta.task{fileIdx}
        case 'COactpas'
            event_names = {...
                'startTime',...
                'endTime',...
                'tgtOnTime',...
                'goCueTime',...
                'bumpTime',...
                };
            trial_meta = {...
                'ctrHold',...
                'tgtDir',...
                'ctrHoldBump',...
                'bumpDir',...
                };
        case 'OOR'
            event_names = {...
                'startTime',...
                'endTime',...
                'startTargOnTime',...
                'startTargHoldTime',...
                'goCueTime',...
                'endTargHoldTime',...
                };
            trial_meta = {...
                'tgtDir',...
                'forceDir',...
                };
            params.meta = struct('forceMag',2.5);
        case 'TRT'
            event_names = {...
                'startTime',...
                'endTime',...
                'goCueTime',...
                'bumpTime',...
                'targetStartTime',...
                'ctHoldTime',...
                'otHoldTime',...
                };
            trial_meta = {...
                'spaceNum',...
                'bumpDir',...
                };
        case 'RW'
            event_names = {...
                'startTime',...
                'endTime',...
                'goCueTime',...
                };
            trial_meta = {...
                'tgtSize',...
                'tgtCtr',...
                };
            if strcmpi(meta.td_taskname,'RWTW')
                params.meta = struct('spaceNum',meta.epochname{fileIdx});
            end
        otherwise
            error('not implemented yet')
    end
    
    % construct cds filename
    cds_filename = fullfile(meta.cdslibrary,sprintf('%s_%s_CDS_%s.mat',meta.monkey,meta.date,meta.taskAlias{fileIdx}));
    
    % get signal info
    signal_info = { ...
        initSignalStruct( ...
            'filename',cds_filename, ...
            'routine',spike_routine, ...
            'params',struct(), ... % can pass arbitrary parameters. Shouldn't need much with CDS
            'name','S1', ... % it gets stored under this name... in case of spikes, this gives S1_spikes
            'type','spikes', ... % which type... see documentation of initSignalStruct
            'label',''), ... % label can be indices or strings
        initSignalStruct( ... % continuous data
            'filename',cds_filename, ...
            'routine',cds_routine, ...
            'params',struct('trial_meta',{trial_meta}), ...
            'name',[...
                cont_signal_names,...
                emg_signal_names,...
                event_names,...
                ], ... % stored in this name, matched to the labels below which correspond to the output of the processing routine
            'type',[...
                repmat({'generic'},1,length(cont_signal_names)),...
                repmat({'emg'},1,length(emg_signal_names)),...
                repmat({'event'},1,length(event_names)),...
                ],...
            'label',[...
                cont_signal_labels,...
                strcat('EMG_',emg_signal_names),...
                event_names,...
                ], ... % can also pass [1 2],[3 4] etc if you know the arrangment of the signals in the data struct
            'operation',[]), ...
        };
    
    % load trial_data (will result in warning for no meta info, but we're
    % taking meta info from the CDS anyway)
    td_cell{fileIdx} = convertDataToTD(signal_info,params);
    
    % add some meta information
    if meta.markered
        % add label names
        td_cell{fileIdx}.marker_names = sort(getMarkerNames());
        td_cell{fileIdx}.joint_names = getJointNames();
        td_cell{fileIdx}.muscle_names = getMuscleNames();
    end
    
    if meta.motorcontrol
        td_cell{fileIdx}.motorcontrol_names = {'MotorControlSho','MotorControlElb'};
    end
    
    % make it pretty
    td_cell{fileIdx} = reorderTDfields(td_cell{fileIdx});
end

%% save structure
% put trial datas together
trial_data = horzcat(td_cell{:});
save(fullfile(meta.tdlibrary,[meta.monkey '_' meta.date '_' meta.td_taskname '_TD.mat']),'trial_data','-v7.3')

% and back it up
winopen(fullfile(meta.remotefolder,'td-library'))
winopen(meta.tdlibrary)

%% Post processing?
% switch(meta.task)  
%     case 'OOR'
%         trial_data = td_cell{1};
%         params.array_alias = {'LeftS1Area2','S1'};
%         % params.exclude_units = [255];
%         params.event_list = {'tgtDir','target_direction';'forceDir','force_direction';'startTargHold','startTargHoldTime';'endTargHoldTime','endTargHoldTime'};
%         params.trial_results = {'R','A','F','I'};
%         td_meta = struct('task','OOR');
%         params.meta = td_meta;
%         
%         trial_data = parseFileByTrial(cds_cell{1},params);
% 
%     case 'CObumpcurl'
%         params.array_alias = {'LeftS1Area2','S1'};
%         params.event_list = {'ctrHoldBump';'bumpTime';'bumpDir';'ctrHold'};
%         td_meta = struct('task',meta.task,'epoch','BL');
%         params.trial_results = {'R','A','F','I'};
%         params.meta = td_meta;
%         trial_data_BL = parseFileByTrial(cds_cell{1},params);
%         params.meta.epoch = 'AD';
%         trial_data_AD = parseFileByTrial(cds_cell{2},params);
%         params.meta.epoch = 'WO';
%         trial_data_WO = parseFileByTrial(cds_cell{3},params);
%         
%         trial_data = cat(2,trial_data_BL,trial_data_AD,trial_data_WO);

%     case 'RW'
%         if any(contains(meta.taskAlias,'DL'))
%             fprintf('Interpreting as RWTW task')
%             params.array_alias = {'LeftS1Area2','S1';'RightCuneate','CN'};
%             params.trial_results = {'R','A','F','I'};
%             td_meta = struct('task',meta.task,'spaceNum',2);
%             params.meta = td_meta;
%             trial_data_DL = parseFileByTrial(cds_cell{contains(meta.taskAlias,'DL')},params);
%             td_meta = struct('task',meta.task,'spaceNum',1);
%             params.meta = td_meta;
%             trial_data_PM = parseFileByTrial(cds_cell{contains(meta.taskAlias,'PM')},params);
%             trial_data = [trial_data_PM trial_data_DL];
%             % match up with TRT
%             for trial = 1:length(trial_data)
%                 trial_data(trial).idx_targetStartTime = trial_data(trial).idx_startTime;
%             end
%             trial_data = reorderTDfields(trial_data);
%         else
%             trial_data = td_cell{1};
%             fprintf('Interpreting as regular random walk')
%             params.array_alias = {'LeftS1Area2','S1'};
%             params.trial_results = {'R','A','F','I'};
%             params.include_ts = true;
%             td_meta = struct('task','RW');
%             params.meta = td_meta;
%             trial_data = parseFileByTrial(cds_cell{1},params);
%             % match up with TRT
%             for trial = 1:length(trial_data)
%                 trial_data(trial).idx_targetStartTime = trial_data(trial).idx_goCueTime(1);
%             end
%             trial_data = reorderTDfields(trial_data);
%         end
%         
%     otherwise
%         error('Unrecognized task')
%         
% end
