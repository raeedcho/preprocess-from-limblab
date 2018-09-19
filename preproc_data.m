%% Set up
meta.lab=6;
meta.ranBy='Raeed';
meta.monkey='Butter';
meta.date='20180522';
meta.task='TRT'; % for the loading of cds
meta.taskAlias={'TRT_001'}; % for the filename (cell array list for files to load and save)
meta.array='RightCuneate'; % for the loading of cds
meta.project='MultiWorkspace'; % for the folder in data-preproc
meta.hyperfolder=fullfile('C:\Users\rhc307\Projects\limblab'); % folder for both data_preproc and data_td
meta.superfolder=fullfile(meta.hyperfolder,'data-preproc',meta.project,meta.monkey); % folder for data dump
meta.folder=fullfile(meta.superfolder,meta.date); % compose subfolder and superfolder

if strcmp(meta.monkey,'Chips')
    meta.mapfile='C:\Users\rhc307\Projects\limblab\data-preproc\Meta\Mapfiles\Chips\left_S1\SN 6251-001455.cmp';
    meta.arrayAlias='area2'; % for the filename
elseif strcmp(meta.monkey,'Han')
    meta.mapfile='C:\Users\rhc307\Projects\limblab\data-preproc\Meta\Mapfiles\Han\left_S1\SN 6251-001459.cmp';
    if startsWith(meta.date,'2017')
        meta.arrayAlias = 'area2EMG'; % for the filename
        altMeta = meta;
        altMeta.array='';
        altMeta.arrayAlias='EMGextra';
        altMeta.neuralPrefix = [altMeta.monkey '_' altMeta.date '_' altMeta.arrayAlias];
    else
        meta.arrayAlias = 'area2'; % for the filename
    end
elseif strcmp(meta.monkey,'Lando')
    meta.mapfile='C:\Users\rhc307\Projects\limblab\data-preproc\Meta\Mapfiles\Lando\left_S1\SN 6251-001701.cmp';
    meta.arrayAlias = 'area2';
    altMeta = meta;
    altMeta.array='RightCuneate';
    altMeta.arrayAlias='cuneate';
    altMeta.neuralPrefix = [altMeta.monkey '_' altMeta.date '_' altMeta.arrayAlias];
    altMeta.mapfile='C:\Users\rhc307\Projects\limblab\data-preproc\Meta\Mapfiles\Lando\right_cuneate\SN 1025-001745.cmp';
elseif strcmp(meta.monkey,'Butter')
    meta.mapfile='C:\Users\rhc307\Projects\limblab\data-preproc\Meta\Mapfiles\Butter\right_CN\SN 6250-001799.cmp';
    meta.arrayAlias='cuneate'; % for the filename
end

meta.neuralPrefix = [meta.monkey '_' meta.date '_' meta.arrayAlias];

%% Move data into subfolder
if ~exist(meta.folder,'dir')
    mkdir(meta.folder)
    movefile(fullfile(meta.superfolder,[meta.monkey '_' meta.date '*']),meta.folder)
end

%% Set up folder structure
if ~exist(fullfile(meta.folder,'preCDS'),'dir')
    mkdir(fullfile(meta.folder,'preCDS'))
    movefile(fullfile(meta.folder,[meta.neuralPrefix '*']),fullfile(meta.folder,'preCDS'))
    if exist('altMeta','var')
        movefile(fullfile(meta.folder,[altMeta.neuralPrefix '*']),fullfile(meta.folder,'preCDS'))
    end
end
if ~exist(fullfile(meta.folder,'preCDS','merging'),'dir')
    mkdir(fullfile(meta.folder,'preCDS','merging'))
    movefile(fullfile(meta.folder,'preCDS',[meta.neuralPrefix '*.nev']),fullfile(meta.folder,'preCDS','merging'))
    if exist('altMeta','var') && ~isempty(altMeta.array)
        movefile(fullfile(meta.folder,'preCDS',[altMeta.neuralPrefix '*.nev']),fullfile(meta.folder,'preCDS','merging'))
    end
end
if ~exist(fullfile(meta.folder,'preCDS','Final'),'dir')
    mkdir(fullfile(meta.folder,'preCDS','Final'))
    movefile(fullfile(meta.folder,'preCDS',[meta.neuralPrefix '*.n*']),fullfile(meta.folder,'preCDS','Final'))
    if exist('altMeta','var')
        movefile(fullfile(meta.folder,'preCDS',[altMeta.neuralPrefix '*.n*']),fullfile(meta.folder,'preCDS','Final'))
    end
end
if ~exist(fullfile(meta.folder,'ColorTracking'),'dir')
    mkdir(fullfile(meta.folder,'ColorTracking'))
    movefile(fullfile(meta.folder,'*_colorTracking_*.mat'),fullfile(meta.folder,'ColorTracking'))
end
if ~exist(fullfile(meta.folder,'ColorTracking','Markers'),'dir')
    mkdir(fullfile(meta.folder,'ColorTracking','Markers'))
end
if ~exist(fullfile(meta.folder,'OpenSim'),'dir')
    mkdir(fullfile(meta.folder,'OpenSim'))
end
if ~exist(fullfile(meta.folder,'OpenSim','Analysis'),'dir')
    mkdir(fullfile(meta.folder,'OpenSim','Analysis'))
end
if ~exist(fullfile(meta.folder,'CDS'),'dir')
    mkdir(fullfile(meta.folder,'CDS'))
end
if ~exist(fullfile(meta.folder,'TD'),'dir')
    mkdir(fullfile(meta.folder,'TD'))
end

%% Merge and strip files for spike sorting
% Run processSpikesForSorting for the first time to combine spike data from
% all files with a name starting with file_prefix.
processSpikesForSorting(fullfile(meta.folder,'preCDS','merging'),meta.neuralPrefix);
if exist('altMeta','var') && ~isempty(altMeta.array)
    processSpikesForSorting(fullfile(altMeta.folder,'preCDS','merging'),altMeta.neuralPrefix);
end

% Now sort in Offline Sorter!

%% Load colorTracking file (and settings if desired) -- NOTE: Can do this simultaneously with sorting, since it takes some time
first_time = 1;
for fileIdx = 1:length(meta.taskAlias)
    colorTrackingFilename = [meta.monkey '_' meta.date '_colorTracking_' meta.taskAlias{fileIdx}];

    fname_load=ls(fullfile(meta.folder,'ColorTracking',[colorTrackingFilename '*']));
    load(deblank(fullfile(meta.folder,'ColorTracking',fname_load)))

    % Run color tracking script
    color_tracker_4colors_script;

    % Save
    markersFilename = [meta.monkey '_' meta.date '_markers_' meta.taskAlias{fileIdx}];
    fname_save=fullfile(meta.folder,'ColorTracking','Markers',[markersFilename '.mat']);
    save(fname_save,'all_medians','all_medians2','led_vals','times');

    if first_time
        fname_save_settings=fullfile(meta.folder,'ColorTracking','Markers',['settings_' meta.monkey '_' meta.date]);
        save(fname_save_settings,'red_elbow_dist_from_blue','red_blue_arm_dist_max',...
            'green_hand_dists_elbow','red_hand_dists_elbow','blue_hand_dists_elbow','yellow_hand_dists_elbow','green_separator',...
            'green_hand_dists_bluearm','red_hand_dists_bluearm','blue_hand_dists_bluearm','yellow_hand_dists_bluearm',...
            'green_hand_dists_redarm', 'red_hand_dists_redarm', 'blue_hand_dists_redarm','yellow_hand_dists_redarm',...
            'green_dist_min','red_keep','green_keep','blue_keep','yellow_keep','marker_inits');
        clearvars -except meta altMeta 'red_elbow_dist_from_blue' 'red_blue_arm_dist_max'...
            'green_hand_dists_elbow' 'red_hand_dists_elbow' 'blue_hand_dists_elbow' 'yellow_hand_dists_elbow' 'green_separator' ...
            'green_hand_dists_bluearm' 'red_hand_dists_bluearm' 'blue_hand_dists_bluearm' 'yellow_hand_dists_bluearm' ...
            'green_hand_dists_redarm'  'red_hand_dists_redarm'  'blue_hand_dists_redarm' 'yellow_hand_dists_redarm' ...
            'green_dist_min' 'red_keep' 'green_keep' 'blue_keep' 'yellow_keep' 'marker_inits'
        first_time = 0;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except meta altMeta

%% Split files and move to Final folder before loading
processSpikesForSorting(fullfile(meta.folder,'preCDS','merging'),meta.neuralPrefix);
if exist('altMeta','var') && ~isempty(altMeta.array)
    processSpikesForSorting(fullfile(altMeta.folder,'preCDS','merging'),altMeta.neuralPrefix);
end

% move into final folder
for fileIdx = 1:length(meta.taskAlias)
    movefile(fullfile(meta.folder,'preCDS','merging',[meta.neuralPrefix '_' meta.taskAlias{fileIdx} '.mat']),...
        fullfile(meta.folder,'preCDS','Final'));
    if exist('altMeta','var') && ~isempty(altMeta.array)
        movefile(fullfile(altMeta.folder,'preCDS','merging',[altMeta.neuralPrefix '_' altMeta.taskAlias{fileIdx} '.mat']),...
            fullfile(altMeta.folder,'preCDS','Final'));
    end
    movefile(fullfile(meta.folder,'preCDS','merging',[meta.neuralPrefix '_' meta.taskAlias{fileIdx} '.nev']),...
        fullfile(meta.folder,'preCDS','Final'));
    if exist('altMeta','var') && ~isempty(altMeta.array)
        movefile(fullfile(altMeta.folder,'preCDS','merging',[altMeta.neuralPrefix '_' altMeta.taskAlias{fileIdx} '.nev']),...
            fullfile(altMeta.folder,'preCDS','Final'));
    end
end

%% Load data into CDS file
% Make CDS files
cds = cell(size(meta.taskAlias));
for fileIdx = 1:length(meta.taskAlias)
    cds{fileIdx} = commonDataStructure();
    cds{fileIdx}.file2cds(fullfile(meta.folder,'preCDS','Final',[meta.neuralPrefix '_' meta.taskAlias{fileIdx}]),...
        ['ranBy' meta.ranBy],['array' meta.array],['monkey' meta.monkey],meta.lab,'ignoreJumps',['task' meta.task],['mapFile' meta.mapfile]);

    % also load second file if necessary
    if exist('altMeta','var')
        if ~isempty(altMeta.array)
            cds{fileIdx}.file2cds(fullfile(altMeta.folder,'preCDS','Final',[altMeta.neuralPrefix '_' altMeta.taskAlias{fileIdx}]),...
                ['ranBy' altMeta.ranBy],['array' altMeta.array],['monkey' altMeta.monkey],altMeta.lab,'ignoreJumps',['task' altMeta.task],['mapFile' altMeta.mapfile]);
        else
            cds{fileIdx}.file2cds(fullfile(altMeta.folder,'preCDS','Final',[altMeta.neuralPrefix '_' altMeta.taskAlias{fileIdx}]),...
                ['ranBy' altMeta.ranBy],['monkey' altMeta.monkey],altMeta.lab,'ignoreJumps',['task' altMeta.task]);
        end
    end
end

%% Load marker file and create TRC
for fileIdx = 1:length(meta.taskAlias)
    markersFilename = [meta.monkey '_' meta.date '_markers_' meta.taskAlias{fileIdx} '.mat'];
    affine_xform = cds{fileIdx}.loadRawMarkerData(fullfile(meta.folder,'ColorTracking','Markers',markersFilename));
    writeTRCfromCDS(cds{fileIdx},fullfile(meta.folder,'OpenSim',[meta.monkey '_' meta.date '_' meta.taskAlias{fileIdx} '_markerData.trc']))
    writeHandleForceFromCDS(cds{fileIdx},fullfile(meta.folder,'OpenSim',[meta.monkey '_' meta.date '_' meta.taskAlias{fileIdx} '_handleForce.mot']))
end

%% Do openSim stuff and save analysis results to analysis folder

% do this in opensim for now

%% Add kinematic information to CDS
for fileIdx = 1:length(meta.taskAlias)
%     % load joint information
%     cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'joint_ang')
% 
%     % load joint velocities
%     cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'joint_vel')
% 
%     % load joint moments
%     % cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'joint_dyn')
% 
%     % load muscle information
%     cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'muscle_len')
% 
%     % load muscle velocities
%     cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'muscle_vel')
%     
%     % load hand positions
%     cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'hand_pos')
%     
%     % load hand velocities
%     cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'hand_vel')
%     
%     % load hand accelerations
%     cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'hand_acc')
% 
%     % load elbow positions
%     cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'elbow_pos')
%     
%     load elbow velocities
%     cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'elbow_vel')
%     
%     load elbow accelerations
%     cds{fileIdx}.loadOpenSimData(fullfile(meta.folder,'OpenSim','Analysis'),'elbow_acc')
end

%% Save CDS
save(fullfile(meta.folder,'CDS',[meta.monkey '_' meta.date '_CDS.mat']),'cds','-v7.3')

%% Make TD
% COactpas
% params.array_alias = {'LeftS1Area2','S1'};
% params.exclude_units = [255];
% params.event_list = {'ctrHoldBump';'bumpTime';'bumpDir';'ctrHold'};
% params.trial_results = {'R','A','F','I'};
% td_meta = struct('task',meta.task);
% params.meta = td_meta;
% trial_data = parseFileByTrial(cds{1},params);
% td_taskname = 'COactpas';

% OOR
% params.array_alias = {'LeftS1Area2','S1'};
% % params.exclude_units = [255];
% params.event_list = {'tgtDir','target_direction';'forceDir','force_direction';'startTargHold','startTargHoldTime';'endTargHoldTime','endTargHoldTime'};
% params.trial_results = {'R','A','F','I'};
% td_meta = struct('task','OOR');
% params.meta = td_meta;
% 
% trial_data = parseFileByTrial(cds{1},params);

% Bumpcurl
% params.array_alias = {'LeftS1Area2','S1'};
% params.event_list = {'ctrHoldBump';'bumpTime';'bumpDir';'ctrHold'};
% td_meta = struct('task',meta.task,'epoch','BL');
% params.trial_results = {'R','A','F','I'};
% params.meta = td_meta;
% trial_data_BL = parseFileByTrial(cds{1},params);
% params.meta.epoch = 'AD';
% trial_data_AD = parseFileByTrial(cds{2},params);
% params.meta.epoch = 'WO';
% trial_data_WO = parseFileByTrial(cds{3},params);
% 
% trial_data = cat(2,trial_data_BL,trial_data_AD,trial_data_WO);
% td_taskname = 'CObumpcurl';

% TRT
params.array_alias = {'LeftS1Area2','S1';'RightCuneate','CN'};
params.event_list = {'bumpTime';'bumpDir';'ctHoldTime';'otHoldTime';'spaceNum';'targetStartTime'};
params.trial_results = {'R','A','F','I'};
td_meta = struct('task',meta.task);
params.meta = td_meta;
trial_data = parseFileByTrial(cds{1},params);
% sanitize?
idxkeep = cat(1,trial_data.spaceNum) == 1 | cat(1,trial_data.spaceNum) == 2;
trial_data = trial_data(idxkeep);
td_taskname = 'TRT';

% RW DL/PM
% params.array_alias = {'LeftS1Area2','S1';'RightCuneate','CN'};
% params.trial_results = {'R','A','F','I'};
% td_meta = struct('task',meta.task,'spaceNum',2);
% params.meta = td_meta;
% trial_data_DL = parseFileByTrial(cds{1},params);
% td_meta = struct('task',meta.task,'spaceNum',1);
% params.meta = td_meta;
% trial_data_PM = parseFileByTrial(cds{2},params);
% trial_data = [trial_data_PM trial_data_DL];
% % match up with TRT
% for trial = 1:length(trial_data)
%     trial_data(trial).idx_targetStartTime = trial_data(trial).idx_startTime;
% end
% trial_data = reorderTDfields(trial_data);
% td_taskname = 'RWTW';

% RW
% params.array_alias = {'LeftS1Area2','S1'};
% params.trial_results = {'R','A','F','I'};
% params.include_ts = true;
% td_meta = struct('task','RW');
% params.meta = td_meta;
% trial_data = parseFileByTrial(cds{1},params);
% % match up with TRT
% for trial = 1:length(trial_data)
%     trial_data(trial).idx_targetStartTime = trial_data(trial).idx_goCueTime(1);
% end
% trial_data = reorderTDfields(trial_data);
% td_taskname = 'RW';

%% Save TD
save(fullfile(meta.folder,'TD',[meta.monkey '_' meta.date '_' td_taskname '_TD.mat']),'trial_data')
copyfile(fullfile(meta.folder,'TD',[meta.monkey '_' meta.date '_' td_taskname '_TD.mat']),fullfile(meta.hyperfolder,'data-td',meta.project))