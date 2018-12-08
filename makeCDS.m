%% Pre-processing parameters
clear meta altMeta
meta.lab=6;
meta.ranBy='Raeed';
meta.monkey='Han';
meta.date='20171116';
meta.task={'COactpas'}; % for the loading of cds
meta.taskAlias={'COactpas_002'}; % for the filename (cell array list for files to load and save)
meta.EMGrecorded = false; % whether or not EMG was recorded
meta.motionTracked = false; % whether or not we have motion tracking
meta.sorted = true; % whether or not the neurons have already been sorted
meta.markered = false; % whether or not the colorTracking has been markered
meta.array='LeftS1Area2'; % for the loading of cds

%% Set up meta fields
meta.localdatafolder=fullfile('C:\Users\rhc307\data\'); % folder with data-td and working data folder
meta.workingfolder=fullfile(meta.localdatafolder,'workspace'); % folder
meta.cdslibrary=fullfile(meta.localdatafolder,'cds-library');
meta.tdlibrary=fullfile(meta.localdatafolder,'td-library');
meta.remotefolder=fullfile('Z:\limblab\User_folders\Raeed');
meta.semirawfolder=fullfile(meta.remotefolder,'semi-raw');
meta.markersfolder=fullfile(meta.semirawfolder,'markers');
meta.sortedfolder=fullfile(meta.semirawfolder,'sorted');
meta.opensimsettingsfolder=fullfile(meta.remotefolder,'opensim-settings');

if strcmp(meta.monkey,'Chips')
    meta.rawfolder=fullfile('Z:\data\Chips_12H1\RAW');
    meta.mapfile=fullfile(meta.localdatafolder,'limblab-data\metadata\mapfiles\Chips\left_S1\SN 6251-001455.cmp');
    meta.arrayAlias='area2'; % for the filename
elseif strcmp(meta.monkey,'Han')
    meta.rawfolder=fullfile('Z:\data\Han_13B1\Raw');
    meta.mapfile=fullfile(meta.localdatafolder,'limblab-data\metadata\mapfiles\Han\left_S1\SN 6251-001459.cmp');
    if meta.EMGrecorded
        meta.arrayAlias = 'area2EMG'; % for the filename
        altMeta = meta;
        altMeta.array='';
        altMeta.arrayAlias='EMGextra';
        altMeta.neuralPrefix = [altMeta.monkey '_' altMeta.date '_' altMeta.arrayAlias];
    else
        meta.arrayAlias = 'area2'; % for the filename
    end
elseif strcmp(meta.monkey,'Lando')
    meta.rawfolder=fullfile('Z:\data\Lando_13B2\Raw');
    meta.mapfile=fullfile(meta.localdatafolder,'limblab-data\metadata\mapfiles\Lando\left_S1\SN 6251-001701.cmp');
    meta.arrayAlias = 'area2';
    altMeta = meta;
    altMeta.array='RightCuneate';
    altMeta.arrayAlias='cuneate';
    altMeta.neuralPrefix = [altMeta.monkey '_' altMeta.date '_' altMeta.arrayAlias];
    altMeta.mapfile=fullfile(meta.localdatafolder,'limblab-data\metadata\mapfiles\Lando\right_cuneate\SN 1025-001745.cmp');
elseif strcmp(meta.monkey,'Butter')
    meta.mapfile=fullfile(meta.hyperfolder,'limblab-data\metadata\mapfiles\Butter\right_CN\SN 6250-001799.cmp');
    meta.arrayAlias='cuneate'; % for the filename
end

meta.neuralPrefix = [meta.monkey '_' meta.date '_' meta.arrayAlias];

%% Copy data into working directory
if length(dir(meta.workingfolder))==2 % check if directory is empty first
    if ispc
        % open windows so that we can move data around
        winopen(meta.workingfolder)
        winopen(meta.rawfolder)
    else
        fprintf('These directories are probably all wrong...')
        %fprintf('Please navigate to %s and move data to %s',meta.rawfolder,meta.workingfolder)
        % copyfile(fullfile(meta.rawfolder,[meta.monkey '_' meta.date '*']),meta.workingfolder)
    end
else
    winopen(meta.workingfolder)
    error('The working directory is not empty!')
end

%% Set up folder structure
if ~exist(fullfile(meta.workingfolder,'preCDS'),'dir')
    mkdir(fullfile(meta.workingfolder,'preCDS'))
    movefile(fullfile(meta.workingfolder,[meta.neuralPrefix '*']),fullfile(meta.workingfolder,'preCDS'))
    if exist('altMeta','var')
        movefile(fullfile(meta.workingfolder,[altMeta.neuralPrefix '*']),fullfile(meta.workingfolder,'preCDS'))
    end
end
if ~meta.sorted
    if ~exist(fullfile(meta.workingfolder,'preCDS','merging'),'dir')
        mkdir(fullfile(meta.workingfolder,'preCDS','merging'))
        movefile(fullfile(meta.workingfolder,'preCDS',[meta.neuralPrefix '*.nev']),fullfile(meta.workingfolder,'preCDS','merging'))
        if exist('altMeta','var') && ~isempty(altMeta.array)
            movefile(fullfile(meta.workingfolder,'preCDS',[altMeta.neuralPrefix '*.nev']),fullfile(meta.workingfolder,'preCDS','merging'))
        end
    end
end
if ~exist(fullfile(meta.workingfolder,'preCDS','Final'),'dir')
    mkdir(fullfile(meta.workingfolder,'preCDS','Final'))
    movefile(fullfile(meta.workingfolder,'preCDS',[meta.neuralPrefix '*.n*']),fullfile(meta.workingfolder,'preCDS','Final'))
    if exist('altMeta','var')
        movefile(fullfile(meta.workingfolder,'preCDS',[altMeta.neuralPrefix '*.n*']),fullfile(meta.workingfolder,'preCDS','Final'))
    end
end
if meta.motionTracked
    if ~exist(fullfile(meta.workingfolder,'ColorTracking'),'dir')
        mkdir(fullfile(meta.workingfolder,'ColorTracking'))
        movefile(fullfile(meta.workingfolder,'*_colorTracking_*.mat'),fullfile(meta.workingfolder,'ColorTracking'))
    end
    if ~exist(fullfile(meta.workingfolder,'ColorTracking','Markers'),'dir')
        mkdir(fullfile(meta.workingfolder,'ColorTracking','Markers'))
    end
    if ~exist(fullfile(meta.workingfolder,'OpenSim'),'dir')
        mkdir(fullfile(meta.workingfolder,'OpenSim'))
    end
    if ~exist(fullfile(meta.workingfolder,'OpenSim','Analysis'),'dir')
        mkdir(fullfile(meta.workingfolder,'OpenSim','Analysis'))
    end
end

%% Merge and strip files for spike sorting
% Run processSpikesForSorting for the first time to combine spike data from
% all files with a name starting with file_prefix.
if ~meta.sorted
    processSpikesForSorting(fullfile(meta.workingfolder,'preCDS','merging'),meta.neuralPrefix);
    if exist('altMeta','var') && ~isempty(altMeta.array)
        processSpikesForSorting(fullfile(altMeta.workingfolder,'preCDS','merging'),altMeta.neuralPrefix);
    end
end

% Now sort in Offline Sorter!

%% Load colorTracking file (and settings if desired) -- NOTE: Can do this simultaneously with sorting, since it takes some time
if meta.motionTracked
    if ~meta.markered
        first_time = 1;
        for fileIdx = 1:length(meta.taskAlias)
            colorTrackingFilename = [meta.monkey '_' meta.date '_colorTracking_' meta.taskAlias{fileIdx}];

            fname_load=ls(fullfile(meta.workingfolder,'ColorTracking',[colorTrackingFilename '*']));
            load(deblank(fullfile(meta.workingfolder,'ColorTracking',fname_load)))

            % Run color tracking script
            color_tracker_4colors_script;

            % Save
            markersFilename = [meta.monkey '_' meta.date '_markers_' meta.taskAlias{fileIdx}];
            fname_save=fullfile(meta.workingfolder,'ColorTracking','Markers',[markersFilename '.mat']);
            save(fname_save,'all_medians','all_medians2','led_vals','times');

            if first_time
                fname_save_settings=fullfile(meta.workingfolder,'ColorTracking','Markers',['settings_' meta.monkey '_' meta.date]);
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
        meta.markered = true;
    else
        fprintf('Please move marker data into working directory')
        winopen(fullfile(meta.workingfolder,'ColorTracking','Markers'))
        winopen(meta.markersfolder)
        pause
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except meta altMeta

%% Split files and move to Final folder before loading
if ~meta.sorted
    % pause to let the user finish sorting
    pause
    
    % split things up
    processSpikesForSorting(fullfile(meta.workingfolder,'preCDS','merging'),meta.neuralPrefix);
    if exist('altMeta','var') && ~isempty(altMeta.array)
        processSpikesForSorting(fullfile(altMeta.workingfolder,'preCDS','merging'),altMeta.neuralPrefix);
    end

    % move into final folder
    for fileIdx = 1:length(meta.taskAlias)
        movefile(fullfile(meta.workingfolder,'preCDS','merging',[meta.neuralPrefix '_' meta.taskAlias{fileIdx} '.mat']),...
            fullfile(meta.workingfolder,'preCDS','Final'));
        if exist('altMeta','var') && ~isempty(altMeta.array)
            movefile(fullfile(altMeta.workingfolder,'preCDS','merging',[altMeta.neuralPrefix '_' altMeta.taskAlias{fileIdx} '.mat']),...
                fullfile(altMeta.workingfolder,'preCDS','Final'));
        end
        movefile(fullfile(meta.workingfolder,'preCDS','merging',[meta.neuralPrefix '_' meta.taskAlias{fileIdx} '.nev']),...
            fullfile(meta.workingfolder,'preCDS','Final'));
        if exist('altMeta','var') && ~isempty(altMeta.array)
            movefile(fullfile(altMeta.workingfolder,'preCDS','merging',[altMeta.neuralPrefix '_' altMeta.taskAlias{fileIdx} '.nev']),...
                fullfile(altMeta.workingfolder,'preCDS','Final'));
        end
    end
    meta.sorted = true;
else
    % open windows to move sorted stuff into place
    fprintf('Please move sorted files into working directory')
    winopen(fullfile(meta.workingfolder,'preCDS','Final'))
    winopen(meta.sortedfolder)
    pause
end

%% Load data into CDS file
% Make CDS files
cds_cell = cell(size(meta.taskAlias));
for fileIdx = 1:length(meta.taskAlias)
    cds_cell{fileIdx} = commonDataStructure();
    cds_cell{fileIdx}.file2cds(...
        fullfile(meta.workingfolder,'preCDS','Final',[meta.neuralPrefix '_' meta.taskAlias{fileIdx}]),...
        ['ranBy' meta.ranBy],...
        ['array' meta.array],...
        ['monkey' meta.monkey],...
        meta.lab,...
        'ignoreJumps',...
        'unsanitizedTimes',...
        ['task' meta.task{fileIdx}],...
        ['mapFile' meta.mapfile]);

    % also load second file if necessary
    if exist('altMeta','var')
        if ~isempty(altMeta.array)
            cds_cell{fileIdx}.file2cds(...
                fullfile(altMeta.workingfolder,'preCDS','Final',[altMeta.neuralPrefix '_' altMeta.taskAlias{fileIdx}]),...
                ['ranBy' altMeta.ranBy],...
                ['array' altMeta.array],...
                ['monkey' altMeta.monkey],...
                altMeta.lab,...
                'ignoreJumps',...
                'unsanitizedTimes',...
                ['task' altMeta.task{fileIdx}],...
                ['mapFile' altMeta.mapfile]);
        else
            cds_cell{fileIdx}.file2cds(...
                fullfile(altMeta.workingfolder,'preCDS','Final',[altMeta.neuralPrefix '_' altMeta.taskAlias{fileIdx}]),...
                ['ranBy' altMeta.ranBy],...
                ['monkey' altMeta.monkey],...
                altMeta.lab,...
                'ignoreJumps',...
                'unsanitizedTimes',...
                ['task' altMeta.task{fileIdx}]);
        end
    end
end

%% Load marker file and create TRC
if meta.markered
    for fileIdx = 1:length(meta.taskAlias)
        markersFilename = [meta.monkey '_' meta.date '_markers_' meta.taskAlias{fileIdx} '.mat'];
        affine_xform = cds_cell{fileIdx}.loadRawMarkerData(fullfile(meta.workingfolder,'ColorTracking','Markers',markersFilename));
        writeTRCfromCDS(cds_cell{fileIdx},fullfile(meta.workingfolder,'OpenSim',[meta.monkey '_' meta.date '_' meta.taskAlias{fileIdx} '_markerData.trc']))
        % writeHandleForceFromCDS(cds_cell{fileIdx},fullfile(meta.workingfolder,'OpenSim',[meta.monkey '_' meta.date '_' meta.taskAlias{fileIdx} '_handleForce.mot']))
    end
end

%% Do openSim stuff and save analysis results to analysis folder
% pause for user to run opensim stuff
if meta.markered
    fprintf('Run OpenSim stuff!')
    winopen(fullfile(meta.workingfolder,'OpenSim'))
    winopen(meta.opensimsettingsfolder)
    pause
end

%% Add kinematic information to CDS
if meta.markered
    for fileIdx = 1:length(meta.taskAlias)
        % load joint information
        cds_cell{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'joint_ang')

        % load joint velocities
        cds_cell{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'joint_vel')

        % load joint moments
        % cds{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'joint_dyn')

        % load muscle information
        cds_cell{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'muscle_len')

        % load muscle velocities
        cds_cell{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'muscle_vel')

        % load hand positions
        cds_cell{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'hand_pos')

        % load hand velocities
        cds_cell{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'hand_vel')

        % load hand accelerations
        cds_cell{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'hand_acc')

        % load elbow positions
        cds_cell{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'elbow_pos')

        % load elbow velocities
        cds_cell{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'elbow_vel')

        % load elbow accelerations
        cds_cell{fileIdx}.loadOpenSimData(fullfile(meta.workingfolder,'OpenSim','Analysis'),'elbow_acc')
    end
end

%% Save CDS
for fileIdx = 1:length(meta.taskAlias)
    cds_name = sprintf('%s_%s_CDS_%s.mat',meta.monkey,meta.date,meta.taskAlias{fileIdx});
    cds = cds_cell{fileIdx};
    save(fullfile(meta.cdslibrary,cds_name),'cds','-v7.3')
end

%% Copy stuff back to server for backup
if meta.markered
    winopen(meta.markersfolder)
    winopen(fullfile(meta.workingfolder,'ColorTracking','Markers'))
end

if meta.sorted
    winopen(meta.sortedfolder)
    winopen(fullfile(meta.workingfolder,'preCDS'))
end

winopen(fullfile(meta.remotefolder,'cds-library'))
winopen(meta.cdslibrary)
