%% set data folder
cdslib = 'C:\Users\rhc307\data\cds-library';

% files = dir(fullfile(cdslib,'*COactpas*'));
% names = vertcat({files.name})';
names = {...
    'Han_20170203_OOR25N_CDS_001.mat',...
    };

load_params = struct(...
    'array_name','S1',...
    'cds_array_name','LeftS1Area2',...
    'cont_signal_names',{{...
        'pos',...
        'vel',...
        'force',...
        }},...
    'event_names',{{...
        'startTime',...
        'endTime',...
        'startTargOnTime',...
        'startTargHoldTime',...
        'goCueTime',...
        'endTargHoldTime',...
        }},...
    'trial_meta',{{...
        'tgtDir',...
        'forceDir',...
        }},...
    'bin_size',0.01);

td_cell = cell(length(names),1);
for i = 1:length(names)
    td_cell{i} = loadTDfromCDS(fullfile(cdslib,names{i}),load_params);
    fprintf('File %d processed\n',i)
end

%%
savedir = 'C:\Users\rhc307\data\project-data\limblab\misc\td-library';

file_info_cell = cell(length(names),6);

for i = 1:length(names)
    file_info_cell(i,:) = strsplit(names{i},{'_','.'});
end

file_info = struct(...
    'monkey',file_info_cell(:,1),...
    'date',file_info_cell(:,2),...
    'task',file_info_cell(:,3),...
    'filenum',file_info_cell(:,5));

for i = 1:length(td_cell)
    trial_data = td_cell{i};
    save(fullfile(savedir,sprintf('%s_%s_%s_TD.mat',file_info(i).monkey,file_info(i).date,file_info(i).task)),'trial_data')
    fprintf('File %d saved\n',i)
end