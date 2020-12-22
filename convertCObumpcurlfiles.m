%% set data folder
if ispc
    dataroot = 'C:\Users\rhc307\data';
else
    dataroot = '/data/raeed';
end
cdslib = fullfile(dataroot,'cds-library');

% files = dir(fullfile(cdslib,'*COactpas*'));
% names = vertcat({files.name})';
names = {...
    'Han_20171206_CObumpcurl_baseline_CDS_001.mat',...
    'Han_20171206_CObumpcurl_adaptation_CDS_002.mat',...
    'Han_20171206_CObumpcurl_washout_CDS_003.mat',...
    };

load_params = struct(...
    'array_name','S1',...
    'cds_array_name','LeftS1Area2',...
    'extract_emg',true,...
    'cont_signal_names',{{...
        'pos',...
	    'vel',...
        'force',...
	    'markers',...
	    'joint_ang',...
	    'joint_vel',...
	    'muscle_len',...
	    'muscle_vel',...
        }},...
    'event_names',{{...
        'startTime',...
        'endTime',...
        'goCueTime',...
        'bumpTime',...
        'tgtOnTime',...
        }},...
    'trial_meta',{{...
        'bumpDir',...
        'tgtDir',...
        'ctrHoldBump',...
        'ctrHold',...
        }});

td_cell = cell(length(names),1);
for i = 1:length(names)
    td_cell{i} = loadTDfromCDS(fullfile(cdslib,names{i}),load_params);
    fprintf('File %d processed\n',i)
end

trial_data = horzcat(td_cell{:});

%%
savedir = fullfile(dataroot,'project-data','limblab','s1-adapt','td-library');
savename = 'H_20171206_CObumpcurl_TD.mat';
save(fullfile(savedir,savename),'trial_data')

% file_info_cell = cell(length(names),6);
% 
% for i = 1:length(names)
%     file_info_cell(i,:) = strsplit(names{i},{'_','.'});
% end
% 
% file_info = struct(...
%     'monkey',file_info_cell(:,1),...
%     'date',file_info_cell(:,2),...
%     'task',file_info_cell(:,3),...
%     'filenum',file_info_cell(:,5));
% 
% for i = 1:length(td_cell)
%     trial_data = td_cell{i};
%     save(fullfile(savedir,sprintf('%s_%s_%s_TD.mat',file_info(i).monkey,file_info(i).date,file_info(i).task)),'trial_data')
%     fprintf('File %d saved\n',i)
% end
