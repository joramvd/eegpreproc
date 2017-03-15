%% Preprocessing EEG I: steps prior to Independent Component Analysis
% Joram van Driel, VU Amsterdam, July 2016

%% setup preliminaries

% Paths
clear, close all, warning('off','all')
addpath(genpath('path\to\eeglab'));
addpath('path\to\fieldtrip\'); ft_defaults; % be aware of toolbox conflicts
 
% Data dirs
basedir = 'path\to\EEG\';
readdir = [basedir 'raw\'];
writdir = [basedir 'processed\'];

% length of time for EEG epoch (in sec) 
epochtime=[ -1.5 4 ]; % trial lasts around 3 sec. including post-response window; adding 1sec buffer zone around it

cd(readdir)
sublist=dir('*.set');
sublist={sublist.name};

%% Loop around subjects

for subno=1:length(sublist)
    
    outfilename=[ writdir sublist{subno}(1:end-4) '_preReject.mat' ]; % this code assumes a sensible basename of the raw .bdf file
    
    % skip if already exist
    if exist(outfilename,'file')
        continue
    end
    
    %% load data
    
    % read bdf file
    fprintf('Loading subject %i of %i\n',subno,length(sublist));
    EEG = pop_loadset(sublist{subno});
    
    % Rereferencing and computing bipolar signals (EOG/EMG) requires some
    % juggling with external electrode labels; this of course depends on
    % your own settings, so please adapt accordingly!
    
    % re-reference to linked mastoids/earlobes
    EEG = pop_reref(EEG,[find(strcmpi('EXG1',{EEG.chanlocs.labels})) find(strcmpi('EXG2',{EEG.chanlocs.labels}))],'refstate','averef');

    % re-reference EOG
    EEG.data(strcmpi({EEG.chanlocs.labels},'EXG3'),:) = squeeze(EEG.data(strcmpi({EEG.chanlocs.labels},'EXG3'),:)-EEG.data(strcmpi({EEG.chanlocs.labels},'EXG4'),:));
    EEG.data(strcmpi({EEG.chanlocs.labels},'EXG5'),:) = squeeze(EEG.data(strcmpi({EEG.chanlocs.labels},'EXG5'),:)-EEG.data(strcmpi({EEG.chanlocs.labels},'EXG6'),:));
    EEG.data(strcmpi({EEG.chanlocs.labels},'EXG7'),:) = squeeze(EEG.data(strcmpi({EEG.chanlocs.labels},'EXG7'),:)-EEG.data(strcmpi({EEG.chanlocs.labels},'EXG8'),:));
    
    EEG.chanlocs(strcmpi({EEG.chanlocs.labels},'EXG3')).labels = 'HEOG';
    EEG.chanlocs(strcmpi({EEG.chanlocs.labels},'EXG5')).labels = 'EMGR';
    EEG.chanlocs(strcmpi({EEG.chanlocs.labels},'EXG7')).labels = 'EMGL';
    
    % remove unnecessary channels
    EEG = pop_select(EEG,'nochannel',[ find(strcmpi({EEG.chanlocs.labels},'EXG1')) find(strcmpi({EEG.chanlocs.labels},'EXG2')) find(strcmpi({EEG.chanlocs.labels},'EXG4')) find(strcmpi({EEG.chanlocs.labels},'EXG6'))  find(strcmpi({EEG.chanlocs.labels},'EXG8')) ]);
    EEG = pop_select(EEG,'nochannel',[ find(strcmpi({EEG.chanlocs.labels},'Erg1')) find(strcmpi({EEG.chanlocs.labels},'Erg2')) find(strcmpi({EEG.chanlocs.labels},'GSR1')) find(strcmpi({EEG.chanlocs.labels},'GSR2')) find(strcmpi({EEG.chanlocs.labels},'Resp')) find(strcmpi({EEG.chanlocs.labels},'Plet')) find(strcmpi({EEG.chanlocs.labels},'Temp')) ]);
    
    % (high-pass) filter; quite fast with eegfiltnew; make sure you have an eeglab version that has this updated function
    EEG = pop_eegfiltnew(EEG,.5,0);
    
    %%
    
    % Again very experiment-specific: the triggers that you want to lock
    % your data to; the condition-names are just for convenience
    % You (probably) shouldn't put all the different triggers that you used
    % in here! For example, suppose a trial consists of multiple stimuli,
    % responses, displays, that you all gave unique triggers; but you want
    % to timelock only to the first display; then you should only put this
    % triggger in the variable below; the other triggers will still be
    % present as 'events' in the epochs that result from locking to the
    % triggers in the below variable.
    triggers={
        'conditionA'    { '1' };
        'conditionB'    { '2' };
        };
    
    % Epoch the data, locked to the above triggers; 'epochtime' is the time
    % window that you specified above, with some time before, and some time
    % after the trigger
    EEG = pop_epoch(EEG,[triggers{:,2}],[epochtime(1) epochtime(2)]);
    
    % baseline-correct
    % the empty brackets result in whole-trial baseline correction
    EEG = pop_rmbase(EEG,[]);
    
    EEG.setname=sublist{subno};
    
    save(outfilename,'EEG')
    
    
end

%% manual trial rejection
clear, close all

% Rekenbeest dirs
basedir = '\path_to\EEG\';
rawdir  = [basedir 'raw\'];        % we need both the raw files
writdir = [basedir 'processed\'];  % as well as the preproceses files from the steps above

cd(writdir)
procfilz = dir('*preReject.mat');
rawfilz  = dir([rawdir '*.bdf']);

%% Loop through raw files to automatically detect artifacts

% Below script makes the trial rejection procedure fully automatic;
% To have it be semi-automatic (i.e. inspect the results of the algorithm
% for each subjects' raw data), uncomment-out the pieces of eeglab code!

for filei = 1:length(procfilz)
    %%
    
    outfilename = ([writdir procfilz(filei).name(1:end-13) 'rejected_trials.mat']); % you should check whether the .name(1:end-13) produces a sensible filename!
    if exist(outfilename,'file'), error('file already exists!'); end
    
%     load(filz(filei).name);
%    
%    % backup
%     origEEG=EEG;
%     ALLEEG=EEG;
%     eeglab redraw
%     pop_eegplot(EEG,1,1,0); % plot EEG to inspect raw data
    
    % use fieldtrip to detect artifacts
    
    cfg=[];
    cfg.filename              = [rawdir rawfilz(filei).name];
    cfg.headerfile            = cfg.filename;
    cfg.trialdef.eventtype    = 'STATUS';
    cfg.trialdef.eventvalue   = 11:67; % change these values according to your own trigger scheme! use the same triggers to which you locked your epochs (i.e. don't use all triggers if you have multiple triggers in a trial)
    cfg.trialdef.prestim      = 0; % latency in seconds
    cfg.trialdef.poststim     = 2; % latency in seconds
    
    cfg = ft_definetrial(cfg);
    trl = cfg.trl;
        
    %% Automatic artifact detection algorithm customized to detect muscle artifacts
    
    filename       = [rawdir rawfilz(filei).name];
    cfg            = [];
    cfg.trl        = trl;
    cfg.datafile   = filename;
    cfg.headerfile = filename;
    cfg.continuous = 'yes';
    
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel     = 1:64;  % or 128, or 32; check your EEG system
    cfg.artfctdef.zvalue.flexible    = 'yes'; % if 'yes', the below value is added to the 'bandwidth' of z-values (i.e. the cut-off will change per subject); if 'no', the below value is the actual cut-off
    cfg.artfctdef.zvalue.cutoff      = 4;     % this value may change depending on the subject; try to stick to one value for entire data set!
    cfg.artfctdef.zvalue.trlpadding  = 0.5;   % how much time around the trial are artifacts still not acceptable?
    cfg.artfctdef.zvalue.fltpadding  = 0.5;   % this removes time due to edge artifacts of the band-pass filter
    cfg.artfctdef.zvalue.artpadding  = 0.1;   % this adds some time around the detected artifact
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.bpfilter    = 'yes';
    cfg.artfctdef.zvalue.bpfreq      = [110 140]; % specified band-pass for EMG muscle noise
    cfg.artfctdef.zvalue.bpfiltord   = 9;
    cfg.artfctdef.zvalue.bpfilttype  = 'but';
    cfg.artfctdef.zvalue.hilbert     = 'yes';
    cfg.artfctdef.zvalue.boxcar      = 0.2;
    
    % make the process interactive ('yes') or not ('no'); when you set this
    % to 'no', and leave the EEGLAB code commented-out, the subject-loop
    % will be fully automatic!
    cfg.artfctdef.zvalue.interactive = 'yes'; % just to get a quick glance at the z-value threshold
    
    [cfg] = ft_artifact_zvalue_x(cfg); % this adapted frieldtrip function (hence the _x) also stores a screenshot of the interactive gui, so you can always have a look a the range of z-values
    [cfg] = ft_rejectartifact(cfg); % do the rejection (this only removes trials from the trl matrix)
    
    rejected_trials = ismember(cfg.trlold(:,1),cfg.trl(:,1))==0; % compare old with new trial matrix to get trial indices instead of trial timings
    
%     % put these into EEGLAB format
%     EEG.reject.rejmanual = rejected_trials;
%     EEG.reject.rejmanualE = zeros(EEG.nbchan,length(trl));
%     pop_eegplot(EEG,1,1,0); % the 1,1,0 settings marks trials that were previously rejected (i.e. with fieldtrip)
%     
%     disp(find(rejected_trials)'); % display the trial indices in the command line
%
%     % in case you made changes to the fieldtrip trial selection:
%     rejected_trials = ALLEEG.reject.rejmanual;

    save(outfilename,'rejected_trials','cfg'); % the cfg variable contains the rejection settings, so these can be traced back (e.g. z-val cutoff)
    
    
end
