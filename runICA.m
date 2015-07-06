%% Preprocessing EEG II: remove detected artifact-trials, perform ICA
% Joram van Driel, VU Amsterdam, July 2016

%% setup preliminaries

clear, close all

addpath(genpath('path\to\eeglab'));
addpath('path\to\ADJUST1.1\');

% dirs
basedir = 'path\to\EEG\';
readdir = [basedir 'processed\'];
writdir = [basedir 'processed\'];

cd(readdir)

%%
filz = dir('*preReject.mat');
trialfilz = dir('*rejected_trials.mat');

%% big loop around subjects

for subno=1:length(filz)
    
    outfilename=[ writdir filz(subno).name(1:end-13) 'withICA.mat' ];
    
    % skip if already exist
    if exist(outfilename,'file')
        continue
    end
    
    fprintf('Loading datafile %i for ICA...\n',subno);
    
    % load data
    load([ readdir filz(subno).name ]);
    
    %% remove trials
    
    % this is the trial-index vector from the (semi-)automatic rejection
    % procedure of the preICA.m script
    load(trialfilz(subno).name);
    
    % now the trials are actually removed from the data
    EEG=pop_select(EEG,'notrial',find(rejected_trials));
    
    disp(' ')
    disp([ 'Subject ' num2str(subno) ' had ' num2str(sum(rejected_trials)) ' trials removed (' num2str(100*sum(rejected_trials)/length(rejected_trials)) '% of trials).' ])

    %% import channels to be interpolated (this is done later, after ICA; ICA not applied on to-be-interpolated channels!)
    
    % Based on visual inspection of the data, or because you know a channel
    % broke down during recording, you may have created a txt file with
    % channel names noted down for each subject that eventually need to be
    % interpolated
    % Because ICA needs to be done on data that are as clean as possible at
    % this stage, and because interpolated channels cause 'reduced ranking'
    % of the ICA matrix, it is best to 
    
    fid=fopen('chans2interp.txt','r');
    chans2interp={};
    while ~feof(fid)
        aline=regexp(fgetl(fid),'\t','split');
        chans2interp{size(chans2interp,1)+1,1}=aline{1};
        for i=2:length(aline)
            chans2interp{size(chans2interp,1),i}=aline{i};
        end
    end
    
    chanind=1:64;
    subject_prefix = filz(subno).name(end-16:end-14);
    chans = chans2interp(strcmpi(chans2interp(:,1),subject_prefix),2:end);
    if ~isempty(chans{1})
        bad_chansidx = zeros(1,length(chans));
        for c=1:length(chans)
            if ~isempty(chans{c}), bad_chansidx(c) = find(strcmpi({EEG.chanlocs.labels},chans{c})); end
        end
        bad_chansidx(bad_chansidx==0)=[];
        chanind(bad_chansidx)=[];
    end
    
    %% run ICA
    EEG=pop_runica(EEG,'icatype','runica','extended', 1,'chanind',chanind); % don't include EOG and bad channels
    
	% write data to disk
    % save ICA weights separately (everything else is same as before ICA;
    % saves space and loading time)
    
    ICAEEG.icawinv = EEG.icawinv;
    ICAEEG.icasphere = EEG.icasphere;
    ICAEEG.icaweights = EEG.icaweights;
    ICAEEG.icachansind = EEG.icachansind;

    save(outfilename,'ICAEEG')
            
%% wait for it...
end




%% manually check ICA result

eeglab % opens the EEGLAB gui

filz = dir('*preReject.mat');
icafilz = dir('*withICA.mat');

subno = 1;


% load data

load([ readdir filz(subno).name ]);
load([ readdir icafilz(subno).name]);

EEG.icachansind = ICAEEG.icachansind;
EEG.icasphere = ICAEEG.icasphere;
EEG.icaweights = ICAEEG.icaweights;
EEG.icawinv = ICAEEG.icawinv;

eeglab redraw

% Now in the gui, click on Tools > Reject data using ICA > Reject
% components by map
% The next preprocessing script assumes you note down the blink-components
% manually in a .txt file





%% Automatically check artifact-components using the ADJUST plugin

filz = dir('*preReject.mat');
icafilz = dir('*withICA.mat');

fid=fopen('chans2interp.txt','r');
chans2interp={};
while ~feof(fid)
    aline=regexp(fgetl(fid),'\t','split');
    chans2interp{size(chans2interp,1)+1,1}=aline{1};
    for i=2:length(aline)
        chans2interp{size(chans2interp,1),i}=aline{i};
    end
end

%% check the components

for subno=1:length(icafilz)
    
    outfilename=[ writdir icafilz(subno).name(1:end-12) 'ADJUST.txt' ];
    
    % skip if already exist
    if exist(outfilename,'file')
        continue
    end
    
    fprintf('Loading datafile %i for ADJUST...\n',subno);
    
    % load data
    
    load([ readdir filz(subno).name ]);
    load([ readdir icafilz(subno).name]);
    
    EEG.icachansind = ICAEEG.icachansind;
    EEG.icasphere = ICAEEG.icasphere;
    EEG.icaweights = ICAEEG.icaweights;
    EEG.icawinv = ICAEEG.icawinv;
    
    %% Mark bad channels
    
    bad_chansidx=0;
    subject_prefix = filz(subno).name(end-16:end-14);
    chans = chans2interp(strcmpi(chans2interp(:,1),subject_prefix),2:end);
    if ~isempty(chans{1})
        bad_chansidx = zeros(1,length(chans));
        for c=1:length(chans)
            if ~isempty(chans{c}), bad_chansidx(c) = find(strcmpi({EEG.chanlocs.labels},chans{c})); end
        end
        bad_chansidx(bad_chansidx==0)=[];
    end
    
    
    %% apply ADJUST
    
    % remove bad channel indices and EOG (i.e. the channels that were not
    % included in the ICA)
    % This removal is not saved, but is necessary for the ADJUST plugin to
    % work
    
    if sum(bad_chansidx)>0
        EEG = pop_select(EEG,'nochannel',[65:66 bad_chansidx]);
    else
        EEG = pop_select(EEG,'nochannel',[65:66]);
    end        

    [art, horiz, vert, blink, disc]=ADJUST(EEG,outfilename);

end


