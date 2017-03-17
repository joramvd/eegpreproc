%% Preprocessing EEG III: do the final cleaning and divide into conditions
% Joram van Driel, VU Amsterdam, July 2016

%% setup preliminaries

clear, close all

addpath(genpath('path\to\eeglab'));
addpath('path\to\ADJUST1.1\');
addpath(genpath('path\to\erplab_4.0.3.1'));

datdir = 'path\to\EEG\processed\';

cd(datdir)

prefilz = dir('*preReject*');
icafilz = dir('*withICA*');
trlfilz = dir('*rejected_trials.mat');

use_adjust = false;

%% import ICs to be removed (alternatively: use ADJUST plugin)

fid=fopen('ICs2remove.txt','r');
comps2remove={};
while ~feof(fid)
    aline=regexp(fgetl(fid),'\t','split');
    comps2remove{size(comps2remove,1)+1,1}=aline{1};
    comps2remove{size(comps2remove,1)  ,2}=sscanf(aline{2},'%g');
end

%% import channels to be interpolated

fid=fopen('chans2interp.txt','r');
chans2interp={};
while ~feof(fid)
    aline=regexp(fgetl(fid),'\t','split');
    chans2interp{size(chans2interp,1)+1,1}=aline{1};
    for i=2:length(aline)
        chans2interp{size(chans2interp,1),i}=aline{i};
    end
end

%%
for subno=1:length(prefilz)
    
    outfilename = [icafilz(subno).name(1:end-11) 'CSD.mat'];
    if exist(outfilename,'file'), continue; end
    
    fprintf('Preprocessing subject %i of %i\n',subno,length(prefilz))
    
    %%
    load(prefilz(subno).name);
    load(icafilz(subno).name);
    
    % save peri-trial info
    
    for ei=2:length(EEG.epoch)-1
        EEG.epoch(ei).nextinfo = EEG.epoch(ei+1);
        EEG.epoch(ei).previnfo = EEG.epoch(ei-1);
    end
    new_epoch=EEG.epoch; % save epoch structure
    
    % remove detectec trials with artifacts
    
    load(trlfilz(subno).name);
    
    EEG=pop_select(EEG,'notrial',find(rejected_trials));
    
    % update trial info in epochs
    new_epoch(rejected_trials)=[];
    EEG.epoch=new_epoch;
    
    disp(' ')
    disp([ 'Subject ' num2str(subno) ' had ' num2str(sum(rejected_trials)) ' trials removed (' num2str(100*sum(rejected_trials)/length(rejected_trials)) '% of trials).' ])
    
    % add ICA weights to EEG structure
    
    EEG.icachansind = ICAEEG.icachansind;
    EEG.icasphere = ICAEEG.icasphere;
    EEG.icaweights = ICAEEG.icaweights;
    EEG.icawinv = ICAEEG.icawinv;
    clear ICAEEG
    
    % if subject contains bad channels:
    % - first remove channels from data, then remove blink ICs, then
    % - interpolate removed channels
    
    bad_chansidx=0;
    subject_prefix = prefilz(subno).name(end-16:end-14);
    chans = chans2interp(strcmpi(chans2interp(:,1),subject_prefix),2:end);
    if ~isempty(chans{1})
        bad_chansidx = zeros(1,length(chans));
        for c=1:length(chans)
            if ~isempty(chans{c}), bad_chansidx(c) = find(strcmpi({EEG.chanlocs.labels},chans{c})); end
        end
        bad_chansidx(bad_chansidx==0)=[];
    end
    
    if sum(bad_chansidx)>0
        EEG2 = pop_select(EEG,'nochannel',[65:66 bad_chansidx]);
    else
        EEG2 = pop_select(EEG,'nochannel',[65:66]);
    end
    
    % remove IC blink components from ADJUST output (it may sometimes
    % generate errors due to negative values...)
    
    if use_adjust
        try
            adjustname=[ datdir icafilz(subno).name(1:end-12) '_ADJUST.txt' ];
            [art, horiz, vert, blink, disc]=ADJUST(EEG2,adjustname,2); % 2 is added to multiply the automatically determined thresholds by 3; otherwise too many components were selected
        catch me2; continue
        end
    else
        blink = comps2remove{subno,2}';
    end
    
    %%
    
    EEG2 = pop_subcomp(EEG2,blink(blink<21)); % only remove detected blink-ICs lower than 20; higher IC# contribute less variance, and are likely to reflect noise in frontal electrodes rather than blinks
    
    %% put IC-cleaned channels back in original EEG structure and interpolate bad ones
    
    good_chans = 1:64;
    if sum(bad_chansidx)>0
        good_chans(bad_chansidx)=[];
        EEG.data(good_chans,:,:) = EEG2.data;
        EEG = eeg_interp(EEG,bad_chansidx);
    else
        EEG.data(good_chans,:,:) = EEG2.data;
    end
    
    clear EEG2
    
    %% horizontal eye-movements
    
    rejected_trials_heog = pop_artstepX(EEG,'channel',66,'flag',1,'threshold',25,'Twindow', [ -50 1050], 'Windowsize',  100, 'Windowstep',  50 );
    save([ prefilz(subno).name(1:end-13) 'rejected_trials_heog.mat'],'rejected_trials_heog');
    EEG = pop_select(EEG,'notrial',find(rejected_trials_heog));
    
%     EEG.reject.rejmanual = rejected_trials_heog;
%     EEG.reject.rejmanualE = zeros(EEG.nbchan,length(EEG.epoch));
% 
%     pop_eegplot(EEG,1,1,0); % the 1,1,0 settings marks trials that were previously rejected (i.e. with pop_artstepX)

    
    %% remove EOG for good
    
    EEG = pop_select(EEG,'nochannel',65:66);
    
    %% Group epochs into conditions
    
    % -- Correct, easy vs diff search        % -- Same for errors
    % SimpleLacc1    .{11-17 } {81}
    % SimpleRacc1    .{21-27 } {81}
    % DistLacc1      .{31-37 } {81}          % DistLacc0      .{31-37 } {80}
    % DistRacc1      .{41-47 } {81}          % DistRacc0      .{41-47 } {80}
    % NonDistLacc1   .{51-57 } {81}          % NonDistLacc0   .{51-57 } {80}
    % NonDistRacc1   .{61-67 } {81}          % NonDistRacc0   .{61-67 } {80}
    
    try
        for ei=1:length(EEG.epoch)
            
            cueonset = find(cell2mat(EEG.epoch(ei).eventlatency)==0);
            if cell2mat(EEG.epoch(ei).eventtype(cueonset)) < 18 && cell2mat(EEG.epoch(ei).eventtype(cueonset)) >10 %&& cell2mat(EEG.epoch(ei).eventtype(2)) == 81
                EEG.epoch(ei).conname = 'SimpleL';
            elseif cell2mat(EEG.epoch(ei).eventtype(cueonset)) < 28 && cell2mat(EEG.epoch(ei).eventtype(cueonset)) >20 %&& cell2mat(EEG.epoch(ei).eventtype(2)) == 81
                EEG.epoch(ei).conname = 'SimpleR';
            elseif cell2mat(EEG.epoch(ei).eventtype(cueonset)) < 38 && cell2mat(EEG.epoch(ei).eventtype(cueonset)) >30 %&& cell2mat(EEG.epoch(ei).eventtype(2)) == 81
                EEG.epoch(ei).conname = 'DistL';
            elseif cell2mat(EEG.epoch(ei).eventtype(cueonset)) < 48 && cell2mat(EEG.epoch(ei).eventtype(cueonset)) >40 %&& cell2mat(EEG.epoch(ei).eventtype(2)) == 81
                EEG.epoch(ei).conname = 'DistR';
            elseif cell2mat(EEG.epoch(ei).eventtype(cueonset)) < 58 && cell2mat(EEG.epoch(ei).eventtype(cueonset)) >50 %&& cell2mat(EEG.epoch(ei).eventtype(2)) == 81
                EEG.epoch(ei).conname = 'NonDistL';
            elseif cell2mat(EEG.epoch(ei).eventtype(cueonset)) < 68 && cell2mat(EEG.epoch(ei).eventtype(cueonset)) >60 %&& cell2mat(EEG.epoch(ei).eventtype(2)) == 81
                EEG.epoch(ei).conname = 'NonDistR';
            end
        end
    catch me, continue;
    end
    
    connames = {'SimpleL','SimpleR','DistL', 'DistR', 'NonDistL', 'NonDistR'};%, 'DistLacc0', 'DistRacc0', 'NonDistLacc0', 'NonDistRacc0'};
    ALLEEG = EEG;
    for condi=1:length(connames)
        ALLEEG(condi) = EEG;
        ALLEEG(condi).data = ALLEEG(condi).data(:,:,(strcmpi({EEG.epoch.conname}, connames{condi})));
        ALLEEG(condi).epoch = ALLEEG(condi).epoch((strcmpi({EEG.epoch.conname}, connames{condi})));
        ALLEEG(condi).trials=size(ALLEEG(condi).data,3);
        ALLEEG(condi).setname=connames{condi};
    end
    
    %% CSD data
    
    load GHmatrices.mat
    for ci=1:length(ALLEEG)
        fprintf('Applying CSD on condition %i of %i\n',ci,length(connames))
        ALLEEG(ci).data = CSD(ALLEEG(ci).data,G,H);
    end
    
    % save CSD-epoched data for this subject!
    save(outfilename, 'ALLEEG');
    
end
