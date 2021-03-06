function eegpreproc(cfg)

%- Preprocessing EEG
%  Joram van Driel
%  Vrije Universiteit Amsterdam
%  March 2017
%
% This function applies the following four steps to Biosemi raw data (bdf files)
% - 1) it imports the raw data, references the data to mastoids/earlobes, imports
%      predefined channel layout, and applies a high-pass filter if recorded on >512 Hz sampling
%      rate, the data are downsampled
% - 2) it epochs the data around [stimulus] triggers, performs [whole-trial] baseline
%      correction, and marks trials contaminated with EMG [muscle] artifacts (note: it does not
%      remove the trials yet);
%      it also saves a figure with snapshots of topographical maps showing electrode variance; this
%      is to inspect bad electrodes that should later be interpolated
% - 3) performs ICA using eeglab
% - 4) does the final cleaning by removing the bad trials, removing occulomotor
%      artifacts captured by ICA, interpolating bad channels  (note: if interested in
%      trial-sequence-specific information [e.g. switch-trials, correct after error trials; i.e.
%      used for single-trial regression] one needs to adapt this part on a study-specific basis, as
%      here the trials marked for rejection are actually removed, thus destroying the original trial
%      sequence)
%
% At every step a datafile is saved; for trial rejection this only contains a binary vector;  for
% ICA (step #3) this only contains the weights. Datafile naming assumes the first four charachters
% correspond to subject number: pp## (where pp refers to Dutch "proefpersoon").
%
% The pipeline "forks" at the first stage: one EEG version is high-pass filtered the requested
% cut-off (e.g. 0.1 Hz) and only saved for the later final cleaning step; another version is
% filtered at a higher cut-off (e.g. 1.5 Hz; 2.5 max, recommendation: between 1-2 Hz), which aids
% data visualization and improves ICA. All choices (trials to reject, channels to interpolate,
% components to remove) are based on the 1.5 Hz version; the actual cleaning is done on the (e.g.)
% 0.1 Hz version. It is not recommended to final analyses on 1.5 Hz filtered data, especially for
% ERP analyses 0.1 or lower is recommended.
%
% If you will do time-frequency analysis, choose a rather wide epoch window to avoid edge artifacts
% (rule of thumb: one second longer at each side); the 1.5 Hz filtered version for artifact
% detection will be done on a smaller window, which corresponds to your final time window of
% interest. This is specified in the cfg structure (see below).
%
% Artifact rejection is done "semi-automatically": it uses an algorithm that checks for
% above-threshold deflections of data filtered in a high-frequency band (110-140 Hz); code is
% borrowed from the Fieldtrip routines. The result is included in the EEGLAB structure and plotted,
% so you can manually go over the automatic detection, and add/remove where you disagree. 
%
% This function runs each step consecutively within one subject, and skips previous steps if those
% were already saved.
% After step #2 you need to manually specifiy which channels to interpolate and after step #3 you
% need to manually specificy which independent components to remove in custom .txt files, which the
% next step read in. The function is paused for this and waits for input.
%
% This function can only run with properly added paths to eeglab and fieldtrip. The eeglab package
% should at least contain the function eegfiltnew (included at least in v12 and up)
%
% This function requires as input a cfg (inspired by how Fieldtrip handles functions)
% This cfg should contain the following study-specific info:
% - paths to eeglab and fieldtrip packages
% - paths to raw files
% - path to output directory
% - path to channel layout
% - channel labels for referencing
% - channel labels that correspond to the EOG
% - new sampling rate if recorded on >512 (preferably a fraction of the original sr)
% - triggers for epoching
% - time window around which to epoch
% - time window for artifact rejection
% - filename of .txt file with channels to interpolate
% - filename of .txt file with components to remove
%
% For example:
%
% %-% part I: loading raw data and filter
% cfg = [];
% cfg.eeglab_path = 'Z:\Toolboxes\eeglab12_0_2_3b';
% cfg.ft_path     = 'Z:\Toolboxes\fieldtrip-20150318';
% cfg.readdir     = 'Z:\TempOrient\EEG\raw\';
% cfg.writdir     = 'Z:\Stuff\eegpreproc\';
% cfg.layout      = [cfg.eeglab_path '\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp'];
% cfg.projectname = 'temporient';
%
% cfg.nchan       = 64;
% cfg.reref       = {'EXG5','EXG6'};
% cfg.veog        = {'EXG1','EXG2'};
% cfg.heog        = {'EXG3','EXG4'};
% cfg.resrate     = 512;
% cfg.highcut     = 0.1; % desired high-pass of final analyses
% cfg.icacut      = 1.5; % temporary high-pass for ICA and cleaning
%
% %-% part II: epoching and trial rejection marking
%
% triggers={
%     'face'        { '10' };
%     'house'       { '20' };
%     'letter'      { '30' };
%     };
%
% cfg.epochtime   = [-1.5 2.5];
% cfg.artiftime   = [-0.5 1.5];
% cfg.artcutoff   = 12;
% cfg.triggers    = triggers;
% cfg.trigger_subtract = []; % depending on physical lab settings, sometimes weird high numbers get added to the trigger values specified in your experiment script
% cfg.inspect_chans = false; % if true, pauses function and shows figures with topomaps of variance, to highlight bad channels
%
% %-% part III: ICA
% cfg.chanfilename = 'chans2interp.txt';
%
% %-% part IV: final cleaning
% cfg.icafilename = 'ICs2remove.txt';
%
% %-% now run it
% eegpreproc(cfg);


%% Paths to packages


path_list_cell = regexp(path,pathsep,'Split');

if (exist(cfg.ft_path,'dir') == 7) && ~any(ismember(cfg.ft_path,path_list_cell))
    addpath(cfg.ft_path,'-begin');
end

if (exist(cfg.eeglab_path,'dir') == 7) && ~any(ismember(cfg.eeglab_path,path_list_cell))
     % addpath(genpath(eeglab_path),'-begin');
     addpath(cfg.eeglab_path);
     eeglab;
     
     path_list_cell = regexp(path,pathsep,'Split');
     % remove conflicting paths
     ft_on_eeglab = fullfile(cfg.eeglab_path,'external','fieldtrip-partial','utilities');
     if any(ismember(ft_on_eeglab,path_list_cell))
         rmpath(genpath(ft_on_eeglab));
     end
end
ft_defaults;
close all
clearvars('-except','cfg')
v2struct(cfg);

%% change dir to files and make list of subject filenames

cd(readdir)
sublist=dir('*.bdf');
sublist={sublist.name};

%% Loop around subjects

for subno=1:length(sublist)
    
    %%
    if ~exist([writdir sublist{subno}(1:4)],'dir')
        mkdir([writdir sublist{subno}(1:4)])
    end
    
    % filenames assume a sensible basename of the raw .bdf file
    outfilename1=[ writdir sublist{subno}(1:4) filesep sublist{subno}(1:4) '_' projectname '_reref_hpfilt.mat' ];
    outfilename1b=[ writdir sublist{subno}(1:4) filesep sublist{subno}(1:4) '_' projectname '_reref_hpfilt_4ica.mat' ];
    outfilename2=[ writdir sublist{subno}(1:4) filesep sublist{subno}(1:4) '_' projectname '_rejectedtrials.mat' ];
    outfilename3=[ writdir sublist{subno}(1:4) filesep sublist{subno}(1:4) '_' projectname '_icaweights.mat' ];
    outfilename4=[ writdir sublist{subno}(1:4) filesep sublist{subno}(1:4) '_' projectname '_epoched_cleaned.mat' ];
    
    % do all steps consecutively for one subject, or go to the step where
    % you left off
    
    reload = true;
    go_on = true;
    
    while go_on
        %% go to preprocessing step
        
        if ~exist(outfilename1,'file')
            
            %% load data
            
            % read bdf file
            fprintf('Loading subject %i of %i for re-referencing and high-pass filtering...\n',subno,length(sublist));
            EEG = pop_biosig(sublist{subno});
            
            %%
            % if malfunctioning/flat channels were replaced with external electrodes during measurement, replace them here
            
            if isfield(cfg,'replaced')
                for k=1:size(replaced,1)
                    EEG.data(find(strcmpi(replaced(k,1),{EEG.chanlocs.labels})),:) = EEG.data(find(strcmpi(replaced(k,2),{EEG.chanlocs.labels})),:);
                end
            end
            
            %%
            % re-reference to linked mastoids/earlobes
            if strcmp(ref,'average')
                EEG = pop_reref(EEG,[],'refstate','averef');
            else
                EEG = pop_reref(EEG,[find(strcmpi(ref(1),{EEG.chanlocs.labels})) find(strcmpi(ref(2),{EEG.chanlocs.labels}))],'refstate','averef');
            end
            
            % re-reference EOG
            veogdat = squeeze(EEG.data(strcmpi({EEG.chanlocs.labels},veog(1)),:)-EEG.data(strcmpi({EEG.chanlocs.labels},veog(2)),:));
            heogdat = squeeze(EEG.data(strcmpi({EEG.chanlocs.labels},heog(1)),:)-EEG.data(strcmpi({EEG.chanlocs.labels},heog(2)),:));
            
            EEG.data(nchan+1,:) = veogdat;
            EEG.data(nchan+2,:) = heogdat;
            
            clear veogdat heogdat
            
            % remove unnecessary channels
            EEG = pop_select(EEG,'nochannel',nchan+3:length(EEG.chanlocs));
            
            % resample if asked for
            try
                if EEG.srate>resrate
                    EEG = pop_resample(EEG,resrate);
                end
            catch me
                disp(['No resampling done. Sampling rate is ' num2str(EEG.srate)])
            end
            
            % (high-pass) filter; eegfiltnew is quite fast
            try
                EEG = pop_eegfiltnew(EEG,highcut,0);
                EEGb = pop_eegfiltnew(EEG,icacut,0); % high edge cut-off for stable ICA
            catch me
                error('eegfiltnew not present in eeglab package!')
            end
            
            EEG=pop_chanedit(EEG,'lookup',layout);
            EEGb=pop_chanedit(EEGb,'lookup',layout);
            EEG.chanlocs(nchan+3:end)=[];
            EEG.chanlocs(nchan+1).labels = 'VEOG';
            EEG.chanlocs(nchan+2).labels = 'HEOG';
            
            
            %%
            EEGm = EEG;
            m=[EEGm.event.type];
            if ~isempty(trigger_subtract)
                m=m-(trigger_subtract);
            end
            
            for i=1:length(m)
                EEGm.event(i).type=m(i);
            end
            
            mtable=tabulate(m);
            mtable(mtable(:,2)==0,:)=[];
            mtable=mtable(:,[1 2]);
            disp(mtable)
            disp('Make sure the triggers and their number of ocurrence in the experiment make sense! If so, type in return and hit enter...')
            keyboard
            
            %%
            EEG.event  = EEGm.event;
            EEGb.event = EEGm.event;
                        
            %%
            disp('Saving data of step I')
            save(outfilename1,'EEG')
            save(outfilename1b,'EEGb')
            reload = false;
            
        elseif ~exist(outfilename2,'file')
            
            if reload
                fprintf('Loading subject %i of %i for epoching and bad trial marking...\n',subno,length(sublist));
                load(outfilename1b)
            else
                fprintf('Epoching and bad trial marking for subject %i of %i...\n',subno,length(sublist));
            end
            
            %% Epoch the data
            
            EEGb = pop_epoch( EEGb, [triggers{:,2}],[artiftime]);
            EEGb = pop_rmbase( EEGb, []);
            
            %% check for bad channels
            
            if inspect_chans
                % divide data into data blocks of 10
                figure;
                nblocks = 10;
                ntrials=floor(EEGb.trials/nblocks);
                for b=1:nblocks
                    subplot(3,ceil(nblocks/3),b)
                    newdata = reshape(squeeze(EEGb.data(1:nchan,:,1+((b-1)*ntrials):b*ntrials)),nchan,size(EEGb.data,2)*ntrials);
                    zstd = std(zscore(newdata),[],2);
                    topoplot(zstd,EEGb.chanlocs(1:nchan),'electrodes','on','maplimits',[min(zstd) max(zstd)]);
                    colorbar
                end
                colormap hot
                saveas(gcf,[writdir sublist{subno}(1:4) filesep sublist{subno}(1:4) '_bad_electrode_check.png'])
                figure
                subplot(211)
                topoplot([],EEGb.chanlocs(1:nchan),'style','empty','electrodes','labels');
                subplot(212)
                topoplot([],EEGb.chanlocs(1:nchan),'style','empty','electrodes','numbers');
                keyboard
            end
            
            %% custom-written artifact rejection based on Fieldtrip routines
            
            FT_EEG = eeglab2ft(EEGb);
            
            cfg = [];
            cfg.channel = {EEGb.chanlocs(1:nchan).labels}';
            dat = ft_selectdata(cfg,FT_EEG);
            dat = dat.trial;
            dat_filt = cell(size(dat));
            
            cfg=[];
            % these are standard fieldtrip settings taken from online tutorial
            cfg.bpfreq = [110 140];
            cfg.bpfiltord = 6;
            cfg.bpfilttype = 'but';
            cfg.hilbert = 'yes';
            cfg.boxcar = 0.2;
            cfg.channels = 1:nchan;
            cfg.cutoff = artcutoff;
            cfg.art_time = artiftime; % window within which artifacts are not allowed; note: also use this window for ICA!
            
            %-% loop over trials
            reverseStr='';
            numtrl=EEGb.trials;
            tidx = dsearchn(FT_EEG.time{1}',cfg.art_time')';
            
            for ei=1:numtrl
                
                % display progress
                msg = sprintf('Filtering trial %i/%i...',  ei,numtrl);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                
                % filter in high-freq band
                tmpdat = ft_preproc_bandpassfilter(dat{ei},EEGb.srate, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype);
                tmpdat = ft_preproc_hilbert(tmpdat,cfg.hilbert);
                nsmp = round(cfg.boxcar*EEGb.srate);
                if ~rem(nsmp,2)
                    % the kernel should have an odd number of samples
                    nsmp = nsmp+1;
                end
                tmpdat = ft_preproc_smooth(tmpdat, nsmp); % better edge behaviour
                dat_filt{ei} = double(tmpdat(:,tidx(1):tidx(2)));
                
                if ei==1
                    sumval = zeros(size(dat_filt{1},1), 1);
                    sumsqr = zeros(size(dat_filt{1},1), 1);
                    numsmp = zeros(size(dat_filt{1},1), 1);
                    numsgn = size(dat_filt{1},1);
                end
                
                % accumulate the sum and the sum-of-squares
                sumval = sumval + sum(dat_filt{ei},2);
                sumsqr = sumsqr + sum(dat_filt{ei}.^2,2);
                numsmp = numsmp + size(dat_filt{ei},2);
            end
            fprintf('\n')
            
            % avg and std
            datavg = sumval./numsmp;
            datstd = sqrt(sumsqr./numsmp - (sumval./numsmp).^2);
            
            zmax = cell(1, numtrl);
            zsum = cell(1, numtrl);
            zindx = cell(1, numtrl);
            
            indvec = ones(1,numtrl);
            for ei = 1:numtrl
                % initialize some matrices
                zmax{ei}  = -inf + zeros(1,size(dat_filt{ei},2));
                zsum{ei}  = zeros(1,size(dat_filt{ei},2));
                zindx{ei} = zeros(1,size(dat_filt{ei},2));
                
                nsmp          = size(dat_filt{ei},2);
                zdata         = (dat_filt{ei} - datavg(:,indvec(ei)*ones(1,nsmp)))./datstd(:,indvec(ei)*ones(1,nsmp));  % convert the filtered data to z-values
                zsum{ei}   = nansum(zdata,1);                   % accumulate the z-values over channels
                [zmax{ei},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
                zindx{ei}      = cfg.channels(ind);                % also remember the channel number that has the largest z-value
                
                zsum{ei} = zsum{ei} ./ sqrt(numsgn);
            end
            cfg.threshold = median([zsum{:}]) + abs(min([zsum{:}])-median([zsum{:}])) + cfg.cutoff;
            
            figure
            plot([zsum{:}])
            hold on
            plot([1 length([zsum{:}])], [cfg.threshold cfg.threshold],'m')
            saveas(gcf,[ writdir sublist{subno}(1:4) filesep sublist{subno}(1:4) '_zvals_rejectedtrials.png'])
            rejected_trials = zeros(1,numtrl);
            for ei=1:numtrl
                if sum(zsum{ei}>cfg.threshold)>0
                    rejected_trials(ei) = 1;
                end
            end
            
            rejeegplottmp = trial2eegplot(  rejected_trials, zeros(EEGb.nbchan,EEGb.trials), EEGb.pnts, [1 1 0.5]);
            eegplot(EEGb.data,'srate', EEGb.srate,'limits', [EEGb.xmin EEGb.xmax]*1000 ,'eloc_file', EEGb.chanlocs,'winrej',rejeegplottmp);

            fprintf(['At this point, do following checks:\n'...
                '- check if rejected trials make sense, also check zsum and threshold figure\n' ...
                '- turn back to topoplots of possible bad channels, and verify (add them to .txt file with name as specfied in cfg.chanfilename)\n'...
                '- bad trials and bad channels need to be removed before accurate ICA can be done!\n'...
                '>> DO NOT close the eegplot figure if you want to change the marked trials!\n'...
                '>> Click on epochs to add/remove marked trials; newly added trials are blue, already marked trials are yellow\n'...
                '>> If satisfied, just type "return" in the command window. The figure will close automatically.\n\n'])
            keyboard

            h=get(gcf);
            tmprej = h.UserData.winrej;
            tmprej = eegplot2trial(tmprej,EEGb.pnts,EEGb.trials);
            close 'Scroll activity -- eegplot()'

            if sum(tmprej)~=sum(rejected_trials)
                fprintf('Manually Added/removed %i trials to the automatically detected trials...\n', sum(tmprej)-sum(rejected_trials))
                rejected_trials = tmprej;
            end
            
            %%
            disp('Saving data of step II')
            save(outfilename2,'rejected_trials','cfg'); % the cfg variable contains the rejection settings, so these can be traced back (e.g. z-val cutoff)
            reload = false;
            
        elseif ~exist(outfilename3,'file')
            
            if reload
                fprintf('Loading subject %i of %i for ICA...\n',subno,length(sublist));
                load(outfilename1b)
                load(outfilename2)
                EEGb = pop_epoch( EEGb, [triggers{:,2}],[artiftime]);
                EEGb = pop_rmbase( EEGb, []);
            else
                fprintf('Running ICA for subject %i of %i...\n',subno,length(sublist));
            end
            
            %% Independent Component Analysis
            
            % remove trials with artifacts detected in previous step
            EEGb = pop_select(EEGb,'notrial',find(rejected_trials));
            
            %% sets bad channels to zero if necessary
            
            fid=fopen([writdir chanfilename],'r');
            chans2interp={};
            while ~feof(fid)
                aline=regexp(fgetl(fid),'\t','split');
                chans2interp{size(chans2interp,1)+1,1}=aline{1};
                for i=2:length(aline)
                    chans2interp{size(chans2interp,1),i}=aline{i};
                end
            end
            
            chanind=1:nchan;
            subject_prefix = sublist{subno}(1:4);
            chans = chans2interp(strcmpi(chans2interp(:,1),subject_prefix),2:end);
            if ~isempty(chans{1})
                bad_chansidx = zeros(1,length(chans));
                for c=1:length(chans)
                    if ~isempty(chans{c}), bad_chansidx(c) = find(strcmpi({EEGb.chanlocs.labels},chans(c))); end
                end
                bad_chansidx(bad_chansidx==0)=[];
                chanind(bad_chansidx)=[];
            end
            
            %% run ICA
            
            EEGb=pop_runica(EEGb,'icatype','runica','dataset',1,'options',{'extended',1},'chanind',chanind);
            ICAEEG.icawinv = EEGb.icawinv;
            ICAEEG.icasphere = EEGb.icasphere;
            ICAEEG.icaweights = EEGb.icaweights;
            ICAEEG.icachansind = EEGb.icachansind;
            
            %%
            disp('Saving data of step III')
            save(outfilename3,'ICAEEG')
            reload = false;
            
        elseif ~exist(outfilename4,'file')
            
            if reload
                fprintf('Loading subject %i of %i for final cleaning...\n',subno,length(sublist));
                load(outfilename1)
                load(outfilename1b)
                load(outfilename2)
                load(outfilename3)
                
                EEGb = pop_epoch( EEGb, [triggers{:,2}],[artiftime]);
                EEGb = pop_rmbase( EEGb, []);
                EEGb = pop_select(EEGb,'notrial',find(rejected_trials));
                EEGb.icachansind = ICAEEG.icachansind;
                EEGb.icasphere = ICAEEG.icasphere;
                EEGb.icaweights = ICAEEG.icaweights;
                EEGb.icawinv = ICAEEG.icawinv;
                EEGb = eeg_checkset(EEGb);

            else
                fprintf('Final cleaning steps for subject %i of %i...\n',subno,length(sublist));
            end
                                    
            %% inspect components
            
            if inspect_ica
                pop_selectcomps(EEGb,[1:20]); % plot first 20 ICs
                pop_eegplot( EEGb, 0, 1, 1);  % plot the component activation
                pop_eegplot( EEGb, 1, 0, 0);
                
                inspecting = true;
                while inspecting
                    prompt = 'Which component do you want to inspect? (type number or "done" when done) --> ';
                    comp2check = input(prompt,'s');
                    if strcmp(comp2check,'done')
                        inspecting = false;
                    else
                        pop_prop( EEGb, 0, str2double(comp2check));
                    end
                end
            end
            clear EEGb
            
            %% Mark bad channels and ICs to remove
            
            fid=fopen([writdir chanfilename],'r');
            chans2interp={};
            while ~feof(fid)
                aline=regexp(fgetl(fid),'\t','split');
                chans2interp{size(chans2interp,1)+1,1}=aline{1};
                for i=2:length(aline)
                    chans2interp{size(chans2interp,1),i}=aline{i};
                end
            end
            
            chanind=1:nchan;
            subject_prefix = sublist{subno}(1:4);
            chans = chans2interp(strcmpi(chans2interp(:,1),subject_prefix),2:end);
            bad_chansidx = [];
            if ~isempty(chans{1})
                bad_chansidx = zeros(1,length(chans));
                for c=1:length(chans)
                    if ~isempty(chans{c}), bad_chansidx(c) = find(strcmpi({EEG.chanlocs.labels},chans(c))); end
                end
                bad_chansidx(bad_chansidx==0)=[];
                chanind(bad_chansidx)=[];
            end
            
            fid=fopen([writdir icafilename],'r');
            comps2remove={};
            while ~feof(fid)
                aline=regexp(fgetl(fid),'\t','split');
                comps2remove{size(comps2remove,1)+1,1}=aline{1};
                comps2remove{size(comps2remove,1)  ,2}=sscanf(aline{2},'%g');
            end
            
            %%
            try
                disp(['Components removed for subject ' num2str(subno) ': ' num2str(comps2remove{subno,2}')])
            catch me
                disp('No components specified for this subject yet! Turn to your .txt file where you mark the components!')
            end
            
            %% do the cleaning:
            %- epoch in wider window
            %- baseline correct
            %- reject trials
            %- add ICA weights
            %- remove EOG
            %- interpolate channels
            %- remove components
            
            EEG = pop_epoch( EEG, [triggers{:,2}],epochtime);
            EEG = pop_rmbase( EEG, []);
            for ei=1:EEG.trials
                EEG.epoch(ei).trialnum = ei;
            end
            EEG = pop_select(EEG,'notrial',find(rejected_trials));

            EEG.icachansind = ICAEEG.icachansind;
            EEG.icasphere = ICAEEG.icasphere;
            EEG.icaweights = ICAEEG.icaweights;
            EEG.icawinv = ICAEEG.icawinv;
            EEG = eeg_checkset(EEG);
            clear ICAEEG
            
            % remove EOG for good
            EEG = pop_select(EEG,'nochannel',[nchan+1 nchan+2]);
            
            if sum(bad_chansidx)>0
                
                EEG2 = pop_select(EEG,'nochannel',bad_chansidx);
                EEG2 = pop_subcomp( EEG2, comps2remove{subno,2}', 0);
                
                %% put IC-cleaned channels back in original EEG structure and interpolate bad ones
                
                good_chans = 1:EEG.nbchan;
                good_chans(bad_chansidx)=[];
                EEG.data(good_chans,:,:) = EEG2.data;
                EEG = eeg_interp(EEG,bad_chansidx);
                
                clear EEG2
            else
                EEG = pop_subcomp( EEG, comps2remove{subno,2}', 0);
            end
            clear EEGb
            
            %%
            disp('Saving data of step IV, done preprocessing for this subject!')
            save(outfilename4,'EEG','comps2remove','bad_chansidx');
            go_on = false; % this gets you out of the while loop, so the subject loop continues
        else
            go_on = false; % this gets you out of the while loop, so the subject loop continues
        end
    end
end
end

