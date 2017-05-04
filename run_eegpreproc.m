%% run brand new eegpreproc function

restoredefaultpath

%%
addpath('Z:\Toolboxes\Git\eegpreproc\')

%-% part I: loading raw data and filter
cfg = [];
cfg.eeglab_path = 'path/to/eeglab12_0_2_3b';
cfg.ft_path     = 'path/to/fieldtrip-20150318';
cfg.readdir     = 'some_dir_to/EEG/raw';
cfg.writdir     = 'some_dir_to/EEG/processed';
cfg.layout      = [cfg.eeglab_path '\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp'];
cfg.projectname = 'a_sensible_name';

cfg.nchan       = 64;
cfg.reref       = {'EXG5','EXG6'};
cfg.veog        = {'EXG1','EXG2'};
cfg.heog        = {'EXG3','EXG4'};
cfg.resrate     = 512;

%-% part II: epoching and trial rejection marking

triggers={
    'conditionA'        { '1' };
    'conditionB'        { '2' };
    };

cfg.epochtime   = [-1.5 7.5];
cfg.artiftime   = [-0.25 5.25];
cfg.artcutoff   = 12;
cfg.triggers    = triggers;
cfg.trigger_subtract = []; % depending on physical lab settings, sometimes weird high numbers get added to the trigger values specified in your experiment script
cfg.inspect_chans = true; % if true, pauses function and shows figures with topomaps of variance, to highlight bad channels

%-% part III: ICA
cfg.chanfilename = 'chans2interp.txt';

%-% part IV: final cleaning
cfg.icafilename = 'ICs2remove.txt';
cfg.inspect_ica = true;

%-% now run it
eegpreproc(cfg);
