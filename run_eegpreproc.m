%% run brand new eegpreproc function

restoredefaultpath
cd('Z:\Stuff\eegpreproc\')

%% part I
cfg = [];
cfg.eeglab_path = 'Z:\Toolboxes\eeglab12_0_2_3b';
cfg.ft_path     = 'Z:\Toolboxes\fieldtrip-20150318';
cfg.readdir     = 'Z:\TempOrient\EEG\raw\';
cfg.writdir     = 'Z:\Stuff\eegpreproc\';
cfg.layout      = [cfg.eeglab_path '\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp'];
cfg.projectname = 'temporient';

cfg.epochtime   = [-1.5 2.5];
cfg.nchan       = 64;
cfg.reref       = {'EXG5','EXG6'};
cfg.veog        = {'EXG1','EXG2'};
cfg.heog        = {'EXG3','EXG4'};
cfg.resample    = 512;

cfg.triggers    = [21:23 31:33 41:43];
cfg.trigger_subtract = [];

eegpreproc(cfg);

%% part II



