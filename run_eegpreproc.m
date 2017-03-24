%% run brand new eegpreproc function

restoredefaultpath

%%
cd('Z:\Stuff\eegpreproc\')

%-% part I: loading raw data and filter
cfg = [];
cfg.eeglab_path = 'Z:\Toolboxes\eeglab12_0_2_3b';
cfg.ft_path     = 'Z:\Toolboxes\fieldtrip-20150318';
cfg.readdir     = 'Z:\TempOrient\EEG\raw\';
cfg.writdir     = 'Z:\Stuff\eegpreproc\';
cfg.layout      = [cfg.eeglab_path '\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp'];
cfg.projectname = 'temporient';

cfg.nchan       = 64;
cfg.reref       = {'EXG5','EXG6'};
cfg.veog        = {'EXG1','EXG2'};
cfg.heog        = {'EXG3','EXG4'};
cfg.resrate     = 512;

%-% part II: epoching and trial rejection marking

triggers={
	'conditionA'	{ '1' };
	'conditionB'	{ '2' };
	};

cfg.epochtime   = [-1.5 2.5];
cfg.artiftime   = [-0.5 1];
cfg.artcutoff   = 12;
cfg.triggers    = triggers;
cfg.trigger_subtract = []; % depending on physical lab settings, sometimes weird high numbers get added to the trigger values specified in your experiment script

%-% part III: ICA
cfg.chanfilename = 'chans2interp.txt';

%-% part IV: final cleaning
cfg.icafilename = 'ICs2remove.txt';
cfg.inspect_ica = true

%-% now run it
eegpreproc(cfg);
