This Matlab code is meant for preprocessing EEG data, and tested on 64 channel Biosemi data.
With some custom adjustments, it may well be suited for other electrophysiological measurement systems as well.

It needs eeglab and fieldtrip. Newest versions of those packages are not necessary, except that eeglab needs eegfiltnew.m.

The code run_eegpreproc.m specified some parameters and path settings, which is then passed on to the function eegpreproc.m. This function does four steps:
1) Re-referencing and high-pass filtering
2) Epoching and trail/channel rejection marking
3) ICA
4) Final cleaning

It does these steps consecutively, saving .mat files in between, and looping over subjects. When stopped and reran, it continues where it left off.

In between steps, the function requires some manual input. It’s interactive. This is because preprocessing should not be done fully automatically; some visual inspection is required to see if all goes well, and to make informed decisions. For example, in this pipeline, which bad channels to interpolate and which independent components to remove from the data is done by visual inspection. However, this is done while running the function; it doesn’t need to be stopped.

Please note that there is no unique, correct way of preprocessing EEG data. The steps presented are a result of roughly 7 years of (personal and collaborative) experience.

The pdf manual that is currently in this repo is a bit old. I’m working on a newer manual.

Please contact joramvandriel@gmail.com for questions and remarks.