%% Set preferences, configuration and load list of subjects.
clear; clc; close all

restoredefaultpath
prefs = get_prefs_alphaNoise('eeglab_all', 0); 
cfg   = get_cfg_alphaNoise;

% ------------------------------------------------------------------------
% **Important**: these variables determine which data files are used as
% input and output. 
suffix_in  = 'final';
suffix_out = ''; 
do_overwrite = true;
% ------------------------------------------------------------------------

subjects = get_list_of_subjects(cfg.dir, do_overwrite, suffix_in, suffix_out);

%% Run across subjects.
nthreads = min([prefs.max_threads, length(subjects)]);
% parfor(isub = 1:length(subjects), nthreads)
for isub = 1:length(subjects)
    
    % --------------------------------------------------------------
    % Load the dataset 
    % --------------------------------------------------------------
    EEG = pop_loadset('filename', subjects(isub).name, 'filepath', subjects(isub).folder);    
     
    % create the dataset(s) we need
    split_datasets_simple(cfg, EEG)

    fprintf('\nFinished with subject %s\n', subjects(isub).name(1:4))
    nfinal_trls = size(EEG.data, 3);
    fprintf('\nKept %i/1008 trials\n', nfinal_trls)
    perc_trls = (1008-nfinal_trls)/1008;
    fprintf('\nThis corresponds to a proportion of %1.2f trials excluded\n', perc_trls)
    
    if perc_trls > .15

        warning('.15 threshold exceeded')

    end

    proportion_rejected(isub) = perc_trls;
    
end
disp('Done.')
