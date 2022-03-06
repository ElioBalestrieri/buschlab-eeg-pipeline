%% Set preferences, configuration and load list of subjects.
clear; clc; close all

restoredefaultpath
prefs = get_prefs_alphaNoise('eeglab_all', 1);
cfg   = get_cfg_alphaNoise;

% ------------------------------------------------------------------------
% **Important**: these variables determine which data files are used as
% input and output. 
suffix_in  = '';
suffix_out = 'import';
do_overwrite = true;
% ------------------------------------------------------------------------

subjects = get_list_of_subjects(cfg.dir, do_overwrite, suffix_in, suffix_out);
% do not consider among subjects the doubel datafile for subj 15
subjects(16) = [];


%% Run across subjects.
nthreads = min([prefs.max_threads, length(subjects)]);
nthreads = nthreads - 5; % buschd30 got stuck with 15 subjects given the huge filesize, so reduce the maxtrhead. Needed only here because further on the data is not that big.
parfor(isub = 1:length(subjects), nthreads) % set nthreads to 0 for normal for loop.
% for isub = 1%:length(subjects)
    
    if ~exist(subjects(isub).outdir, 'dir')
        mkdir(subjects(isub).outdir)
    end
    
    % --------------------------------------------------------------
    % Import Biosemi raw data.
    % --------------------------------------------------------------
    if strcmp(subjects(isub).name(1:4), '15TF') % recording stopped due to alarm

        EEG1 = func_import_readbdf(cfg.dir, subjects(isub).name);
        EEG2 = func_import_readbdf(cfg.dir, '15TF_2.bdf');
        EEG = pop_mergeset(EEG1, EEG2);

        swap = EEG;

        for ievent = 1:length(EEG.event)

            swap.event(ievent).type = str2double(EEG.event(ievent).type);
            itemurevent = EEG.urevent(ievent).type;
            
            if isnumeric(itemurevent)

                swap.urevent(ievent).type = itemurevent;

            else

                swap.urevent(ievent).type = nan;

            end
                

        end

        test_events = cell2mat({swap.event(:).type});
        test_urevents = cell2mat({swap.urevent(:).type});

        unique(test_events)
        unique(test_urevents)

        swap.event(isnan(test_events)) = [];
        swap.urevent(isnan(test_urevents)) = [];

        EEG = swap;

    else
        EEG = func_import_readbdf(cfg.dir, subjects(isub).name);
    end
    
    % ---------------------------------------------------------------------
    % add exceptions: 
    % ---------------------------------------------------------------------
    if strcmp(subjects(isub).name(1:4), '02LT')

        % channels B14 and B16 were swapped during the
        % recording. Correct this immediately
        B14_chan = EEG.data(46, :);
        B16_chan = EEG.data(48, :);
        EEG.data(46, :) = B16_chan;
        EEG.data(48, :) = B14_chan;  
        
    elseif strcmp(subjects(isub).name(1:4), '06KD')
        
        % A28 substituted with ext2 for noise
        new_A28_chan = EEG.data(66, :);
        EEG.data(28, :) = new_A28_chan;

    elseif strcmp(subjects(isub).name(1:4), '20NB')
        
        % in this subject there was a mysterious glitch with trigger pins
        % the 8 pin was always active (???)
        % apply correction here

        vect_trigs = cell2mat({EEG.event(:).type})';
        trl_onsets = vect_trigs == 74;
        trgt_onsets_idxs = find(trl_onsets) +1;
        trgt_trigs = vect_trigs(trgt_onsets_idxs);
        
        vect_trigs(trgt_onsets_idxs) = trgt_trigs - 8;
        vect_trigs(trl_onsets) = 66;
        vect_trigs(vect_trigs == 28) = 20;

        % after this there was an additional "8" somewhere. Make it NaN
        vect_trigs(vect_trigs==8) = nan;

        for itrig = 1:length(vect_trigs)
            EEG.event(itrig).type = vect_trigs(itrig);
            EEG.urevent(itrig).type = vect_trigs(itrig);
        end

        whosnan = find(isnan(vect_trigs)); 
        EEG.event(whosnan) = [];
        EEG.urevent(whosnan) = [];

    end


    % --------------------------------------------------------------
    % Select data channels.
    % --------------------------------------------------------------
    EEG = func_import_selectchans(EEG, cfg.chans);
    
    % --------------------------------------------------------------
    % Biosemi is recorded reference-free. We apply rereferencing in
    % software.
    % --------------------------------------------------------------
    EEG = func_import_reref(EEG, cfg.prep);
    
    % --------------------------------------------------------------
    % Compute VEOG and HEOG.
    % --------------------------------------------------------------
    EEG = func_import_eyechans(EEG, cfg.chans);
    
    % --------------------------------------------------------------
    % Filter the data.
    % --------------------------------------------------------------
    % We want to keep the VEOG/HEOG data unfiltered to make sure they
    % are not distorted by the filter. We keep a copy here and then put
    % it back after filtering.
    tmp = EEG.data;
    EEG = func_import_filter(EEG, cfg.prep);
    EEG.data(cfg.chans.VEOGchan,:) = tmp(cfg.chans.VEOGchan,:);
    EEG.data(cfg.chans.HEOGchan,:) = tmp(cfg.chans.HEOGchan,:);
    
    %---------------------------------------------------------------
    % Remove all events from non-configured trigger devices
    %---------------------------------------------------------------
    EEG = func_import_remove_triggers(EEG, cfg.epoch);
    
    % --------------------------------------------------------------
    % Import Eyetracking data.
    % --------------------------------------------------------------
    eyetrack = cfg.eyetrack;
    eyetrack.eye_startEnd = [];
    eyetrack.eye_startEnd = [EEG.event(1).type, EEG.event(end).type];   

    EEG = func_import_importEye(EEG, subjects(isub).namestr, cfg.dir, eyetrack);    
    
    % --------------------------------------------------------------
    % Downsample data if required. IMPORTANT: use resampling only after
    % importing the eye tracking data, or else the ET data will not be in
    % sync with EEG data.
    % --------------------------------------------------------------
    EEG = func_import_downsample(EEG, cfg.prep);
               
    % -------------------------------------------------------------
    % Epoch data.
    % --------------------------------------------------------------
    EEG = func_import_epoch(EEG, cfg.epoch, cfg.eyetrack.coregister_Eyelink);
    %     EEG = pop_rmbase(EEG, [], []);
        
    % --------------------------------------------------------------
    % Import behavioral data.
    % --------------------------------------------------------------
%    EEG = func_importBehavior(EEG, subjects(isub).namestr, cfg.dir, cfg.epoch);
    [EEG, Trial] = func_read_csv_logfile(cfg, EEG, subjects(isub).name(1:4));

    % adjust trials in excess
    EEG = first_coregistercheck(EEG);

    % --------------------------------------------------------------
    % Save the new EEG file in EEGLAB format.
    % --------------------------------------------------------------
    EEG = func_saveset(EEG, subjects(isub));
    
end

disp('Done.')
