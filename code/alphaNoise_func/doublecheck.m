function doublecheck(EEG)


%% check EEG triggers

% stimulus
trigs_strings = {'1', '2', '3'};
cell_events = {EEG.event(:).type};

mask_trigs_stim = false(length(cell_events), 1);

for ievent = trigs_strings

    temp_ = ismember(cell_events, ievent)';
    mask_trigs_stim = mask_trigs_stim | temp_;

end

stim_on = cell_events(mask_trigs_stim);
stim_on = cellfun(@(x) str2double(x), stim_on);

%% check final consistency with beh
stim_on_beh = cell2mat({EEG.trialinfo.signal}) + 1;

out_ = stim_on == stim_on_beh;
perfect_correspondence = all(out_);

if perfect_correspondence

    disp('SUCCESS!!!')
    disp('the data passed the final sanity check for coregistration!!!')

else

    error('PORCODIO NO CORRESPONDENCE')


end






