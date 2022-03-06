function EEG = first_coregistercheck(EEG)


cell_typevent = {EEG.urevent(:).type};

acc_idxs = 0; this_mat_evs = nan(1500, 4); 
for iev = 1:length(cell_typevent)

    this_ev = cell_typevent{iev};

    if isnumeric(this_ev)

        if (this_ev==1) || (this_ev==2) || (this_ev==3)
    
            acc_idxs = acc_idxs +1;
    
            this_mat_evs(acc_idxs, 1) = this_ev;
            try
    
                this_mat_evs(acc_idxs, 2) = cell_typevent{iev+1};
                this_mat_evs(acc_idxs, 3) = cell_typevent{iev+2};
                this_mat_evs(acc_idxs, 4) = cell_typevent{iev+3};
     
            catch
                disp('anomaly')
            end

        end
    
    end

end
this_mat_evs(acc_idxs+1:end, :) = [];

%% prune entries with failed fixcontrol
sus2 = this_mat_evs(:, 2)==66;
pruned_mat = this_mat_evs(~sus2,:);

%% check congruency

stimtype = cell2mat({EEG.trialinfo(:).signal})+1;
congruent_stimtype = pruned_mat(:, 1) == stimtype';

if sum(congruent_stimtype)==1008

    badidxs = find(sus2);
    EEG = pop_select( EEG, 'notrial', badidxs);
    disp('Success!')
    
else

    warning('no coreg compatibility ')

end


end