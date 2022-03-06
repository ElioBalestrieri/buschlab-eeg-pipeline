function EEG = sanity_check_coregister(EEG)

% store as first thing the original tinfo
EEG.full_trialinfo = EEG.trialinfo;


% collect responses in a vector
urev_cell = {EEG.urevent(:).type}; 
% chop char event
all_chars = cellfun(@(x) ischar(x), urev_cell);
urev_cell(all_chars) = [];
urev_array = cell2mat(urev_cell);
trigon = ismember(urev_array, [1, 2, 3]);
responses_array = urev_array(find(trigon)+1);


% simply reject from trialinfo the trials erased so far
if isfield(EEG, 'rejected_trials')

    EEG.trialinfo(EEG.rejected_trials) = [];
        
end

if isfield(EEG, 'rejected_trials_afterICA')

    EEG.trialinfo(EEG.rejected_trials_afterICA) = [];
        
end


% collect the events of interest
events_cell = cellfun(@(x) str2double(x(1)), {EEG.event(:).type}, 'UniformOutput', false);
maskempty = cellfun(@(x) isempty(x), events_cell);
eventriggers = cell2mat(events_cell(~maskempty));
eventriggers = eventriggers(~isnan(eventriggers));

beh_sig = cellfun(@(x) x, {EEG.trialinfo(:).signal}) + 1; % + 1, and "signal" specific for alphaNoise
beh_resp = cellfun(@(x) x, {EEG.trialinfo(:).bin_resp})*10 + 10;

try

    if all(eventriggers == beh_sig)
        
        disp('SUCCESS')
        disp('Perfect correspondence between behavioural data and triggercdes')
        
    else
           
        warning('sorry, no correspondence between trigger codes and behavioural data')
        fprintf('\n\nLack of correspondence in %i out of %i\n\n', ...
            sum(eventriggers == beh_sig), n_Trials)
        
        error('cant go on')
        
    end

catch

    if strcmp(EEG.urname(1:2), '13')

        EEG = myfixer(EEG, eventriggers, responses_array, ...
            beh_sig, beh_resp, EEG.urname);

    else

        EEG = myfixer(EEG, eventriggers, responses_array, beh_sig, beh_resp);

    end
        
    n_discEEG = length(EEG.final_EEG_rejections);
    n_discBEH = length(EEG.final_beh_rejections);

    fprintf('fixing procedure ended:\n%i trials were further discarded from EEG\n%i trials were further discarded from behaviour', n_discEEG, n_discBEH);
    if (n_discEEG>50) || (n_discBEH>50)

        error('too many trials discarded')

    end
    
end

end



function EEG = myfixer(EEG, eventriggers, responses_array, beh_sig, beh_resp, varargin)

    fprintf('\nDetected mismatch in coregistering EEG and behavioural data, attempting to fix this...\n')
    
    itisnotfixed = true;

    swap_beh = beh_sig;
    swap_eeg = eventriggers;

    cellcorrections = {};


    if ~isempty(varargin)

        if strcmp(varargin{1}, '13RZ.bdf')

            % weird triggers on 13
            badapples = [714, 742, 835];
            badapples_EEG = [badapples, 836];
        
            % eraser   
            beh_resp(badapples) = [];
            swap_beh(badapples) = [];
            swap_eeg(badapples_EEG) = [];
            responses_array(badapples_EEG) = [];

            cellcorrections = {[714, 714], [741, 741], [833, 833], [nan, 833]};

        end

    end


    celltrigs = {swap_beh, swap_eeg};
    cellresps = {beh_resp, responses_array};
    
    cycles_on = 0;
    while itisnotfixed
        
        cycles_on = cycles_on +1;

        L_arrays = cellfun("length", celltrigs);

        [map_trig, map_resp] = deal(false(min(L_arrays),1));

        for itrial = 1:min(L_arrays)

            map_trig(itrial) = celltrigs{1}(itrial) == celltrigs{2}(itrial);
            map_resp(itrial) = cellresps{1}(itrial) == cellresps{2}(itrial);

        end

        map_tot = map_trig & map_resp;
        problem_idx = find(map_tot==0, 1, 'first');

        if isempty(problem_idx)

            break

        else

            vectcorr = nan(1, 2);

            if L_arrays(1) ~= L_arrays(2)
            
                whichlonger = L_arrays==max(L_arrays);
                celltrigs{whichlonger}(problem_idx) = [];
                cellresps{whichlonger}(problem_idx) = [];
                vectcorr(whichlonger) = problem_idx;

            else
                
                for ientry=1:length(celltrigs)  
                    celltrigs{ientry}(problem_idx) = [];
                    cellresps{ientry}(problem_idx) = [];
                    vectcorr(ientry) = problem_idx;
                end

            end
            
            cellcorrections{cycles_on} = vectcorr;

        end


    end

    EEG.cellcorrections = cellcorrections;

    EEG = local_apply_corrections(EEG);


end



function EEG = local_apply_corrections(EEG)

final_EEG_rejections = [];
final_beh_rejections = [];

for icorrections = 1:length(EEG.cellcorrections)

    this_vect = EEG.cellcorrections{icorrections};

    if ~isnan(this_vect(1))
        final_beh_rejections(icorrections) = this_vect(1);
        EEG.trialinfo(this_vect(1)) = [];
    end

    if ~isnan(this_vect(2))
        final_EEG_rejections(icorrections) = this_vect(2);
        vect_false = false(size(EEG.data, 3), 1);
        vect_false(this_vect(2)) = true;
        EEG = pop_rejepoch(EEG, vect_false, 0);
    end


end

EEG.final_beh_rejections = final_beh_rejections;
EEG.final_EEG_rejections = final_EEG_rejections;


end




%% sandbox

% this chunk was coming from before the badgaze trials were erased in the import
% script

% if length(responses_array) > length(EEG.trialinfo)
%     
%     % it means that a bad trig (or more than one) was passed into the data.
%     % this happens when people break fixation after stimulus presentation,
%     % for instance.
%     % remove the instances generated in this way and correct the trial to
%     % be rejected
% 
%     fulldata_wrongtrigs = (responses_array ~= 10) & (responses_array ~= 20);
%     whoiswrong = find(fulldata_wrongtrigs);
% 
%     if (length(responses_array)-length(whoiswrong) ~= length(EEG.trialinfo)) 
% 
%         % hence, exception reached in the definition of bad trials
%         whoswrongright = responses_array(whoiswrong)~=66;
%         whoiswrong(whoswrongright) = [];
% 
%     end
% 
%     if isfield(EEG, 'rejected_trials')
% 
%         fulldata_wrongtrigs(EEG.rejected_trials) = [];
%         responses_array(EEG.rejected_trials) = [];
% 
%         for wrongtrigs = whoiswrong
%             lgcl_abovewrong = EEG.rejected_trials > wrongtrigs;
%             EEG.rejected_trials = EEG.rejected_trials - lgcl_abovewrong;
%         end
%     end
% 
%     if isfield(EEG, 'rejected_trials_afterICA')
% 
%         fulldata_wrongtrigs(EEG.rejected_trials_afterICA) = [];
%         responses_array(EEG.rejected_trials_afterICA) = [];
% 
%         for wrongtrigs = whoiswrong
%             lgcl_abovewrong = EEG.rejected_trials_afterICA > wrongtrigs;
%             EEG.rejected_trials_afterICA = EEG.rejected_trials_afterICA - lgcl_abovewrong;
%         end
%     end
% 
%     % finally get rid of this last fucking bastard
%     EEG = pop_rejepoch(EEG, fulldata_wrongtrigs, 0);
%     responses_array(fulldata_wrongtrigs) = [];
% 
% end




