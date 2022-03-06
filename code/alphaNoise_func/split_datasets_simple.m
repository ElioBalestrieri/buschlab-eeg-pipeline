function split_datasets_simple(cfg, EEG)
% alphaNoise specific
% the function splits the dataset into:
% 1. A post-stimulus only, python easy readable data format to perform decoding analysis
% 2. Full time window simplified mat version with energy maps 

% % add foof folder
% addpath('/data5/Elio/fooof_mat/fooof_mat')

%% post stim data
% internal cfg for poststim def
base_tw = [-.2 0];
new_tw = [-.2, 1.5] * 1e3; % times in ms
channels = 1:64;

% parse subject name and define output folder
subjname = EEG.urname(1:4);
outfolder = fullfile(cfg.dir.main, 'Preprocessing', 'preprocessed', ...
                     'decode_postim', subjname);

if ~exist(outfolder, 'dir')
    mkdir(outfolder)
end
                 
% baseline correction and time window/channel selection
baseline_corrected = pop_rmbase(EEG, base_tw);
mask_time = ((baseline_corrected.times >= min(new_tw)) &...
             (baseline_corrected.times <= max(new_tw)));
mat_data = baseline_corrected.data(channels, mask_time, :);
time_vect = baseline_corrected.times(mask_time);
logfile_table = struct2table(baseline_corrected.trialinfo);

% save files
save(fullfile(outfolder, [subjname, '_mat_data.mat']), 'mat_data');
save(fullfile(outfolder, 'time_vect.mat'), 'time_vect');
writetable(logfile_table, fullfile(outfolder, [subjname, '_logfile.csv']))


%% full data

SUBJ = [];

% general
SUBJ.name = subjname;
SUBJ.beh = EEG.trialinfo;
SUBJ.chanlocs = EEG.chanlocs(channels);
SUBJ.srate = EEG.srate;

% poststim window (ERP)
SUBJ.poststim = [];
SUBJ.poststim.baselinecorrected = 'yes';
SUBJ.poststim.data = mat_data;
SUBJ.poststim.time = time_vect;

% full window
SUBJ.fullwin = [];
SUBJ.fullwin.data = EEG.data(channels, :, :);
SUBJ.fullwin.time = EEG.times;

% energy maps
SUBJ.energymaps = [];
path_to_maps = fullfile(cfg.dir.main, 'Preprocessing', 'preprocessed', 'EnergyMaps', subjname, [subjname '.mat']);
load(path_to_maps, 'energymap')

if isfield(EEG, 'rejected_trials')
    energymap(:, :, EEG.rejected_trials) = [];
end

if isfield(EEG, 'rejected_trials_afterICA')
    energymap(:, :, EEG.rejected_trials_afterICA) = [];
end

if isfield(EEG, 'final_beh_rejections')

    for ireject = EEG.final_beh_rejections
        energymap(:, :, ireject) = [];
    end

end

SUBJ.energymaps.data = energymap;
SUBJ.energymaps.dimord = 'freqXangle';
SUBJ.energymaps.freqspace = [1, 3];
SUBJ.energymaps.anglespace = [-pi/6, pi/6];


path_to_finalfile = fullfile(cfg.dir.main, 'Preprocessing', 'preprocessed', 'fulldata_final');

if ~isfolder(path_to_finalfile)
    mkdir(path_to_finalfile);
end
    
save(fullfile(path_to_finalfile, [subjname '.mat']), 'SUBJ');

end

%% ############# LOCAL FUNCTION

function [alphaComps, pow_spctr, avg_spctr] = local_selectAlphaComp(EEG, prestim)

% convert prestim in ms
prestim = prestim*1e3;
alpharange = [8, 13]; 
expfit_freqrange = [3, 30];
cumsum_limit = .95;

%% 1. subselect time window of interest
lgcl_time = EEG.times>=min(prestim) &  EEG.times<max(prestim);
comp_prestim = EEG.icaact(:, lgcl_time, :);

%% component selection based on alpha band
% concatenate along trials
splitTRLS = num2cell(comp_prestim, [1, 2]);
comps_cat = horzcat(splitTRLS{:});

% normalize
zscored_comps = zscore(comps_cat, [], 2);

% reshape
recat_trls = reshape(zscored_comps, size(comp_prestim));

% define signal length
lsign = sum(lgcl_time); length_fft = 2^nextpow2(lsign); % for speed, if necessary

% apply fft
cmplx = fft(recat_trls, length_fft, 2);

% compute power
pow_spctr = (abs(cmplx)/length_fft).^2;    
pow_spctr = pow_spctr(:, 1:(ceil(length_fft/2)+1), :);
fin_signal_length = size(pow_spctr, 2);

xfreq = linspace(0, EEG.srate/2, fin_signal_length);
avg_spctr = mean(pow_spctr, 3);

%% exponential decay function fit
out = local_fit_expmodel(xfreq, avg_spctr, expfit_freqrange);

%% sort components based on local peaks

ncomps = size(avg_spctr, 1);
peakmagnitude = zeros(ncomps, 1);

for icomp = 1:ncomps
    
    [pks, locs] = findpeaks(out.detrended(:, icomp), out.pruned_xfreqs, ...
                            'MinPeakDistance',2);

    lgcl_mask_peakpresent = locs>=min(alpharange) & locs<=max(alpharange);
    
    if any(lgcl_mask_peakpresent)
        
        peakmagnitude(icomp) = max(pks(lgcl_mask_peakpresent));
        
    end
    
end

%% get final component scores
peakmagnitude(peakmagnitude<0) = 0;
[sorted_peaks ,p2] = sort(peakmagnitude,'descend');

norm_cum = cumsum(sorted_peaks) ./ sum(sorted_peaks);

mask_keep_comps = norm_cum <= cumsum_limit;

alphaComps.idxs = p2(mask_keep_comps);
alphaComps.Zpeaks = sorted_peaks(mask_keep_comps);
alphaComps.alpharange = alpharange;
alphaComps.expFreqRange = expfit_freqrange;
alphaComps.prestim = prestim;
alphaComps.xfreq = xfreq;

end

function out = local_fit_expmodel(xfreq, avg_spctr, freqrange)


lgcl_freqmask = (xfreq>=min(freqrange)) & (xfreq<=max(freqrange));

Xfit = xfreq(lgcl_freqmask)';
Yfit = double(avg_spctr(:, lgcl_freqmask)');

anon_expfunc = @(v, x) v(1).*exp(-x./v(2)) + v(3);


ncomps = size(Yfit, 2);

[mat_detrended_spectra, fitted_curves] = deal(nan(sum(lgcl_freqmask), ncomps));
vect_adjR2 = nan(ncomps, 1);

for icomp = 1:ncomps

    vectY = Yfit(:, icomp);
    fo = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower',[0,0,0],...
                    'Upper',[Inf, Inf, Inf],...
                    'StartPoint',[vectY(1), 1, vectY(end)]);
    
    ft = fittype('a.*exp(-x./b) + c','options',fo);
    
    [this_expcurve, gof] = fit(Xfit,vectY,ft);

    rep_curve = anon_expfunc(coeffvalues(this_expcurve), Xfit);

    mat_detrended_spectra(:, icomp) = vectY - rep_curve;
    fitted_curves(:, icomp) = rep_curve;

    vect_adjR2(icomp) = gof.adjrsquare;


end

out.adjR2 = vect_adjR2;
out.detrended = mat_detrended_spectra;
out.fitted_exp = fitted_curves;
out.pruned_xfreqs = Xfit;

end


% function ICAfield = local_selectAlphaComp_old(EEG, prestim)
% 
% % the algorithm follows 2 criterion for selecting the best component:
% % 1. highest power in the alpha band (8-13) AND
% % 2. "bumpiness", aka local maxima, identified with the lowest absolute
% % derivative
% % 3. ? iclabel ?
% % value of the derivative in the amplitude.
% % All of this is referred to the prestimulus window defined in prestim.
% 
% % convert prestim in ms
% prestim = prestim*1e3;
% alpharange = [8, 13]; 
% 
% % 1. subselect time window of interest
% lgcl_time = EEG.times>=min(prestim) &  EEG.times<max(prestim);
% comp_prestim = EEG.icaact(:, lgcl_time, :);
% 
% % compute power spectra, frequencies and average across trials
% % demean
% demeaned_signal = comp_prestim - mean(comp_prestim, 2);
% 
% % define signal length
% lsign = sum(lgcl_time); length_fft = 2^nextpow2(lsign); % for speed, if necessary
% 
% % apply fft
% cmplx = fft(demeaned_signal, length_fft, 2);
% 
% % compute power
% pow_spctr = (abs(cmplx)/length_fft).^2;    
% pow_spctr = pow_spctr(:, 1:(ceil(length_fft/2)+1), :);
% fin_signal_length = size(pow_spctr, 2);
% 
% xfreq = linspace(0, EEG.srate/2, fin_signal_length);
% 
% avg_spctr = mean(pow_spctr, 3);
% 
% % start our ranking algorithm...
% % 1. highest power in the alpha range
% lgcl_alpha = xfreq>=min(alpharange) & xfreq<=max(alpharange);
% alpha_powspctr = mean(pow_spctr(:, lgcl_alpha),2);
% 
% % 2. local peak in alpha range?
% ncomps = size(avg_spctr, 1);
% 
% peakmagnitude= zeros(ncomps, 1);
% 
% for icomp = 1:ncomps
%     
%     [pks, locs] = findpeaks(avg_spctr(icomp, :), xfreq);
% 
%     lgcl_mask_peakpresent = locs>=min(alpharange) & locs<=max(alpharange);
%     
%     if any(lgcl_mask_peakpresent)
%         
%         peakmagnitude(icomp) = max(pks(lgcl_mask_peakpresent));
%         
%     end
%     
% end
% 
% %% get final component scores
% [~,p1] = sort(alpha_powspctr,'descend');
% [~,p2] = sort(peakmagnitude,'descend');
% 
% 
% % this is not working since the list of component probability is still
% % referred to all components before rejection...
% % brainpos = ismember(EEG.etc.ic_classification.ICLabel.classes, 'Brain');
% % [~,p3] = sort(EEG.etc.ic_classification.ICLabel.classifications(:, brainpos), ...
% %               'descend'); % sort the probability of being brain
% 
% 
% [r1, r2] = deal((1:ncomps)');
% 
% r1(p1) = r1;
% r2(p2) = r2;
% r2(peakmagnitude==0)=ncomps;
% % r3(p3) = r3;
% 
% fin_score = r1 + r2;
% 
% % select the first 3 components in the scoring
% [~, p_final] = sort(fin_score); idxs_comps_kept = p_final(1:3);
% 
% %% generate ICAfield
% 
% ICAfield.full = [];
% ICAfield.alpha = [];
% 
% % normalize weights map in 0-1 range for each component
% absval_map = abs(EEG.icawinv);
% rescaled = (absval_map - min(absval_map))./(max(absval_map) - min(absval_map));
% 
% %% evaluate anterior-posterior gradient of components
% topo_comps = rescaled(1:64, idxs_comps_kept);
% 
% xyzpos = [cell2mat({EEG.chanlocs(:).X})', cell2mat({EEG.chanlocs(:).Y})',...
%           cell2mat({EEG.chanlocs(:).Z})'];
% 
% corr_comps = corr(topo_comps);
% 
% relative_pos = topo_comps' * xyzpos;
% 
% [~,postAntGrad] = sort(relative_pos(:,1),'ascend');
% 
% % resort components based on the gradient
% idxs_comps_kept = idxs_comps_kept(postAntGrad);
% 
% % ICAfield alpha
% ICAfield.alpha.idxs_comps_kept = idxs_comps_kept;
% ICAfield.alpharange = alpharange;
% ICAfield.topomaps_rescaled = rescaled(:, idxs_comps_kept);
% ICAfield.alpha.FFT.avg_spectra = avg_spctr(idxs_comps_kept, :);
% ICAfield.alpha.FFT.freqs = xfreq;
% 
% singletrial_FFT_range = pow_spctr(idxs_comps_kept, :, :);
% singletrial_FFT_range = squeeze(mean(singletrial_FFT_range, 2));
% ICAfield.alpha.FFT.singletrial = singletrial_FFT_range;
% 
% % wavelet sub sub field
% CFG = []; 
% CFG.min_freq = 2;
% CFG.max_freq = 35;
% CFG.num_frex = 40;
% CFG.frex = linspace(CFG.min_freq, CFG.max_freq, CFG.num_frex);
% CFG.fwhm_t = 2 .* 1./CFG.frex;
% CFG.nextpow2 = true;
% CFG.ITC = false;
% CFG.keeptrials = true;
% CFG.resample = true;
% 
% selected_comp_data = EEG.icaact(idxs_comps_kept, :, :);
% wv = optimized_wavelet(selected_comp_data, EEG.srate, CFG);
% 
% ICAfield.alpha.wv = wv;
% 
% new_srate = ceil(2.5*CFG.max_freq); %(1/2 over nyquist)
% str_t = EEG.times(1); end_t = EEG.times(end);
% new_time = linspace(str_t, end_t, ceil(new_srate*(end_t-str_t)/1e3))';
% ICAfield.alpha.wv.times = new_time;
% 
% % ICAfield full
% ICAfield.full.icaact = EEG.icaact;
% ICAfield.full.icawinv = EEG.icawinv;
% 
% 
% 
% 
% 
% end