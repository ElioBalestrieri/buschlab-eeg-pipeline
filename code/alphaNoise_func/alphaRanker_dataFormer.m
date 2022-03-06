function [vect_min_comp] = alphaRanker_dataFormer(cfg, EEG)

cumsum_thresh = .8; % selection of components accounting for .8

%% compute spectra for the components

FFT_cfg = [];
FFT_cfg.alpharange = [8, 13];
FFT_cfg.expfit_freqrange = [3, 30];
FFT_cfg.timewin = [-1, 0] * 1e3;

[full_spectra, ...
 singletrial_alpharange, ...
 pruned_spectra] = local_FFT(FFT_cfg, EEG);

%% fit exponential decay function
% and determine first ranking criterion: peak amplitude (and presence) in
% the alpha range

out = local_fit_expmodel(pruned_spectra);

ncomps = size(out.detrended, 2);
peakmagnitude = zeros(ncomps, 1);

for icomp = 1:ncomps
    
    [pks, locs] = findpeaks(out.detrended(:, icomp), out.pruned_xfreqs, ...
                            'MinPeakDistance',2);

    lgcl_mask_peakpresent = locs>=min(FFT_cfg.alpharange) & ...
                            locs<=max(FFT_cfg.alpharange);
    
    if any(lgcl_mask_peakpresent) 
        
        this_peak = max(pks(lgcl_mask_peakpresent)); % added to exclude those bad fits where the peak was below the detrending line
        if this_peak < 0            
            this_peak = 0;
        end
        peakmagnitude(icomp) = this_peak;
        
    end
    
end

% assign the rank based on the peakmagnitude
[~, order] = sort(peakmagnitude, 'descend');
rank1 = (1:ncomps)';
rank1(order) = rank1;


%% rank components for posteriority gradient

% normalize weights map in 0-1 range for each componen
absval_map = abs(EEG.icawinv);
rescaled = (absval_map - min(absval_map))./(max(absval_map) - min(absval_map));

topo_comps = rescaled(1:64, :);

xyzpos = [cell2mat({EEG.chanlocs(:).X})', cell2mat({EEG.chanlocs(:).Y})',...
          cell2mat({EEG.chanlocs(:).Z})'];

relative_pos = topo_comps' * xyzpos;

[~,postAntGrad] = sort(relative_pos(:,1),'ascend');

rank2 = (1:ncomps)';
rank2(postAntGrad) = rank2;

%% create matrix

mat_comps = [(1:ncomps)', peakmagnitude, rank1, rank2, (rank1 + rank2)/2];
srtd_mat = sortrows(mat_comps, 5, 'ascend');

norm_peaks = srtd_mat(:, 2)/sum(srtd_mat(:, 2));
cumsum_alpha_pow = cumsum(norm_peaks);
mask_comp_sel = cumsum_alpha_pow<cumsum_thresh;

% reduce matrix
srtd_mat = srtd_mat(mask_comp_sel, :);

% weights for further selection
W1 = norm_peaks(mask_comp_sel); W1 = W1/sum(W1);
W2 = 1 ./ (1:sum(mask_comp_sel));

%%

idxs_kept = srtd_mat(:, 1);
cross_trials_corr = singletrial_alpharange.crosscorr(idxs_kept, idxs_kept);

st_alpha_selection = singletrial_alpharange.pow(idxs_kept, :)';
norm_st_alpha = zscore(st_alpha_selection, [], 1);

plot_summary = false;

if plot_summary

    figure; subplot(1, 2, 1) %#ok<UNRCH> 
    plot(pruned_spectra.freqs, pruned_spectra.pow(idxs_kept, :)');
    title('original spectra')

    subplot(1, 2, 2)
    plot(pruned_spectra.freqs, out.detrended(:, idxs_kept))
    title('detrended spectra')
    

end


%% regression model with all the alpha st powers

stimEn = EEG.trialinfo.energy_probed;
Ytrgt = abs(EEG.trialinfo.bin_resp - 2);

anyNAN = isnan(stimEn); stimEn(anyNAN) = []; Ytrgt(anyNAN) = [];
stimEn = zscore(stimEn); norm_st_alpha(anyNAN, :) = [];

logReg_mat = [stimEn, norm_st_alpha];

[B, ~, stats] = mnrfit(logReg_mat, Ytrgt);

%% compute wavelet transform on components selected

CFG = []; 
CFG.min_freq = 2;
CFG.max_freq = 35;
CFG.num_frex = 40;
CFG.frex = linspace(CFG.min_freq, CFG.max_freq, CFG.num_frex);
CFG.fwhm_t = 2 .* 1./CFG.frex;
CFG.nextpow2 = true;
CFG.ITC = false;
CFG.keeptrials = true;
CFG.resample = true;

selected_comp_data = EEG.icaact(idxs_kept, :, :);
wv = optimized_wavelet(selected_comp_data, EEG.srate, CFG);

%% create components based structure

ICAcomps = [];
ICAcomps.idxs_kept = idxs_kept;

ICAcomps.FFT.cfg = FFT_cfg;
ICAcomps.FFT.pruned_spectra = pruned_spectra.pow(idxs_kept, :)';
ICAcomps.FFT.detrended = out.detrended(:, idxs_kept);
ICAcomps.FFT.fullFreqs = full_spectra.freqs;
ICAcomps.FFT.prunedFreqs = pruned_spectra.freqs;
ICAcomps.FFT.cross_trials_corr = cross_trials_corr;
ICAcomps.FFT.srtd_mat = srtd_mat;

% single trials 
ICAcomps.FFT.ST.all = st_alpha_selection;

if size(st_alpha_selection, 2) >= 3
    ranked_subselect = tiedrank(st_alpha_selection(:, 1:3));
else
    ranked_subselect = tiedrank(st_alpha_selection(:, 1:end));
end

zscored_rnkd_sub = zscore(ranked_subselect, [], 1);
zscored_rnkd_sub = mean(zscored_rnkd_sub, 2);

% repeat the logreg
stimEn = EEG.trialinfo.energy_probed;
Ytrgt = abs(EEG.trialinfo.bin_resp - 2);

anyNAN = isnan(stimEn); stimEn(anyNAN) = []; Ytrgt(anyNAN) = [];
stimEn = zscore(stimEn); zscored_rnkd_sub(anyNAN) = [];

logReg_mat_ranksub = [stimEn, zscored_rnkd_sub];

[B_ranksub, ~, stats_ranksub] = mnrfit(logReg_mat_ranksub, Ytrgt);

ICAcomps.FFT.ST.fin_select.zscrd_rnkd = zscored_rnkd_sub;
ICAcomps.FFT.ST.fin_select.B = B_ranksub(3);
ICAcomps.FFT.ST.fin_select.p = stats_ranksub.p(3);
ICAcomps.FFT.ST.W1 = W1;
ICAcomps.FFT.ST.W2 = W2;


disp(B_ranksub(3))
disp(stats_ranksub.p(3))

pause(5)

%%

reduced_B = B(3:end); reduced_p = stats.p(3:end);
array_logreg = [idxs_kept, reduced_B, reduced_p]; % remove intercept and stimulus energy
ICAcomps.FFT.tableLogReg = array2table(array_logreg, 'VariableNames', ...
    {'orig_idx', 'beta', 'p'});

disp(ICAcomps.FFT.tableLogReg); 

ICAcomps.WV = wv;
ICAcomps.WV.cfg = CFG;

ICAcomps.topo = rescaled(1:64, idxs_kept);

% pipe out ouput to have a population summary
idxs_minor = reduced_B == min(reduced_B);
B_minor = reduced_B(idxs_minor);
p_minor = reduced_p(idxs_minor);

vect_min_comp = [find(idxs_minor), B_minor, p_minor];

%% create final structure

% parse subject name and define output folder
subjname = EEG.urname(1:4);

SUBJ = [];

% append ICA comps
SUBJ.ICAcomps = ICAcomps;

% general
SUBJ.name = subjname;
SUBJ.beh = EEG.trialinfo;

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

SUBJ.chanlocs = EEG.chanlocs(1:64);

path_to_finalfile = fullfile(cfg.dir.main, 'Preprocessing', 'preprocessed', 'prestim_comps');

if ~isfolder(path_to_finalfile)
    mkdir(path_to_finalfile);
end
    
save(fullfile(path_to_finalfile, [subjname '.mat']), 'SUBJ');

end


%% ################### LOCAL FUNCTIONS #####################

function [full, st_alpha, pruned_1overFfit] = local_FFT(FFT_cfg, EEG)

prestim = FFT_cfg.timewin;

% subselect time window of interest
lgcl_time = EEG.times>=min(prestim) &  EEG.times<max(prestim);
comp_prestim = EEG.icaact(:, lgcl_time, :);

% define signal length
lsign = sum(lgcl_time); length_fft = 2^nextpow2(lsign); % for speed, if necessary

% apply fft
cmplx = fft(comp_prestim, length_fft, 2);

% compute power
pow_spctr = (abs(cmplx)/length_fft).^2;    
pow_spctr = pow_spctr(:, 1:(ceil(length_fft/2)+1), :);
fin_signal_length = size(pow_spctr, 2);

xfreq = linspace(0, EEG.srate/2, fin_signal_length);

% assign output structures

% full
full.freqs = xfreq;
full.avg_spctr = mean(pow_spctr, 3);

% single trial alpha
lgcl_alpha = (xfreq >= min(FFT_cfg.alpharange)) & ...
                       (xfreq <= max(FFT_cfg.alpharange));

st_alpha.freqs = xfreq(lgcl_alpha);
st_alpha.range = FFT_cfg.alpharange;
st_alpha.pow = squeeze(mean(pow_spctr(:, lgcl_alpha, :), 2));
st_alpha.crosscorr = corr(st_alpha.pow');

% pruned to for 1/f fit

lgcl_exprange = (xfreq>=min(FFT_cfg.expfit_freqrange)) & ...
                (xfreq<=max(FFT_cfg.expfit_freqrange));

pruned_1overFfit.freqs = xfreq(lgcl_exprange);
pruned_1overFfit.pow = mean(pow_spctr(:, lgcl_exprange, :), 3);

end


function out = local_fit_expmodel(pruned_spectra)

Xfit = pruned_spectra.freqs';
Yfit = double(pruned_spectra.pow');

anon_expfunc = @(v, x) v(1).*exp(-x./v(2)) + v(3);

ncomps = size(Yfit, 2);
nfreqs = length(Xfit);

[mat_detrended_spectra, fitted_curves] = deal(nan(nfreqs, ncomps));
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



