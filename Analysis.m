%% Data analysis of colonic EMG recordings from Ripple
clear
close all
clc

dirHome = [ ];
dirAllData= [dirHome];
addpath(dirAllData)
Matlab_functions =[dirHome];
addpath(Matlab_functions);
NSlibrary= [ ];
addpath(NSlibrary)

% Load NEV file
fileName = ; 
[ns_RESULT, hFile] = ns_OpenFile(fileName);
[ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
numEnt = length(hFile.Entity);
fs_clock = 30000;

%% Enter experiment details
t_drug  = [];
t_start = 0;
t_end   = 1800;

expName = ;
figName = ;
folderPath = ;

if ~exist(folderPath, 'dir'); mkdir(folderPath); end
fileNameTxt     = strcat('results_',figName,'.txt');
fullFilePathTxt = fullfile(folderPath, fileNameTxt);

elecs = [2, 3, 4, 5, 6, 7];

%% Extract raw data for chosen electrodes
for i=1:length(elecs)
    estr_raw(i) = strcat('raw',{' '},num2str(elecs(i)));
end

idx = [];
for ii = 1:numEnt
    if ismember(hFile.Entity(ii).EntityType,{'Analog'}) && ismember(hFile.Entity(ii).Label, estr_raw)
        idx(ii,:) = 1;
    else
        idx(ii,:) = 0;
    end
end
idx = find(idx == 1);

[~, raw_nsAnalogInfo] = ns_GetAnalogInfo(hFile, idx(1));
fs = raw_nsAnalogInfo.SampleRate;

raw = [];
for i = 1:length(idx)
    [~, ~, rawData] = ns_GetAnalogData(hFile, idx(i), 1, 1e8);
    raw(:, i) = rawData;
end

t = (0:(size(raw, 1) - 1)) / fs;
raw_original = raw;

%% EMG-band filtering (1–50 Hz, 2nd-order, zero-phase)
[b,a]=butter(2,[1 50]/(fs/2),'bandpass');
spkband = zeros(size(raw));
for i=1:size(raw,2)
    spkband(:,i)= filtfilt(b,a,raw(:,i));
end

%% Envelope power (250 ms) and 50 ms step (20 Hz)
env_win = round(0.25*fs);
sig_power = movmean(abs(spkband).^2, env_win, 1, 'omitnan');
hop = round(0.05*fs);
sig_power = sig_power(1:hop:end, :);
t_ds      = t(1:hop:end);
fs_ds     = fs / hop;

%% Peak-rate: baseline μ,σ (0–180 s) and kσ thresholds; 1 s bins
base_mask = (t_ds >=  & t_ds <= );
mu = mean(sig_power(base_mask,:), 1);
sg = std( sig_power(base_mask,:), 0, 1);
k_list  = [3 5 10]; % Choose k value
minDist = round(0.100 * fs_ds);
analysis_mask = (t_ds >= t_start & t_ds <= t_end);

elecs_used_env = elecs;
[~, col_map_env] = ismember(elecs_used_env, elecs);

edges = t_start:1:t_end;
if edges(end) < t_end, edges = [edges t_end]; end
nBins = numel(edges)-1;

peak_times_env = cell(numel(elecs_used_env), numel(k_list));
for ci = 1:numel(elecs_used_env)
    col = col_map_env(ci);
    x = sig_power(:, col);
    x(~analysis_mask) = -inf;
    for kk = 1:numel(k_list)
        thr = mu(col) + k_list(kk) * sg(col);
        [~, locs] = findpeaks(x, 'MinPeakHeight', thr, 'MinPeakDistance', minDist);
        peak_times_env{ci, kk} = t_ds(locs);
    end
end

file_time_env = fullfile(folderPath, sprintf('ENV_time_binned_counts_%s.xlsx', figName));
if exist(file_time_env,'file'); delete(file_time_env); end
for kk = 1:numel(k_list)
    CountsMat = nan(numel(elecs_used_env), nBins);
    for ci = 1:numel(elecs_used_env)
        pt = peak_times_env{ci, kk};
        CountsMat(ci,:) = histcounts(pt, edges);
    end
    CH = []; B0 = []; B1 = []; CNT = [];
    for ci = 1:numel(elecs_used_env)
        cnt = CountsMat(ci,:);
        CH  = [CH;  repmat(elecs_used_env(ci), nBins, 1)];
        B0  = [B0;  edges(1:end-1)'];
        B1  = [B1;  edges(2:end)'];
        CNT = [CNT; cnt(:)];
    end
    T_time_env = table();
    T_time_env.Channel   = CH;
    T_time_env.BinStartS = B0;
    T_time_env.BinEndS   = B1;
    T_time_env.Count     = CNT;
    T_time_env.SigmaK    = repmat(k_list(kk), numel(CH), 1);
    thr_vec = (mu(col_map_env) + k_list(kk)*sg(col_map_env));
    ThrCol = zeros(size(CH));
    BMCol  = zeros(size(CH));
    BSCol  = zeros(size(CH));
    for ci = 1:numel(elecs_used_env)
        mask = (CH == elecs_used_env(ci));
        ThrCol(mask) = thr_vec(ci);
        BMCol(mask)  = mu(col_map_env(ci));
        BSCol(mask)  = sg(col_map_env(ci));
    end
    T_time_env.Thr_Power   = ThrCol;
    T_time_env.BaseMeanPow = BMCol;
    T_time_env.BaseStdPow  = BSCol;
    writetable(T_time_env, file_time_env, 'Sheet', sprintf('k%ds', k_list(kk)));
end

%% PSD (Welch; 2 s, 50% overlap, NFFT=2^nextpow2(N)) on spkband; 0–10 Hz and band powers
baseline_win = [ ];   % [t0 t1] rows
stim_windows = [ ];   % [t0 t1] rows

win   = fs * 2;
ovlp  = fs * 1;
nfft  = 2^nextpow2(win);

X = spkband;
tvec = (0:size(X,1)-1)/fs;

nBase = size(baseline_win, 1);
nStim = size(stim_windows, 1);

all_base_each = [];
all_stim      = [];

for ch = 1:size(X,2)
    x = X(:, ch);
    tmp_base = zeros(nfft/2+1, nBase);
    for b = 1:nBase
        mask = (tvec >= baseline_win(b,1)) & (tvec <= baseline_win(b,2));
        [Pxx, F] = pwelch(x(mask), win, ovlp, nfft, fs, "onesided");
        tmp_base(:, b) = Pxx;
    end
    all_base_each(:, :, ch) = tmp_base;
    for s = 1:nStim
        mask = (tvec >= stim_windows(s,1)) & (tvec <= stim_windows(s,2));
        [Pxx, ~] = pwelch(x(mask), win, ovlp, nfft, fs, "onesided");
        all_stim(:, s, ch) = Pxx;
    end
end

mean_base_lin = mean(all_base_each, 3, "omitnan");
mean_base_lin = mean(mean_base_lin, 2, "omitnan");
mean_stim_lin = mean(all_stim, 3, "omitnan");

fmask = (F >= 0) & (F <= 10);
F10   = F(fmask);
base10 = mean_base_lin(fmask);
stim10 = mean_stim_lin(fmask, :);

bands = [0 2; 2 5; 5 10];
band_labels = ["P_0_2", "P_2_5", "P_5_10"];

BaseBand = zeros(1, size(bands,1));
for bi = 1:size(bands,1)
    bm = (F10 >= bands(bi,1)) & (F10 < bands(bi,2));
    BaseBand(bi) = trapz(F10(bm), base10(bm));
end

StimBand = zeros(nStim, size(bands,1));
for s = 1:nStim
    for bi = 1:size(bands,1)
        bm = (F10 >= bands(bi,1)) & (F10 < bands(bi,2));
        StimBand(s, bi) = trapz(F10(bm), stim10(bm, s));
    end
end

T_psd = table(F10, 10*log10(base10 + eps), 'VariableNames', {'F_Hz','Baseline_dBperHz'});
for s = 1:nStim
    T_psd.(['Stim' num2str(s) '_dBperHz']) = 10*log10(stim10(:,s) + eps);
end

T_band = table();
T_band.Metric = ["Baseline"; strcat("Stim", string((1:nStim)'))];
for bi = 1:size(bands,1)
    col = [BaseBand(bi); StimBand(:,bi)];
    T_band.(band_labels(bi)) = col;
end

outPath_psd  = fullfile(folderPath, sprintf('PSD_0_10Hz_%s.xlsx', figName));
outPath_band = fullfile(folderPath, sprintf('BandPowers_0_10Hz_%s.xlsx', figName));
if exist(outPath_psd,  'file'); delete(outPath_psd);  end
if exist(outPath_band, 'file'); delete(outPath_band); end
writetable(T_psd,  outPath_psd,  'Sheet','PSD_0_10Hz');
writetable(T_band, outPath_band, 'Sheet','BandPowers');

fprintf('Saved: %s\n', outPath_psd);
fprintf('Saved: %s\n', outPath_band);

%% ===== Latency to first 3σ burst =====
k_target = 3; 
k_idx = find(k_list == k_target, 1);

stim_groups = {'STIM1','STIM2','STIM3'};
stim_onsets = [];
stim_tags   = strings(0);
for g = 1:numel(stim_groups)
    name = stim_groups{g};
    if exist('SEG_RANGES','var') && isfield(SEG_RANGES,name)
        ranges = SEG_RANGES.(name);
        for sidx = 1:numel(ranges)
            ts = ranges{sidx}(1);
            stim_onsets(end+1) = ts; %#ok<SAGROW>
            stim_tags(end+1)   = sprintf('%s_%d', name, sidx); %#ok<SAGROW>
        end
    end
end

latency_mat = nan(numel(elecs_used_env), numel(stim_onsets));
for ci = 1:numel(elecs_used_env)
    pt = peak_times_env{ci, k_idx};
    for jj = 1:numel(stim_onsets)
        ts = stim_onsets(jj);
        idx = find(pt >= ts, 1, 'first');
        if ~isempty(idx)
            latency_mat(ci, jj) = pt(idx) - ts;
        end
    end
end

file_latency = fullfile(folderPath, sprintf('ENV_latency_firstPeak_%dsig_%s.xlsx', k_target, figName));
if exist(file_latency,'file'); delete(file_latency); end

T_lat = table();
T_lat.Channel = elecs_used_env(:);
for jj = 1:numel(stim_onsets)
    T_lat.(sprintf('Latency_%s_s', stim_tags(jj))) = latency_mat(:, jj);
end
T_lat.MedianLatency_s = median(latency_mat, 2, 'omitnan');
T_lat.MeanLatency_s   = mean(latency_mat, 2, 'omitnan');

writetable(T_lat, file_latency, 'Sheet', 'latency_firstPeak');
fprintf('Latency(3σ, first burst) saved: %s\n', file_latency);
