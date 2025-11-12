%% ==============================================================
%  EMG Signal Processing Pipeline
%  Author: Camille Cunin
%  Date:   30/10/25
%  Description:
%     This script processes NEV data files by:
%       - Loading raw data
%       - Notch filtering (60 Hz + harmonics)
%       - Removing artefacts and bursts
%       - Extracting EMG band (1–50 Hz)
%       - Computing time-dependent EMG power
%       - Calculating EMG peak rates
%       - Power Spectral Density (PSD) analysis
%
%  Usage:
%     1. Set directory paths and filename in "User Parameters".
%     2. Run the script.
%
%  Dependencies:
%     - Neural Signal Processing (NS) library
%     - Custom functions (e.g., read_ripple_data.m)
%
% ==============================================================

clear; close all; clc;

%% -------------------- USER PARAMETERS ------------------------
% Directory settings (replace with your own paths)
dirHome = 'YOUR_HOME_DIRECTORY';
dirAllData = 'YOUR_DATA_DIRECTORY';
addpath(dirAllData)

% Data file name
fileName = 'rawdata.nev';
full_data_file_path = fullfile(dirAllData, fileName);

% MATLAB functions and NS library directories
Matlab_functions = dirHome;
addpath(Matlab_functions);
NSlibrary = dirHome;
addpath(NSlibrary)

%% ------------------- FOLDER STRUCTURE ------------------------
% Processed data directory
processedDataPath = 'YOUR_PROCESSED_DATA_DIRECTORY';

% Base name for folder structure (based on NEV file name)
baseFolderName = strrep(fileName, '_estim.nev', '');
baseFolderName = strrep(baseFolderName, '.nev', '');

% Create main folder for the experiment
mainFolderPath = fullfile(processedDataPath, baseFolderName);
if ~exist(mainFolderPath, 'dir')
    mkdir(mainFolderPath);
end

% Subfolders for processed data
subFolders = {'Raw data', 'PSD', 'Filtered data', 'Spectrogram', ...
              'Sigma thresholds', 'Rates'};

for i = 1:length(subFolders)
    subFolderPath = fullfile(mainFolderPath, subFolders{i});
    if ~exist(subFolderPath, 'dir')
        mkdir(subFolderPath);
    end
end

% Experiment display and results file
expName = baseFolderName; 
figName = baseFolderName;

% Define results text file path
fileNameTxt = strcat('results_', figName, '.txt');
fullFilePathTxt = fullfile(mainFolderPath, fileNameTxt);

% Specify electrodes to analyze
elecs = 1:7;

% Specify x- and y-axis boundaries
x_lim = [0 1200];   % seconds
y_lim = [-500 500]; % µV

%% ------------------- LOAD RAW DATA ---------------------------
[ns_RESULT, hFile] = ns_OpenFile(full_data_file_path);
[ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
[raw, t, fs]=read_ripple_data(hFile, elecs);
[numSamples, numChannels] = size(raw);

%% ------------------- NOISE FILTERING -------------------------------
% Remove 60 Hz and harmonics
harmonics = 1:7;                     % Harmonics to remove
notch_frequencies = 60 * harmonics;  % Harmonic frequencies
Q = 5;                               % Quality factor
Fs = fs;                             % Sampling frequency

rawWithout60HzFolder = fullfile(mainFolderPath, 'Filtered data', 'Raw without 60Hz');
if ~exist(rawWithout60HzFolder, 'dir')
    mkdir(rawWithout60HzFolder);
end

notch_filtered_data = zeros(size(raw));  

for ch = 1:numChannels
    signal = raw(:, ch);  % Get the data for the current channel
    
    % Notch filter each harmonic, including the 60Hz
    for i = 1:length(notch_frequencies)
        f0 = notch_frequencies(i);  % Current harmonic frequency
        
        % Design the notch filter for the current harmonic using designfilt
        notchFilter = designfilt('bandstopiir', 'FilterOrder', 2, ...
            'HalfPowerFrequency1', f0 - 1, 'HalfPowerFrequency2', f0 + 1, ...  % Adjust the frequency range if necessary
            'DesignMethod', 'butter', 'SampleRate', Fs);
        
        % Apply the notch filter
        signal = filter(notchFilter, signal);  % Apply the filter to the signal
    end
    
    % Store the filtered signal back in notch_filtered_data
    notch_filtered_data(:, ch) = signal;
end


%% ------------------ PULSE BURSTS REMOVAL ------------------ %%

raw_original = notch_filtered_data;  % or raw if you want to keep the 60Hz signal
raw_all_bursts_removed = raw_artefacts_removed;

burst_windows = [
    599, 606;   % Burst 1 (600)
    659, 666;   % Burst 2 (660)
    719, 726;   % Burst 3 (720)
    779, 786;   % Burst 4 (780)
    839, 846;   % Burst 5 (840)
    907, 914;   % Burst 6 (909)
    959, 966;   % Burst 7 (960)
    1019, 1026; % Burst 8 (1020)
    1079, 1086; % Burst 9 (1080)
    1139, 1146; % Burst 10 (1140)
    1199, 1206; % Burst 11 (1200)
];

disp('Manual burst windows:');
disp(burst_windows);

% ----- Zero out data within defined burst windows -----
for b = 1:size(burst_windows, 1)
    t_start_burst = burst_windows(b, 1);
    t_end_burst   = burst_windows(b, 2);

    % Convert time to sample indices
    idx_start = max(1, round(t_start_burst * fs));
    idx_end   = min(length(t), round(t_end_burst * fs));

    % Set those regions to zero across all channels
    raw_all_bursts_removed(idx_start:idx_end, :) = 0;
end


%% ------------------- EMG BAND ------------------------
b1 = 1; 
b2 = 50;
[b, a] = butter(2, [b1 b2] / (fs/2), 'bandpass');

emg_band = zeros(size(raw_all_bursts_removed));
for ch = 1:size(raw_all_bursts_removed,2)
    emg_band(:, ch) = filtfilt(b, a, raw_all_bursts_removed(:, ch));
end

emg_band = emg_band;  % [samples x channels]

emg_folder = fullfile(mainFolderPath, 'Filtered data', 'EMG');
if ~exist(emg_folder, 'dir'); mkdir(emg_folder); end

numChannels = length(elecs);
t = (0:size(emg_band,1)-1)/fs;

figure('Units','normalized','OuterPosition',[0 0 1 1]);
for ch = 1:numChannels
    subplot(numChannels,1,ch);
    plot(t, emg_band(:,ch), 'k'); hold on;

    for b = 1:size(burst_windows,1)
        burst_start = burst_windows(b,1);
        burst_end   = burst_windows(b,2);
        patch([burst_start burst_end burst_end burst_start], ...
              [-1000 -1000 1000 1000], 'k','FaceAlpha',1,'EdgeColor','none');
    end

    xlim(x_lim); ylim([-200 200]);
    set(gca, 'FontSize',14);
    ylabel(sprintf('Ch %d', elecs(ch)), 'FontSize',14);
    if ch == numChannels
        xlabel('Time (s)','FontSize',14);
    else
        set(gca,'XTickLabel',[]);
    end
end
sgtitle('EMG Band (1-50 Hz)','FontSize',16);

emg_band_path = fullfile(emg_folder, 'emg_band.png');
exportgraphics(gcf, emg_band_path,'Resolution',600);
close;


%% ------------------- EMG POWER OVER TIME ------------------------
emg_rect = abs(emg_band);

emg_pow_folder = fullfile(emg_folder, 'EMG Power');
if ~exist(emg_pow_folder,'dir'); mkdir(emg_pow_folder); end

window_RMS = round(0.25*fs);
step_size = round(0.05*fs);

emg_power = zeros(size(emg_rect));
for ch = 1:numChannels
    emg_power(:,ch) = movmean(emg_rect(:,ch).^2, window_RMS, 'Endpoints','discard');
end

t_power = (0:size(emg_power,1)-1)/fs;
down_idx = 1:step_size:length(t_power);
emg_power_ds = emg_power(down_idx,:);
t_power_ds = t_power(down_idx);

channels_to_plot = 2:7;
figure('Units','normalized','OuterPosition',[0 0 1 1]);
for i = 1:length(channels_to_plot)
    ch_idx = channels_to_plot(i);
    subplot(length(channels_to_plot),1,i);
    plot(t_power_ds, emg_power_ds(:,i),'k'); hold on;

    patch([t_drug_start t_drug_end t_drug_end t_drug_start], ...
          [-2000 -2000 5000 5000],'k','FaceAlpha',1,'EdgeColor','none');

    xlim(x_lim); ylim([0 2000]);
    set(gca,'FontSize',20,'XTickLabel',[],'YTickLabel',[]);
end
emg_power_path = fullfile(emg_pow_folder,'emg_power');
print(gcf,[emg_power_path,'.png'],'-dpng','-r600');
print(gcf,[emg_power_path,'.svg'],'-dsvg');


%% ------------------- EMG PEAK RATES ------------------------
win_sec = 1;
samples_per_win = win_sec*fs;
num_win = floor(size(emg_power,1)/samples_per_win);
sigma_thresh = [3 5 10];
t_full = (0:size(emg_power,1)-1)/fs;

all_rate_data = cell(numChannels,1);
global_max_rate = 0;

for c = 1:numChannels
    sm_power = emg_power(:,c);
    baseline = sm_power(t_full>=0 & t_full<=180);
    baseline_mean = mean(baseline);
    baseline_std  = std(baseline);
    
    rate_matrix = zeros(num_win,length(sigma_thresh));
    t_axis = zeros(1,num_win);

    for w = 1:num_win
        idx_start = (w-1)*samples_per_win + 1;
        idx_end   = idx_start + samples_per_win -1;
        if idx_end > length(sm_power); continue; end

        window_data = sm_power(idx_start:idx_end);
        t_axis(w) = mean(t_full(idx_start:idx_end));

        for s = 1:length(sigma_thresh)
            threshold = baseline_mean + sigma_thresh(s)*baseline_std;
            [~,locs] = findpeaks(window_data,'MinPeakHeight',threshold);
            rate_matrix(w,s) = length(locs)/win_sec;
            global_max_rate = max(global_max_rate, rate_matrix(w,s));
        end
    end
    all_rate_data{c} = struct('channel',elecs(c),'t_axis',t_axis,'rate_matrix',rate_matrix);
end

channels_to_plot = 2:7;
selected_data = all_rate_data(channels_to_plot);

combined_y_max = 0;
for c = 1:length(selected_data)
    combined_y_max = max(combined_y_max, max(selected_data{c}.rate_matrix(:)));
end
combined_y_max = combined_y_max*1.1;

figure('Units','normalized','OuterPosition',[0 0 1 1]);
for c = 1:length(selected_data)
    subplot(length(selected_data),1,c); hold on;
    data = selected_data{c};
    t_axis = data.t_axis;
    rate_matrix = data.rate_matrix;

    patch([t_drug_start t_drug_end t_drug_end t_drug_start], [-2000 -2000 5000 5000],'k','FaceAlpha',1,'EdgeColor','none');
    
    for s = 1:size(rate_matrix,2)
        plot(t_axis, rate_matrix(:,s),'LineWidth',1.5,'DisplayName',sprintf('%.0fσ',sigma_thresh(s)));
    end

    ylabel(sprintf('Ch %d',c),'FontSize',20);
    ylim([0 combined_y_max]); xlim([0 1800]);
    grid on; set(gca,'FontSize',18);
    if c == length(selected_data); xlabel('Time (s)','FontSize',20);
    else; set(gca,'XTickLabel',[]); end
end

sgtitle('EMG Peak Rate Over Time','FontSize',24,'FontWeight','bold');
emg_rate_path = fullfile(emg_folder,'emg_peak_rate');
print(gcf,[emg_rate_path,'.png'],'-dpng','-r600');
print(gcf,[emg_rate_path,'.svg'],'-dsvg');


%% ------------------ SPECTROGRAM 0-50 Hz ------------------ %%
ylim_vals = [0 50];
clim_vals = [-10 20];
channels_to_plot = [1 2 3 4 5 6];
x_start = 0;
x_end   = 1799;
fs_ds = 200;
ds_factor = round(30000 / fs_ds);

signal_ds = downsample(raw_all_bursts_removed, ds_factor);
signal_ds = signal_ds(:, channels_to_plot);
[numSamples, numChannels_plot] = size(signal_ds);

t_ds = (0:numSamples-1)/fs_ds;
idx_range = (t_ds >= x_start) & (t_ds <= x_end);

signal_window = signal_ds(idx_range,:);
t_window = t_ds(idx_range);

window   = round(0.5*fs_ds);
noverlap = round(0.45*fs_ds);
nfft     = 1024;

spectrogramFolder = fullfile(mainFolderPath,'Spectrogram');
if ~exist(spectrogramFolder,'dir'); mkdir(spectrogramFolder); end

fig = figure('Units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(numChannels_plot,1,'Padding','compact','TileSpacing','compact');

for idx = 1:numChannels_plot
    nexttile;
    signal = signal_window(:,idx);
    [S,F,T,P] = spectrogram(signal,window,noverlap,nfft,fs_ds);
    T = T + t_window(1);

    imagesc(T,F,10*log10(P+1e-12));
    axis xy; colormap('parula'); caxis(clim_vals); ylim(ylim_vals);

    ax = gca; ax.FontSize = 18; ax.LabelFontSizeMultiplier = 1.2;

    cb = colorbar;
    cb.Label.String = 'Power (dB)';
    cb.Label.FontSize = 16;
    cb.Label.FontWeight = 'bold';
    cb.FontSize = 16;

    for b = 1:size(burst_windows,1)
        burst_start = max(burst_windows(b,1),x_start);
        burst_end   = min(burst_windows(b,2),x_end);
        if burst_end < x_start || burst_start > x_end
            continue;
        end
        patch([burst_start burst_end burst_end burst_start],[ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)],...
            'red','FaceAlpha',0.9,'EdgeColor','none');
    end

    xlim([x_start x_end]);
    ylabel(sprintf('Ch %d', channels_to_plot(idx)),'FontSize',20);
    if idx == numChannels_plot
        xlabel('Time (s)','FontSize',20);
    else
        set(gca,'XTickLabel',[]);
    end
end
sgtitle('Spectrogram (0-50 Hz)','FontSize',24,'FontWeight','bold');

filename_base = fullfile(spectrogramFolder,'spectrogram_0-50Hz');
print(fig,[filename_base,'.png'],'-dpng','-r600');
print(fig,[filename_base,'.svg'],'-dsvg');

%% ------------------ SPECTROGRAM 0-15 kHz ------------------ %%
ylim_vals = [0 15000];
clim_vals = [-10 20];
channels_to_plot = [1 2 3 4 5 6];
x_start = 0;
x_end   = 1799;
fs_ds = 30000;

signal_window = raw_all_bursts_removed(:,channels_to_plot);
t_window = (0:size(signal_window,1)-1)/fs_ds;

idx_range = (t_window >= x_start) & (t_window <= x_end);
signal_window = signal_window(idx_range,:);
t_window = t_window(idx_range);

window   = round(1*fs_ds);
noverlap = round(0.5*fs_ds);
nfft     = 1024;

fig = figure('Units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(length(channels_to_plot),1,'Padding','compact','TileSpacing','compact');

for idx = 1:length(channels_to_plot)
    nexttile;
    signal = signal_window(:,idx);
    [S,F,T,P] = spectrogram(signal,window,noverlap,nfft,fs_ds);
    T = T + t_window(1);

    imagesc(T,F,10*log10(P+1e-12));
    axis xy; colormap('parula'); caxis(clim_vals); ylim(ylim_vals);

    ax = gca; ax.FontSize = 18; ax.LabelFontSizeMultiplier = 1.2;

    cb = colorbar;
    cb.Label.String = 'Power (dB)';
    cb.Label.FontSize = 16;
    cb.Label.FontWeight = 'bold';
    cb.FontSize = 16;

    for b = 1:size(burst_windows,1)
        burst_start = max(burst_windows(b,1),x_start);
        burst_end   = min(burst_windows(b,2),x_end);
        if burst_end < x_start || burst_start > x_end
            continue;
        end
        patch([burst_start burst_end burst_end burst_start],[ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)],...
            'red','FaceAlpha',0.9,'EdgeColor','none');
    end

    xlim([x_start x_end]);
    ylabel(sprintf('Ch %d', channels_to_plot(idx)),'FontSize',20);
    if idx == length(channels_to_plot)
        xlabel('Time (s)','FontSize',20);
    else
        set(gca,'XTickLabel',[]);
    end
end
sgtitle('Spectrogram (0-15 kHz)','FontSize',24,'FontWeight','bold');

filename_base = fullfile(spectrogramFolder,'spectrogram_0-15kHz');
print(fig,[filename_base,'.png'],'-dpng','-r600');
print(fig,[filename_base,'.svg'],'-dsvg');
