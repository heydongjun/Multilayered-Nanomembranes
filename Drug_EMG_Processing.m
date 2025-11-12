%% ==============================================================
%  EMG Signal Processing Pipeline
%  Author: Camille Cunin
%  Date:   30/10/25
%  Description:
%     This script processes NEV data files by:
%       - Loading raw data
%       - Notch filtering (60 Hz + harmonics)
%       - Removing drug injection window
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
fileName = 'rawdata.nev'; %raw data file name
full_data_file_path = fullfile(dirAllData, fileName);

% MATLAB functions and NS library directories
Matlab_functions = dirHome;
addpath(Matlab_functions);
NSlibrary = dirHome;
addpath(NSlibrary)

%% ------------------- FOLDER STRUCTURE ------------------------
% Base path for processed data
processedDataPath = 'YOUR_PROCESSED_DATA_DIRECTORY';

% Base name for folder structure (based on NEV file name)
baseFolderName = strrep(fileName, '_drug.nev', '');
baseFolderName = strrep(baseFolderName, '.nev', '');

% Create main folder for the experiment
mainFolderPath = fullfile(processedDataPath, baseFolderName);
if ~exist(mainFolderPath, 'dir')
    mkdir(mainFolderPath);
end

% Create subfolders for processed data
subFolders = {'Raw data', 'PSD', 'Filtered data', 'Spectrogram', ...
              'Sigma thresholds', 'Rates'};

for i = 1:length(subFolders)
    subFolderPath = fullfile(mainFolderPath, subFolders{i});
    if ~exist(subFolderPath, 'dir')
        mkdir(subFolderPath);
    end
end

% Experiment display and results file setup
expName = baseFolderName; 
figName = baseFolderName;

% Define results text file path
fileNameTxt = strcat('results_', figName, '.txt');
fullFilePathTxt = fullfile(mainFolderPath, fileNameTxt);

% Specify electrodes to analyze
elecs = 1:7;

% Specify x- and y-axis boundaries
x_lim = [0 600];
y_lim = [-500 500];

%% ------------------- LOAD RAW DATA ---------------------------
[ns_RESULT, hFile] = ns_OpenFile(full_data_file_path);
[ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
[raw, t, fs]=read_ripple_data(hFile, elecs);
[numSamples, numChannels] = size(raw);


%% ------------------- NOISE FILTERING -------------------------------
harmonics = 1:7;
notch_frequencies = 60 * harmonics;
Q = 5;
Fs = fs;

rawWithout60HzFolder = fullfile(mainFolderPath, 'Filtered data', 'Raw without 60Hz');
if ~exist(rawWithout60HzFolder, 'dir')
    mkdir(rawWithout60HzFolder);
end

notch_filtered_data = zeros(size(raw));

for ch = 1:numChannels
    signal = raw(:, ch);
    for i = 1:length(notch_frequencies)
        f0 = notch_frequencies(i);
        notchFilter = designfilt('bandstopiir', 'FilterOrder', 2, ...
            'HalfPowerFrequency1', f0 - 1, 'HalfPowerFrequency2', f0 + 1, ...
            'DesignMethod', 'butter', 'SampleRate', Fs);
        signal = filter(notchFilter, signal);
    end
    notch_filtered_data(:, ch) = signal;
end



%% ------------------- DRUG INJECTION WINDOW REMOVAL ------------------------
t_drug_start = 600;
t_drug_end = 605;

notch_filtered_data_with_drug_zeroed = notch_filtered_data;
drug_idx_start = find(t >= t_drug_start, 1, 'first');
drug_idx_end = find(t <= t_drug_end, 1, 'last');
notch_filtered_data_with_drug_zeroed(drug_idx_start:drug_idx_end, :) = 0;

channels_to_plot = elecs;
[~, ch_indices] = ismember(channels_to_plot, elecs);

fig = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tlo = tiledlayout(numChannels, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tlo, 'Raw Data without 60Hz and with Drug Window Zeroed — Channels 1, 3, 6', ...
      'FontSize', 16, 'FontWeight', 'bold');

for i = 1:length(ch_indices)
    ch = ch_indices(i);
    ax = nexttile;
    plot(t, notch_filtered_data_with_drug_zeroed(:, ch), 'k'); hold on;
    patch([t_drug_start, t_drug_end, t_drug_end, t_drug_start], ...
          [-1000, -1000, 1000, 1000], 'black', 'FaceAlpha', 1, 'EdgeColor', 'none');
    xlim(x_lim);
    ylim(y_lim);
    set(ax, 'FontSize', 14);
    ylabel(sprintf('Ch %d', elecs(ch)), 'FontSize', 14);
    if i == length(ch_indices)
        xlabel('Time (s)', 'FontSize', 14);
    else
        set(ax, 'XTickLabel', []);
    end
end

outputFolder = fullfile(mainFolderPath, 'Raw data');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

filename_selected = fullfile(outputFolder, 'raw_data_drug_zeroed');
print(fig, '-dpng', '-r600', [filename_selected, '.png']);
print(fig, '-dsvg', [filename_selected, '.svg']);
close(fig);


%% ------------------- EMG BAND ------------------------
emg_band = [];
[b, a] = butter(2, [1 50] / (fs / 2), 'bandpass');
for i = 1:size(notch_filtered_data_with_drug_zeroed, 2)
    emg_band(i, :) = filtfilt(b, a, notch_filtered_data_with_drug_zeroed(:, i));
end
emg_band = emg_band';

emg_folder = fullfile(mainFolderPath, 'Filtered data', 'EMG');
if ~exist(emg_folder, 'dir')
    mkdir(emg_folder);
end

figure;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

for ch = 1:numChannels
    ax = subplot(numChannels, 1, ch);
    plot(t, emg_band(:, ch), 'k'); hold on;
    patch([t_drug_start, t_drug_end, t_drug_end, t_drug_start], ...
          [-1000, -1000, 1000, 1000], 'black', ...
          'FaceAlpha', 1, 'EdgeColor', 'none');
    xlim(x_lim);
    ybound = 100;
    ylim([-ybound ybound]);
    ylabel(sprintf('Ch %d', ch), 'FontSize', 20);
    set(ax, 'FontSize', 20);
    if ch == numChannels
        xlabel('Time (s)', 'FontSize', 20);
    else
        set(ax, 'XTickLabel', []);
    end
end

sgtitle('EMG Band (1–50 Hz)', 'FontSize', 24, 'FontWeight', 'bold');

baseFilename = fullfile(emg_folder, sprintf('EMG_band_%d', ybound));
print(gcf, [baseFilename, '_highres.png'], '-dpng', '-r600');
print(gcf, [baseFilename, '.svg'], '-dsvg');
close;


%% ------------------- EMG POWER OVER TIME ------------------------
emg_rectified = abs(emg_band);

emg_pow_folder = fullfile(emg_folder, 'EMG Power');
if ~exist(emg_pow_folder, 'dir'); mkdir(emg_pow_folder); end

window_size_sec = 0.25;
window_size_RMS = round(window_size_sec * fs);
step_size = round(0.05 * fs);
num_channels = size(emg_rectified, 2);

emg_power = cell(num_channels,1);
for ch = 1:num_channels
    emg_power{ch} = movmean(emg_rectified(:, ch).^2, window_size_RMS, 'Endpoints', 'discard');
end

emg_power = cat(2, emg_power{:});
n_samples = size(emg_power, 1);
t_power = (0:n_samples-1) / fs;

downsample_idx = 1:step_size:n_samples;
emg_power_ds = emg_power(downsample_idx, :);
t_power_ds = t_power(downsample_idx);

fig = figure;
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
channels_to_plot = [2, 3, 4, 5, 6, 7];

for ch_idx = 1:length(channels_to_plot)
    subplot(length(channels_to_plot), 1, ch_idx);
    plot(t_power_ds, emg_power_ds(:, ch_idx), 'k'); hold on;
    patch([t_drug_start, t_drug_end, t_drug_end, t_drug_start], ...
          [-2000, -2000, 5000, 5000], 'black', ...
          'FaceAlpha', 1, 'EdgeColor', 'none');
    xlim(x_lim);
    ylim([0, 2000]);
    set(gca, 'FontSize', 20, 'XTickLabel', [], 'YTickLabel', []);
end

filename_all = fullfile(emg_pow_folder, 'emg_power');
print(fig, [filename_all, '.png'], '-dpng', '-r600');
print(fig, [filename_all, '.svg'], '-dsvg');


%% ------------------- EMG PEAK RATES ------------------------
win_sec = 1;
samples_per_win = win_sec * fs;
num_win = floor(length(emg_power(:,1)) / samples_per_win);
sigma_thresholds = [3, 5, 10];
t_full = (0:length(emg_power)-1) / fs;
global_max_rate = 0;
all_rate_data = cell(length(elecs), 1);
fig = figure;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

for c = 1:length(elecs)
    ch = elecs(c);
    smoothed_power_vals = emg_power(:, c);
    baseline_idx = t_full >= 0 & t_full <= 180;
    baseline_power = smoothed_power_vals(baseline_idx);
    baseline_mean = mean(baseline_power);
    baseline_std = std(baseline_power);
    rate_matrix = zeros(num_win, length(sigma_thresholds));
    t_axis = zeros(1, num_win);
    for i = 1:num_win
        idx_start = (i - 1) * samples_per_win + 1;
        idx_end = idx_start + samples_per_win - 1;
        if idx_end > length(smoothed_power_vals)
            continue;
        end
        window_data = smoothed_power_vals(idx_start:idx_end);
        t_axis(i) = mean(t_full(idx_start:idx_end));
        for t_idx = 1:length(sigma_thresholds)
            threshold = baseline_mean + sigma_thresholds(t_idx) * baseline_std;
            [~, locs] = findpeaks(window_data, 'MinPeakHeight', threshold);
            rate = length(locs) / win_sec;
            rate_matrix(i, t_idx) = rate;
            if rate > global_max_rate
                global_max_rate = rate;
            end
        end
    end
    all_rate_data{c} = struct('channel', ch, 't_axis', t_axis, 'rate_matrix', rate_matrix);
end

y_shade_min = 0;
y_shade_max = global_max_rate * 1.1;
channels_to_plot = 2:7;
selected_data = {};
for i = 1:length(all_rate_data)
    if ismember(all_rate_data{i}.channel, channels_to_plot)
        selected_data{end+1} = all_rate_data{i};
    end
end

fig = figure;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
combined_y_min = 0;
combined_y_max = 0;
for c = 1:length(selected_data)
    rate_matrix = selected_data{c}.rate_matrix;
    combined_y_max = max(combined_y_max, max(rate_matrix(:)));
end
combined_y_max = combined_y_max * 1.1;

for c = 1:length(selected_data)
    ax = subplot(length(selected_data), 1, c);
    hold on;
    data = selected_data{c};
    t_axis = data.t_axis;
    rate_matrix = data.rate_matrix;
    patch([t_drug_start, t_drug_end, t_drug_end, t_drug_start], [-2000, -2000, 5000, 5000], 'black', 'FaceAlpha', 1, 'EdgeColor', 'none');
    for t_idx = 1:size(rate_matrix, 2)
        plot(t_axis, rate_matrix(:, t_idx), 'LineWidth', 1.5, 'DisplayName', sprintf('%.0fσ', sigma_thresholds(t_idx)));
    end
    ylabel(sprintf('Ch %d', c), 'FontSize', 20);
    ylim([combined_y_min combined_y_max]);
    xlim([0 1800]);
    grid on;
    set(ax, 'FontSize', 18);
    if c == length(selected_data)
        xlabel('Time (s)', 'FontSize', 20);
    else
        set(ax, 'XTickLabel', []);
    end
end

sgtitle('EMG Peak Rate Over Time', 'FontSize', 24, 'FontWeight', 'bold');
filename_all = fullfile(emg_folder, 'emg_peak_rate');
print(fig, [filename_all, '.png'], '-dpng', '-r600');
print(fig, [filename_all, '.svg'], '-dsvg');


%% ------------------ SPECTROGRAM 0-50 Hz ------------------ %%
ylim_vals = [0 50];
clim_vals = [-10 20];
expName = '';
x_start = 0;
x_end   = 1799;
channels_to_plot = [2, 3, 4, 5, 6, 7];
display_labels = 1:6;
t_drug_start = 600;  
t_drug_end   = 605;  
fs_ds = 200;
ds_factor = round(30000 / fs_ds);
signal_ds = downsample(notch_filtered_data_with_drug_zeroed, ds_factor);
[numSamples, numChannels] = size(signal_ds);
window = round(0.5 * fs_ds);
noverlap = round(0.45 * fs_ds);
nfft = 1024;

spectrogramFolder = fullfile(mainFolderPath, 'Spectrogram');
if ~exist(spectrogramFolder, 'dir')
    mkdir(spectrogramFolder);
end

fig = figure;
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
numChannels_plot = length(channels_to_plot);
tiledlayout(numChannels_plot, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

for idx = 1:numChannels_plot
    ch = channels_to_plot(idx);
    nexttile;
    signal = signal_ds(:, idx);
    [S, F, T, P] = spectrogram(signal, window, noverlap, nfft, fs_ds);
    imagesc(T, F, 10*log10(P + 1e-12));
    axis xy;
    colormap('parula');
    caxis(clim_vals);
    ylim(ylim_vals);
    ax = gca;
    ax.FontSize = 18;
    ax.LabelFontSizeMultiplier = 1.2;
    cb = colorbar;
    cb.Label.String = 'Power (dB)';
    cb.Label.FontSize = 16;
    cb.Label.FontWeight = 'bold';
    cb.FontSize = 16;
    y_min = ylim_vals(1);
    y_max = ylim_vals(2);
    patch([t_drug_start t_drug_end t_drug_end t_drug_start], ...
          [y_min y_min y_max y_max], 'red', ...
          'FaceAlpha', 0.9, 'EdgeColor', 'none');
    xlim([x_start, x_end]);
    ylabel(sprintf('Ch %d', idx), 'FontSize', 20);
    set(gca, 'FontSize', 20);
    if idx == numChannels_plot
        xlabel('Time (s)', 'FontSize', 20);
    else
        set(gca, 'XTickLabel', []);
    end
end

sgtitle('Spectrogram (0–50 Hz)', 'FontSize', 20, 'FontWeight', 'bold');
filename_base = fullfile(spectrogramFolder, 'Spectrogram_0-50Hz');
print(fig, '-dpng', '-r300', [filename_base '.png']);
print(fig, '-dsvg', [filename_base '.svg']);


%% ------------------ SPECTROGRAM 0-15 kHz ------------------ %%
ylim_vals = [0 15000];
clim_vals = [-10 20];
expName = '';

channels_to_plot = [2, 3, 4, 5, 6, 7];
display_labels = 1:6;

x_start = 0;
x_end   = 1799;

fs_ds = 30000;
signal_ds = notch_filtered_data_with_drug_zeroed(:, 1:6);
[numSamples, numChannels] = size(signal_ds);

numChannels_plot = numChannels;

if exist('display_labels','var') && ~isempty(display_labels)
    if length(display_labels) ~= numChannels_plot
        warning('display_labels length (%d) does not match number of selected channels (%d). Using channel numbers instead.', ...
                length(display_labels), numChannels_plot);
        display_labels_use = channels_to_plot;
    else
        display_labels_use = display_labels;
    end
else
    display_labels_use = channels_to_plot;
end

t_ds = (0:numSamples-1) / fs_ds;

idx_range = (t_ds >= x_start) & (t_ds <= x_end);
if ~any(idx_range)
    error('Selected time window [%.1f %.1f] s contains no samples. Check x_start/x_end.', x_start, x_end);
end
signal_ds_window = signal_ds(idx_range, :);
t_window = t_ds(idx_range);
numSamples_window = size(signal_ds_window, 1);

window   = round(1 * fs_ds);
noverlap = round(0.5 * fs_ds); 
nfft     = 1024;

spectrogramFolder = fullfile(mainFolderPath, 'Spectrogram');
if ~exist(spectrogramFolder, 'dir')
    mkdir(spectrogramFolder);
end

fig = figure;
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tiledlayout(numChannels_plot, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

for idx = 1:numChannels_plot
    ch = channels_to_plot(idx);
    nexttile;
    signal = signal_ds_window(:, idx);
    [S, F, T, P] = spectrogram(signal, window, noverlap, nfft, fs_ds);
    imagesc(T, F, 10*log10(P + 1e-12));
    axis xy;
    colormap('parula');
    caxis(clim_vals);
    ylim(ylim_vals);
    ax = gca;
    ax.FontSize = 18;
    ax.LabelFontSizeMultiplier = 1.2;
    cb = colorbar;
    cb.Label.String = 'Power (dB)';
    cb.Label.FontSize = 16;
    cb.Label.FontWeight = 'bold';
    cb.FontSize = 16;
    y_min = ylim_vals(1);
    y_max = ylim_vals(2);
    patch([t_drug_start t_drug_end t_drug_end t_drug_start], ...
          [y_min y_min y_max y_max], 'red', ...
          'FaceAlpha', 0.9, 'EdgeColor', 'none');
    xlim([x_start, x_end]);
    ylabel(sprintf('Ch %d', idx), 'FontSize', 20);
    set(gca, 'FontSize', 20);
    if idx == numChannels_plot
        xlabel('Time (s)', 'FontSize', 20);
    else
        set(gca, 'XTickLabel', []);
    end
end

sgtitle('Spectrogram (0-15 kHz)', 'FontSize', 24, 'FontWeight', 'bold');

filename_base = fullfile(spectrogramFolder, 'Spectrogram_0-15kHz');
print(fig, [filename_base, '.png'], '-dpng', '-r600');
print(fig, [filename_base, '.svg'], '-dsvg');
