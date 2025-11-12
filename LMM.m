%% LMM per-frequency for Power
clear; close all; clc;

%% 0) Settings
fname   = ;            
sheet   = 1;
alphaPV = 0.05;
doLogPower = true;


saveFolder = ;        
if ~exist(saveFolder,'dir'); mkdir(saveFolder); end

%% 1) Load data & basic types 
T = readtable(fname, 'Sheet', sheet);

needVars = {'Animal','Electrode','Condition','Frequency','Power'};
missing  = setdiff(needVars, T.Properties.VariableNames);
assert(isempty(missing), "Required columns are missing in the Excel file: %s", strjoin(missing, ', '));

T.Animal    = categorical(string(T.Animal));
T.Electrode = categorical(string(T.Electrode));
T.Condition = categorical(string(T.Condition));
T.Frequency = double(T.Frequency);
T.Power     = double(T.Power);


if ismember('Segment', T.Properties.VariableNames)
    T = groupsummary(T, {'Animal','Electrode','Condition','Frequency'}, 'mean', 'Power');
    T.Properties.VariableNames(strcmp(T.Properties.VariableNames,'mean_Power')) = {'Power'};
end

T.Frequency = round(T.Frequency, 3);

% Reorder Condition levels: Baseline first, then Stim* in numeric order if present
lev = categories(T.Condition);
iBase = find(strcmpi(lev,'baseline'), 1);
assert(~isempty(iBase), 'A "Baseline" level is required in Condition.');
stimLev = lev(~strcmpi(lev,'baseline'));
nums = regexp(stimLev, '\d+', 'match', 'once');
numv = cellfun(@(s) iff(isempty(s), Inf, str2double(s)), nums);
[~, idx] = sort(numv);
stimLev = stimLev(idx);
wantLev = [{'Baseline'}, stimLev(:)'];
T.Condition = reordercats(T.Condition, wantLev);

% Numeric condition code if absent
if ~ismember('CondNum', T.Properties.VariableNames)
    T.CondNum = double(T.Condition) - 1;
end

% Drop rows with missing key values
T = rmmissing(T, 'DataVariables', {'Power','CondNum'});

%% 2) Power transform (optional log)
p = double(T.Power);
if doLogPower
    if any(~isfinite(p))
        T = T(isfinite(p), :);
        p = double(T.Power);
    end
    minP = min(p);
    posP = p(p > 0);
    hasPos = ~isempty(posP);
    if minP <= 0
        if hasPos, delta = abs(minP) + 0.5*min(posP); else, delta = abs(minP) + 1; end
        T.PowerX = log(p + delta);
    else
        T.PowerX = log(p);
    end
    bad = ~isfinite(T.PowerX) | ~isreal(T.PowerX);
    if any(bad), T(bad, :) = []; end
else
    T.PowerX = p;
end

%% 3) Per-frequency LMM (no FDR)
freqs = unique(T.Frequency);
nF    = numel(freqs);
allRows = {};
overallP = nan(nF,1);

for i = 1:nF
    f  = freqs(i);
    Tf = T(T.Frequency==f, :);
    form = 'PowerX ~ Condition + (1|Animal) + (1|Animal:Electrode)';
    try
        lme = fitlme(Tf, form, 'FitMethod','REML');

        a = anova(lme, 'DFMethod','satterthwaite');
        r = find(strcmp(a.Term, 'Condition'));
        if ~isempty(r), overallP(i) = a.pValue(r); end

        coe = dataset2table(lme.Coefficients);  
        sel = startsWith(coe.Name, 'Condition_');
        if any(sel)
            sub = coe(sel, {'Name','Estimate','SE','DF','tStat','pValue'});
            sub.Frequency_Hz = repmat(f, height(sub),1);
            sub.Level = erase(string(sub.Name), 'Condition_');
            sub = movevars(sub, {'Frequency_Hz','Level'}, 'Before', 'Name');
            allRows{end+1} = sub; %#ok<SAGROW>
        end
    catch ME
        warning('LMM failed at %.3f Hz: %s', f, ME.message);
    end
end

if ~isempty(allRows)
    pairwiseTbl = vertcat(allRows{:});
else
    pairwiseTbl = table();
end

%% 5) Estimated means ±95% CI + SD 
emmRows = {};
for i = 1:nF
    f  = freqs(i);
    Tf = T(T.Frequency==f, :);
    try
        lme = fitlme(Tf, 'PowerX ~ Condition + (1|Animal) + (1|Animal:Electrode)', 'FitMethod','REML');
        conds = categorical(wantLev(:), wantLev);
        nC    = numel(conds);

        tmp = table();
        tmp.Animal    = repmat(mode(Tf.Animal),    nC, 1);
        tmp.Electrode = repmat(mode(Tf.Electrode), nC, 1);
        tmp.Condition = conds;
        tmp.CondNum   = (0:nC-1)';

        [yhat, ci] = predict(lme, tmp);
        if doLogPower, yhat = exp(yhat); ci = exp(ci); end

        sdVals = nan(nC,1);
        for c = 1:nC
            thisCond = wantLev{c};
            vals = Tf.Power(Tf.Condition==thisCond);
            if numel(vals)>1
                sdVals(c) = std(vals, 0, 'omitnan');
            else
                sdVals(c) = NaN;
            end
        end

        tt = table(repmat(f,nC,1), string(conds), yhat, ci(:,1), ci(:,2), sdVals, ...
                   'VariableNames', {'Frequency_Hz','Condition','Mean','Low95','High95','SD'});
        emmRows{end+1} = tt;
    catch ME
        warning('Frequency %.3f Hz: %s', f, ME.message);
    end
end
emmTbl = vertcat(emmRows{:});

% Save per-frequency results
outPair = fullfile(saveFolder, 'LMM_Power_Pairwise_byFreq.xlsx');
outEMM  = fullfile(saveFolder, 'LMM_Power_EMMeans_byFreq.xlsx');
writetable(pairwiseTbl, outPair, 'Sheet','Pairwise');
writetable(table(freqs, overallP, 'VariableNames',{'Frequency_Hz','ANOVA_p_Condition'}), ...
           outPair, 'Sheet','Overall_ANOVA', 'WriteMode','overwritesheet');
writetable(emmTbl, outEMM, 'Sheet','EMM');
fprintf('Saved:\n  %s\n  %s\n', outPair, outEMM);

%% 10) Spectra by condition (Mean ± SEM) — optional QC figure
sumTbl = grpstats(T, {'Condition','Frequency'}, {'mean','std','numel'}, 'DataVars','Power');
sumTbl.SEM = sumTbl.std_Power ./ sqrt(sumTbl.numel_Power);

figure('Color','w'); hold on;
conds = wantLev; nC = numel(conds);
for c = 1:nC
    C = conds{c};
    rowC = sumTbl.Condition==C;
    fC   = sumTbl.Frequency(rowC);
    [fC, ord] = sort(fC);
    mC   = sumTbl.mean_Power(rowC); mC = mC(ord);
    sC   = sumTbl.SEM(rowC);        sC = sC(ord);
    errorbar(fC, mC, sC, 'o-','LineWidth',1.2);
end
xlabel('Frequency (Hz)'); ylabel('Power (mean ± SEM)');
title('Power spectra by condition'); grid on; box on; legend(conds, 'Location','best');

% Mark significant frequencies (raw p < alpha, no FDR)
if ~isempty(allRows)
    yl = ylim;
    for c = 2:nC
        lvl = conds{c};
        sub = vertcat(allRows{:});
        sel = strcmp(string(sub.Level), lvl);
        if any(sel)
            sub = sub(sel,:);
            sigF = sub.Frequency_Hz(sub.pValue < alphaPV);
            for k = 1:numel(sigF)
                f0 = sigF(k);
                r  = sumTbl.Condition==lvl & sumTbl.Frequency==f0;
                if any(r), y0 = sumTbl.mean_Power(r);
                else, y0 = yl(2)*0.95; end
                text(f0, y0, '*', 'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', 'FontSize',12);
            end
        end
    end
    ylim(yl);
end

%% 11) Band averages + LMM (Stim vs Baseline) — for band power
bands = [0 2; 2 5; 5 10];
bandNames = arrayfun(@(i)sprintf('%g–%g Hz', bands(i,1), bands(i,2)), (1:size(bands,1))', 'uni',0);
rows = [];
bandStatsRows = {};

figure('Color','w'); 
tiledlayout(1, size(bands,1), 'Padding','compact','TileSpacing','compact');

for b = 1:size(bands,1)
    f1 = bands(b,1); f2 = bands(b,2);
    idx = T.Frequency>=f1 & T.Frequency<f2;
    Tb = T(idx,:);
    if isempty(Tb), continue; end

    hasPX = ismember('PowerX', Tb.Properties.VariableNames);
    respVar = 'Power'; if hasPX, respVar = 'PowerX'; end
    form = sprintf('%s ~ Condition + (1|Animal) + (1|Animal:Electrode)', respVar);

    try
        lmeB = fitlme(Tb, form, 'FitMethod','REML');
        coeB = dataset2table(lmeB.Coefficients);
        selB = startsWith(coeB.Name,'Condition_');
        subB = coeB(selB, {'Name','Estimate','SE','DF','tStat','pValue'});
        subB.Band = repmat(string(bandNames{b}), height(subB),1);
        subB.Level = erase(string(subB.Name),'Condition_');
        rows = [rows; subB]; %#ok<AGROW>

        Sb = grpstats(Tb, 'Condition', {'mean','std','numel'}, 'DataVars','Power');
        Sb.SEM = Sb.std_Power ./ sqrt(Sb.numel_Power);

        nexttile; 
        bar(categorical(Sb.Condition, wantLev), Sb.mean_Power, 'FaceAlpha',0.7); hold on;
        errorbar(categorical(Sb.Condition, wantLev), Sb.mean_Power, Sb.SEM, ...
                 'k','LineStyle','none','LineWidth',1.2);
        title(bandNames{b}); ylabel('Power (mean ± SEM)'); box on; grid on;

        for c = 2:numel(wantLev)
            lvl = wantLev{c};
            pRow = rows(string(rows.Band)==bandNames{b} & rows.Level==lvl, :);
            if ~isempty(pRow)
                p = pRow.pValue(1);
                if p<0.001, star='***';
                elseif p<0.01, star='**';
                elseif p<0.05, star='*';
                else, star='ns';
                end
                y = Sb.mean_Power(Sb.Condition==lvl) + Sb.SEM(Sb.Condition==lvl);
                text(c, y*1.05, star, 'HorizontalAlignment','center','FontSize',11);
            end
        end
    catch ME
        warning('Band %s LMM failed: %s', bandNames{b}, ME.message);
    end
end

outBands = fullfile(saveFolder, 'LMM_Power_Band_Average.xlsx');
if ~isempty(rows)
    rows = rows(:, {'Band','Level','Estimate','SE','DF','tStat','pValue'});
    writetable(rows, outBands, 'Sheet','Band_LMM');
end
if ~isempty(bandStatsRows)
    bandStatsTbl = vertcat(bandStatsRows{:});
    writetable(bandStatsTbl, outBands, 'Sheet','Band_Summary');
end
fprintf('Band-level saved to:\n%s\n', outBands);

%% ===== Local functions =====
function o = iff(tf, a, b)
    if tf, o = a; else, o = b; end
end
