%% Preprocessing pipeline for early iTEPs in single pulse and tripulse data

% DRCMR 2025 - XCT & LS - BrainStim Methods group

% Preprocessing pipeline OPTION 1 - BASIC


% --> GENERAL INFO:
% The current preprocessing pieplpile is intended to provide a simple and fast eeg
% preprocessing focusing on the essential steps. The preprocessing might be
% sufficient to explore iTEPS in datasets non contaminated with large
% muscular (or other type of) artefacts.


% Only edit Section 0 to set the dataset you want to load. 
% All other sections are preconfigured—run them block by block; 
% they will progress semi-automatically and require no further changes.


% ---> STEPS:
% 0) Load data
% 1) Epoch
% 2) Remove bad trials 

% ----------------------> If we want to conserve the TMS pulse in the data
% 2.1) Baseline, filtering, demean and plot. Save data

% ------------------ ----> If we want to remove the Tms pule from the data (recommended),
% continue here
% 3) Demean data
% 4) Remove TMS arifact
% 5) Remove bad channels
% 6) Interpolate TMS removed data
% 7) Baseline correction
% 8) ICA: remove decay, osclar and muscular contractions only.
% 9) Final Baseline correction, demean and filtering [0.2 2000]
% 10) Save preprocessed dataset. The final datasethas the TMS pulse removed,% is baseline corrected, referenced to original mastoids and has not been filtered.




% --> EXTERNAL TOOLS TO INSTALL
% EEGlab and Fieltrip are NOT distributed with this code. Install them
% on your system and add the paths to matlab.

% EEGLAB install:  https://eeglab.org/tutorials/01_Install/Install.html
% -- Once EEglab installed, you will need to download the TESA plugin. From
% the eeglab interface go to : FILE-->MANAGE EEGLAB EXTENSIONS--> SEARCH TESA TOOLBOX AND INSTALL.

% FIELDTRIP install: https://www.fieldtriptoolbox.org/download/
% Fieltrip in only used for future analyses and not for the present
% preprocessing. IThere is no need to install if you are only runing this
% preprocessing code.


%% ADDPATHS 
clear all;
close all;
clc

%path to external functions.
addpath(genpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/'))

%path eto eeglab
addpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/external_functions/eeglab2025.0.0/')

%path to fieltrip
addpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/external_functions/fieldtrip-20240110/')


%% 0. LOAD DATA

% Selec file:---------> 
eeglab;
 % --------------------------->  Modify subject condition to load
Subject = 'XEB' % XPB or XEB or XTC
Coil = 'coilB35' % coilB35 or coilB65
Target = 'SFG' % SFG or SPL
Orientation = 'LM' % LM or APPA
Intensity = '100' % 70 or 80 or 90 or 100 
MT_MSO = 'MSO'% RMT or MSO
Paradigm = 'tripulse' % singlepulse or tripulse
% <--------------------------------------------------------------




% ---- build dataset name and path ----
name_dataset = sprintf('%s_%s_%s_%s_%s%s_%s.vhdr', ...
    Subject, Coil, Target, Orientation, Intensity, MT_MSO, Paradigm);

path_dataset = fullfile('/mnt/projects/PARADISE/PARADISE_1', Subject, 'EEG_clean');
filePath     = fullfile(path_dataset, name_dataset);

name = [Subject,' ',Coil,' ',Target,' ',Orientation,' ',Intensity,MT_MSO,' ',Paradigm];
%name_dataset = [Subject,'_',Coil,'_',Target,'_',Orientation,'_',Intensity,MT_MSO,'_',Paradigm,'_',pulse,'_cleaned_pipeline_',processing_pipeline,'.set']; 


% ---- sanity checks + conditional load ----
if ~isfolder(path_dataset)
    error('Folder not found: %s', path_dataset);
end

if isfile(filePath)
    fprintf('Loading: %s\n', filePath);
    EEG = pop_loadbv(path_dataset, name_dataset);
    eeglab redraw
else
    warning('File does not exist: %s', filePath);
    % Show helpful candidates in the same folder
    cand = dir(fullfile(path_dataset, sprintf('%s_*.vhdr', Subject)));
    if isempty(cand), cand = dir(fullfile(path_dataset, '*.vhdr')); end
    if isempty(cand)
        fprintf('No .vhdr files in %s\n', path_dataset);
    else
        fprintf('Available .vhdr files in %s:\n', path_dataset);
        for k = 1:min(numel(cand), 20)
            fprintf('  - %s\n', cand(k).name);
        end
        if numel(cand) > 20
            fprintf('  ... (%d more)\n', numel(cand)-20);
        end
    end
end

%% 1. EPOCH

%Epoch: assuming that onliny the first pulseof the tripulse burst is marked as event
EEG = pop_epoch( EEG, {  'R  8'  }, [-0.5         0.5]); % TMS pulses are stamped on the EEG as 'R  8' markers. Modify your marker if necessary for epoching.
eeglab redraw


% Figures
figure; pop_timtopo(EEG, [-5  10], [2.3         3         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 2. Remove bad trials
%% Trials rejection - (identify bad trials) --> Select the bad trials
TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'

 
%% Trials rejection - (remove bad trials)
if ~isempty(TMPREJ)
    [trialrej elecrej]=eegplot2trial(TMPREJ,size(EEG.data,2),size(EEG.data,3));
else
    trialrej=[];
end

EEG.BadTr =find(trialrej==1);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 

eeglab redraw



% Figures
figure; pop_timtopo(EEG, [-5  15], [2.3         3         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 15]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

EEG_trials_corrected = EEG;

%% 2.1 Baseline, filtering ( 0.2 5000Hz) demean. Later Save data with Pulse

% Data
EEGr = EEG;

% Baseline correct,Filter and  Demean
EEGr = pop_rmbase( EEGr, [-110 -10] ,[]);
%EEGr = pop_tesa_filtbutter(EEGr, 0.2, 5000, 2, 'bandpass');
EEG = pop_tesa_filtbutter( EEGr, [], 5000, 2, 'lowpass' ); %zero-phase, 4th-order low pass butterworth filter allowing frequencies below 25000 Hz
EEGr = pop_rmbase( EEGr, [-500 499.9] ,[]); % De meaning assuming a epoch from -500 499.9

% Mean of Channels of itnerest
% Figure
figure;
set(gcf, 'Position', [100, 100, 2000, 600]); % Manually adjust figure size [x y width height]


% Plot individual sensors (mean over trials), in grey
h1 = plot(EEGr.times, mean(EEGr.data, 3)', 'Color', [0.5 0.5 0.5]); 

hold on
% Plot mean of selected ROI channels, in red
%h2 = plot(EEGr.times, mean(mean_roi_raw, 3)', 'r'); 

% Axes settings
xlim([-5 20]);
ylim([-80 80]);
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('Raw EEG Signal 0.2-5000Hz', 'FontSize', 24, 'FontName', 'Arial');

% Legend
% Create dummy lines for legend
dummy1 = plot(NaN, NaN, 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5);
%dummy2 = plot(NaN, NaN, 'r', 'LineWidth', 2.5);

% Add legend using dummies
info = ['All sensors ',name]
legend([dummy1], {info}, 'FontSize', 18);
set(gca,'FontName','Arial','fontsize',20,'FontWeight','bold','LineWidth',1.5)

hold off;




%{
% In case you want topographies: modify data file to plot accordingly:
% Define the time points you want to plot (in ms)
time_points = [2.3, 3.5, 4.8]; % Example time points

% Number of time points
num_plots = length(time_points);

% Create a figure
figure;
set(gcf, 'Position', [100, 100, 2000, 900]); % Adjust figure size to fit both plots

% Create a subplot for topoplots (this will take up the top half of the figure)
subplot(2, 1, 1); % 2 rows, 1 column (top row)
% Loop over each time point for topography plots
for i = 1:num_plots
    time_point = time_points(i); % Get the current time point in ms

    % Find the nearest time index
    [~, time_idx] = min(abs(EEG.times - time_point)); % Find the closest index
    
    % Create a subplot for the current time point (topographies on the left)
    % We're using `subplot(1, num_plots, i)` for all topoplots together in one big subplot.
    subplot(1, num_plots, i); % In the 1 row, 'num_plots' columns (placing topoplots in one row)
    
    % Check if the selected index is valid and there is data at that time
    if ~isempty(time_idx) && ~all(isnan(EEG.data(:, time_idx))) 
        % Plot topography
        topoplot(EEG.data(:, time_idx), EEG.chanlocs, 'electrodes', 'on');
        title([num2str(EEG.times(time_idx)) ' ms'],'FontSize', 15); % Title with the actual time point
        % Add colorbar with custom settings
        colorbar; % Add colorbar
        caxis([-20 20]); % Set color limits (adjust based on your data range)
    else
        % If no valid data, plot a message
        topoplot([], EEG.chanlocs, 'electrodes', 'off'); % Empty topography
        title('No Data'); % Indicate no data for this time point
    end
end
%}

%% SAVE DATASET With TMS pulse?  Then execute this cell to save
%save_name
name_save = [Subject,'_',Coil,'_',Target,'_',Orientation,'_',Intensity,MT_MSO,'_',Paradigm];
%fname = strjoin(name_save, '');  % joins: "XEB_coilB65_SFG_APPA_70RMT_tripulse"

% Add your suffix + extension
fname = name_save + "_with_pulse_cleaned_pipeline_basic.set";  % still a string scalar

% Save (cast to char for EEGLAB)
EEG = pop_saveset(EEG, ...
    'filename', char(fname), ...
    'filepath', char(path_dataset));


% Finish here your preprocessing if you want the TMS pulse on your data. If
% you want to remove it, continue.

%% 3. De mean the data - substract the mean value of the data from each data point to remove DC offset and prepare data for filtering)

EEG = pop_rmbase( EEG, [-500 499.9] ,[]); % De meaning assuming a epoch from -500 499.9


% Figures
figure; pop_timtopo(EEG, [-5  15], [2.3         3         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 15]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 4. REMOVE TMS ARTIFACT
% TMS artefact lengths in ms

Paradigm = string(Paradigm);  % ensure string scalar

if Paradigm == "singlepulse"
    TMSremoval_range = [-0.5 2];
    EEG = pop_tesa_removedata(EEG, TMSremoval_range);

elseif Paradigm == "tripulse"
    TMSremoval_range1 = [-0.5 6.4];
    EEG = pop_tesa_removedata(EEG, TMSremoval_range1);
end
eeglab redraw


% Figures
figure; pop_timtopo(EEG, [-5  15],[], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 15]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 5. REMOVE BAD CHANNELS

% To identify bad channels we run the following section applying a Wiener filter to automatically detect bad channels:
EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')

%%
% labeling the very worst channels to not affect the ICA run
badC = find(sigmas > (median(sigmas) + 5*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);

EEG = pop_select( EEG, 'nochannel', badC);
badC
eeglab redraw

% Figures
figure; pop_timtopo(EEG, [-5  15], [2.3         3         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 15]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

% Run wiener again to verify the bad channel has been removed:
EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')


%{
% ---------------> Alternatively do it manually

% Indentify bad channel numbers
badChannels = [20];

% Remove channels from the data
EEG = pop_select( EEG, 'rmchannel',badChannels);
eeglab redraw

%}


%% 6. INTERPOLATE REMOVED WINDOW

if Paradigm == "singlepulse"
 EEG = pop_tesa_interpdata( EEG, ['cubic'], [abs(TMSremoval_range(1)) TMSremoval_range(2)] ); % Cubic interpolation. Alternatively do linear interpolation, or substitute tms segment with baseline data window

elseif Paradigm == "tripulse"
 EEG = pop_tesa_interpdata( EEG, ['cubic'], [abs(TMSremoval_range1(1)) TMSremoval_range1(2)] ); % Cubic interpolation. Alternatively do linear interpolation, or substitute tms segment with baseline data window

end
eeglab redraw


% Figures
figure; pop_timtopo(EEG, [-5  15], [ ], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 15]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


%% 7. BASELINE CORRECTION

EEG = pop_rmbase( EEG, [-110 -10] ,[]);

% Figures
figure; pop_timtopo(EEG, [-5  15], [2.3         3         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 15]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% Run ICA to remove ocular moments, muscular contractions and decay artifacts

% Resample for ica:
EEG = pop_resample( EEG, 5000); % maybe 25000? 5000 is enough for the iTEPS and faster

%% 8. Prepare data for ICA

EEG = pop_select(EEG, 'nochannel', 1);   % removes first channel and keeps EEG.data/chanlocs consistent

badC = find(sigmas > (median(sigmas) + 5*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);

EEG2ICA = pop_select( EEG, 'nochannel', badC);
badC

tmpdata=reshape(EEG2ICA.data,[size(EEG2ICA.data,1) size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]); 
tmpdata=tmpdata-repmat(mean(tmpdata,2),[1,size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]);
tempCpca=svd(tmpdata);
th=0.01;% this threshold is arbitrary but it works in most cases
figure,semilogy(tempCpca,'.-')% plot of the eigenvalues of the data matrix
a=ylim;
hold on
Cpca=max(find(tempCpca>th));
line([Cpca Cpca],[a(1) a(2)],'Color','k','LineStyle','--');
hold off
'check whether the vertical dashed line correctly separates the eigenvalues close to zero from the ones much higher than zero, otherwise change the threshold th'

% if the vertical dotted line correctly separates the eigenvalues close to
% zero from the ones much higher than zero, Cpca correctly corresponds to
% the total number of independent components that can be estimated from the
% data matrix. if this is not the case, you have to change the threshold above or to manually correct Cpca


%% 8. ICA - First round

% Set the number of dimentions accordingly
% You might need to go one below the previous value
% Here we only remove eye-blinks and horizontal movements: to make your
% choise, if unsure check: https://www-sciencedirect-com.proxy.insermbiblio.inist.fr/science/article/pii/S1053811916305845

EEG = pop_tesa_pcacompress( EEG, 'compVal', 55, 'plot','on' ); % Set the number of dimentions here
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
 

%% Plot the independent components following ICA. Here the user can verify the results of ICA
% As well as manually classify components. For details on classifying components, see
% https://nigelrogasch.gitbook.io/tesa-user-manual/remove_minimise_tms_muscle_activity/auto_comp_select
% Modify channel type name for runing sound

for i = 1:length(EEG.chanlocs)
EEG.chanlocs(i).type = 'EEG';
end 

% Pop automatic preselected components
EEG = pop_tesa_compselect( EEG,'plotTimeX',[-200,499],'muscleFreqIn',[7,70],'tmsMuscleWin',[-2,30],'blinkElecs',{'AF3', 'AF4'},'tmsMuscleThresh',8);


% Alternatively do it manually
% EEG = pop_tesa_compplot( EEG,'figSize','medium','plotTimeX',[-50 250],'plotFreqX',[1 500], 'freqScale','log', 'saveWeights','off');


%% Plot data
figure; pop_timtopo(EEG, [-5  15], [6.4], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 15]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


%% 9. BASELINE CORRECTIOn, demean and filtering
% baseline correction
EEG = pop_rmbase( EEG, [-110 -10] ,[]);
% demean
EEG = pop_rmbase( EEG, [-500 499] ,[]); % De meaning assuming a epoch from -500 499.9
%filtering
EEG = pop_tesa_filtbutter(EEG, 0.2, 2000, 2, 'bandpass');  %----->> Do we want this?

%% 9.1. PLOT FINAL DATA and inspect any error
% Plot data again, full data, high fqreuency data iTEP and low frqeuency TEP

figure; pop_timtopo(EEG, [-5  15], [5 6 6.6 7.4 8.4], 'ERP data and scalp maps 0.2-2000Hz');
% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 15]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('EEG Signal 0.2-2000Hz', 'FontSize', 24, 'FontName', 'Arial');
figure; pop_timtopo(EEG, [-50  150], [5 6 6.6 7.4 8.4], 'ERP data and scalp maps 0.2-2000Hz');


% Make some more plots isolating high freqeuncy (aka iTEPS) and low
%frqeuency (aka TEP) time series

% filtering to isolate High frqeuency iTEP time series
EEG2 = pop_tesa_filtbutter(EEG, 400, 1000, 2, 'bandpass');
figure; pop_timtopo(EEG2, [-5  15], [5 6 6.6 7.4 8.4], 'ERP data and scalp maps 400-1000Hz');
% Figure for paper
figure;
plot(EEG2.times, mean(EEG2.data, 3)', 'b'); 
xlim([-5 15]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('EEG Signal 400-1000Hz', 'FontSize', 24, 'FontName', 'Arial');
figure; pop_timtopo(EEG2, [-50  150], [5 6 6.6 7.4 8.4], 'ERP data and scalp maps 400-1000Hz');



% filtering to isolate Low frqeuency TEP time series
EEG3 = pop_tesa_filtbutter(EEG, 2, 100, 2, 'bandpass');
figure; pop_timtopo(EEG3, [-5  15], [5 6 6.6 7.4 8.4], 'ERP data and scalp maps 2-100Hz');
% Figure for paper
figure;
plot(EEG3.times, mean(EEG3.data, 3)', 'b'); 
xlim([-5 15]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('EEG Signal 2-100Hz', 'FontSize', 24, 'FontName', 'Arial');
figure; pop_timtopo(EEG3, [-50  150], [5 6 6.6 7.4 8.4], 'ERP data and scalp maps 2-100Hz');

%%  Add  finaldetrneding If necessary

%{
close all
% Detrend
EEG = pop_tesa_detrend(EEG, 'linear', [-2 499]);
%}

%% 9. SAVE DATASETS


%save name
name_save = [Subject,'_',Coil,'_',Target,'_',Orientation,'_',Intensity,MT_MSO,'_',Paradigm];
fname = strjoin(name_save, '');  % joins: "XEB_coilB65_SFG_APPA_70RMT_tripulse"

% Add your suffix + extension
fname = fname + "_no_pulse_cleaned_pipeline_basic.set";  % still a string scalar

% Save (cast to char for EEGLAB)
EEG = pop_saveset(EEG, ...
    'filename', char(fname), ...
    'filepath', char(path_dataset));


%% END
close all
clc
