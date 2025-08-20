%% Preprocessing pipeline for early iTEPs in single pulse and tripulse data

% DRCMR 2025 - XCT& LS  - BrainStim Methods group

% Preprocessing pipeline OPTION 2 - ADVANCED

% The preprocessing is inspired in the following references:
%https://www.sciencedirect.com/science/article/pii/S2666166725000280?via%3Dihub
%https://www.sciencedirect.com/science/article/pii/S1053811921005486

% --> GENERAL INFO:
% The current preprocessing pieplpile is intended to provide a elaborated and meticulous eeg
% preprocessing focusing on removing backgrond noise, TMS induced artifacts, muscular artifacts and eye movements. The preprocessing might be
% sufficient to explore iTEPS in datasets contaminated with muscular (or other type of) artifacts.

% Only edit Section 0 to set the dataset you want to load. 
% All other sections are preconfigured—run them block by block; 
% they will progress semi-automatically and require no further changes.

% ---> STEPS:
% 0) Load data
% 1) Epoch
% 2) De-mean data
% 3) Remove 50Hz noise
% 4) Remove bad trials
% 5) Remove TMS artifact
% 6) Interpolate removed TMS signal
% 7) Remove recharging artefact
% 8) First roud ICA - remove only eyes blinks and movements
% 9) High pass filter from 2Hz
% 10) Baseline correct the data
% 11) Remove recording noise and interpolate bad electrodes with SOUND
% 12) De-mean data
% 13) Supres time-locked muscular artifacts with SSP-SIR
% 14) Second round of ICA - remove persistent residual artifacts
% 15) Remove noisiest trials
% 16) Bandpass filter, baseline correction
% 17) Plot final data
% 18) Save datasets



%% ADDPATHS 
clear all;
close all;
clc

% Path to external functions
addpath(genpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/'))

% Path to eeglab
addpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/external_functions/eeglab2025.0.0/')

% Path to fieltrip
addpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/external_functions/fieldtrip-20240110/')

%
%% 0. LOAD DATA

% Open EEGlab 
eeglab;
 % --------------------------->  Modify subject condition to load
Subject = 'XEB' % XPB or XEB or XTC
Coil = 'coilB65' % coilB35 or coilB65
Target = 'SFG' % SFG or SPL
Orientation = 'APPA' % LM or APPA
Intensity = '70' % 70 or 80 or 90 or 100 
MT_MSO = 'RMT'% RMT or MSO
Paradigm = 'tripulse' % singlepulse or tripulse
% <-----------------------------





% Dataset name structure ----->  [Subject]_[Coil]_[Target]_[Orientation]_[Intensity]_[Paradigm].ext
% Define file to load 
name = [Subject,' ',Coil,' ',Target,' ',Orientation,' ',Intensity,MT_MSO,' ',Paradigm];
name_dataset = [Subject,'_',Coil,'_',Target,'_',Orientation,'_',Intensity,MT_MSO,'_',Paradigm,'.vhdr']; 
path_dataset = ['/mnt/projects/PARADISE/PARADISE_1/',Subject,'/EEG_clean/'];

%load file
EEG = pop_loadbv(path_dataset, name_dataset);
eeglab redraw


%% 1. EPOCH

% Epoch: assuming that onliny the first pulseof the tripulse burst is
% marked as event

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

%% 2. DE MEAN

EEG = pop_rmbase( EEG, [-500  499.9]);

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 3. REMOVE LINE NOISE 50hZ FITTING A 50HZ OSCILLATION

EEG.data = fit_and_remove_line_noise_from_trials(EEG.data, EEG.srate, 50);

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 4. Remove bad trials
%% Trials rejection - 4.1 (identify bad trials) --> Select the bad trials

TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'

 
%% Trials rejection - 4.2 (remove bad trials)

if ~isempty(TMPREJ)
    [trialrej elecrej]=eegplot2trial(TMPREJ,size(EEG.data,2),size(EEG.data,3));
else
    trialrej=[];
end

EEG.BadTr =find(trialrej==1);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 

eeglab redraw


%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


%% 5. REMOVE TMS ARTIFACT

Paradigm = string(Paradigm);  % ensure string scalar

if Paradigm == "singlepulse"
    TMSremoval_range = [-0.5 1.5];
    EEG = pop_tesa_removedata(EEG, TMSremoval_range);

elseif Paradigm == "tripulse"
    TMSremoval_range1 = [-0.5 5.9];
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

%% 6. TMS Interpolation

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

%% 7. Recharging artefact if necessary

EEG = EEGLAB_remove_and_interpolate_recharging_artifact(EEG, 1); % 1 + single pulse, 2 = repetitive pulse

%% 8.*******First round of ICA --> identify oscular movements


% Resample for ica:
EEG = pop_resample( EEG, 5000); %25000 if necessary
eeglab redraw


%% Prepare the data for ICA: First remove for a while the very worst channels!

EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')


%% Labeling the very worst channels to not affect the ICA run

badC = find(sigmas > (median(sigmas) + 5*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);

EEG2ICA = pop_select( EEG, 'nochannel', badC);
badC


%% Before applying ICA to reduce artifacts, check which is the maximum number of independent components via singular value decomposition


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


%% ICA - First round

% Set the number of dimentions accordingly
% You might need to go one below the previous value
% Here we only remove eye-blinks and horizontal movements: to make your
% choise, if unsure check: https://www-sciencedirect-com.proxy.insermbiblio.inist.fr/science/article/pii/S1053811916305845

EEG = pop_tesa_pcacompress( EEG, 'compVal', 60, 'plot','on' ); % Set the number of dimentions here
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
 

%% Plot the independent components following ICA. Here the user can verify the results of ICA
% as well as manually classify components. For details on classifying components, see
% https://nigelrogasch.gitbook.io/tesa-user-manual/remove_minimise_tms_muscle_activity/auto_comp_select
EEG = pop_tesa_compplot( EEG,'figSize','medium','plotTimeX',...
    [-50 250],'plotFreqX',[1 500], 'freqScale','log', 'saveWeights','off');


%% Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 9. 2Hz High pass filter

EEG = pop_tesa_filtbutter( EEG, 0.2, [], 4, 'highpass' );

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 10. Baseline correction

EEG = pop_rmbase( EEG, [-110 -10] ,[]);

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


%% 11. Deal with channel noise if necessary with SOUND
%make a copy
data_before_sound = EEG;

% Modify channel type name for runing sound
for i = 1:length(EEG.chanlocs)
EEG.chanlocs(i).type = 'EEG';
end 

EEG = pop_select(EEG, 'nochannel', 1);   % removes first channel and keeps EEG.data/chanlocs consistent


% OPTION 1: USE SOUND
EEG = pop_tesa_sound(EEG, 'lambdaValue', 0.1, 'iter', 10 ); % spherical 3-layer reference model


%% Plot data
figure; pop_timtopo(data_before_sound, [-5  15], [2.3         4.5         4.8], 'ERP data and scalp maps BEFORE SOUND');
figure; pop_timtopo(EEG, [-5  15], [2.3         4.5         4.8], 'ERP data and scalp maps AFTER SOUND');

% Figure for paper
figure;
plot(data_before_sound.times, mean(data_before_sound.data, 3)', 'k'); 
hold on
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('EEG Signal BEFORE SOUND-BLACK  -  AFTER SOUND-BLUE', 'FontSize', 24, 'FontName', 'Arial');

% Figure for paper
figure;
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

plot(data_before_sound.times, mean(data_before_sound.data, 3)', 'r'); 
hold on
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('EEG Signal BEFORE SOUND-BLACK  -  AFTER SOUND-BLUE', 'FontSize', 24, 'FontName', 'Arial');
%legend({'Before SOUND', 'After SOUND'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');
hold off


% Figure for paper
figure;
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

plot(data_before_sound.times, mean(data_before_sound.data, 3)', 'k'); 
hold on
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
%xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('EEG Signal BEFORE SOUND-BLACK  -  AFTER SOUND-BLUE', 'FontSize', 24, 'FontName', 'Arial');
%legend({'Before SOUND', 'After SOUND'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');
hold off



%
%  OPTION 2: INTERPOLATE IT WITH EEGLAB. In the GUI go to TOOLS-->INTERPOLATE
% Interpolate electrode : check manually previous outsanding electrode and change it
%{
EEG = pop_interp(EEG, [7], 'spherical'); % Example to interpolate electrode numb 7

clear data

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%}

%% 12. De-mean the data

EEG = pop_rmbase( EEG, [-500 499] ,[]);

%% 13. SSP-SiR for muscular evoked activity

%make a copy
data_before_sspsir = EEG;

%
[EEG] = pop_tesa_sspsir(EEG, 'artScale', 'automatic');

%% Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(data_before_sspsir.times, mean(data_before_sspsir.data, 3)', 'k'); 
hold on
plot(EEG.times, mean(EEG.data, 3)', 'b');xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


% Figure for paper
figure;
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

plot(data_before_sspsir.times, mean(data_before_sspsir.data, 3)', 'k'); 
hold on
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');
%legend({'Before SSP-SIR', 'After SSP-SIR'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');
hold off


% Figure for paper
figure;
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

plot(data_before_sspsir.times, mean(data_before_sspsir.data, 3)', 'k'); 
hold on
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
%xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');
%legend({'Before SSP-SIR', 'After SSP-SIR'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');
hold off



% Figure for paper
figure;
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

plot(data_before_sspsir.times, mean(data_before_sspsir.data, 3)', 'r'); 
hold on
plot(data_before_sspsir.times, mean(data_before_sspsir.data, 3)', 'k'); 

plot(EEG.times, mean(EEG.data, 3)', 'b'); 
%xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');
%legend({'Before SSP-SIR','Before SSP-SIR', 'After SSP-SIR'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');
hold off


%% 14.*******Second round of ICA --> indentify residual muscular and persisting artifacts
%% Prepare the data for ICA: First remove for a while the very worst channels!

EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')


%% Labeling the very worst channels to not affect the ICA run

badC = find(sigmas > (median(sigmas) + 5*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);

EEG2ICA = pop_select( EEG, 'nochannel', badC);
badC


%% Before applying ICA to reduce artifacts, check which is the maximum number of independent components via singular value decomposition

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


%% ICA - second round

% Set the number of dimentions accordingly
% You might need to go one below the previous value
% Here we  remove persistent muscular artifacts and noise
% choise, if unsure check: https://www-sciencedirect-com.proxy.insermbiblio.inist.fr/science/article/pii/S1053811916305845

EEG = pop_tesa_pcacompress( EEG, 'compVal', 55, 'plot','on' ); % Set the number of dimentions here
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
 

%% Plot the independent components following ICA. Here the user can verify the results of ICA
% as well as manually classify components. For details on classifying components, see
% https://nigelrogasch.gitbook.io/tesa-user-manual/remove_minimise_tms_muscle_activity/auto_comp_select
EEG = pop_tesa_compplot( EEG,'figSize','medium','plotTimeX',...
    [-500 499],'plotFreqX',[1 100], 'freqScale','log', 'saveWeights','off' );



%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');



%% 15. Remove noisiest trials
%% Trials rejection - 15.1 (identify bad trials) --> Select the bad trials
TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'

 
%% Trials rejection - 15.2 (remove bad trials)
if ~isempty(TMPREJ)
    [trialrej elecrej]=eegplot2trial(TMPREJ,size(EEG.data,2),size(EEG.data,3));
else
    trialrej=[];
end

EEG.BadTr =find(trialrej==1);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 

eeglab redraw



% Figures
figure; pop_timtopo(EEGall, [-50  300], [3         4         5 9 19 30 50 80 130 220 260]);

figure; pop_timtopo(EEGall, [-5  15], [3         4         6  9], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 16.BANDPASS FILTER  rereferencing and baseline correction
% baseline correction
EEG = pop_rmbase( EEG, [-110 -10] ,[]);
% demean
EEG = pop_rmbase( EEG, [-500 499] ,[]); % De meaning assuming a epoch from -500 499.9
%filtering
EEG = pop_tesa_filtbutter(EEG, 0.2, 2000, 2, 'bandpass');



%% 17. PLOT FINAL DATA. SEGMENT FAST OSCILLATION around iTEP and plot

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
title('EEG Signal 2-100Hz', 'FontSize', 24, 'FontName', 'Arial');
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


%% 18. SAVE DATASETS
EEG = pop_saveset( EEG, 'filename',[name_dataset,'_no_pulse_cleaned_pipeline_advanced.set'],'filepath',[path_dataset]);


%% END
