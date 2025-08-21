%% Analyses pipeline for early iTEPs in single pulse and tripulse data

% DRCMR 2025 - XCT& LS  - BrainStim Methods group

% Analyses pipeline 

% Only edit Section 0 to set the dataset you want to load. 
% All other sections are preconfigured—run them block by block; 
% they will progress semi-automatically and require no further changes.


% --> GENERAL INFO:
% The current analyses pieplpile is intended to extract from a PREPORCESSED DATSET the following results:
% 1. iTEPS (aka. time series signal grand average across trials in time windows close to TMS pulse
% 1.1. TEPS (aka. time series signal grand average across trials in the entire epoch window
% 2. GMFP (aka. standand deviation across all sensors in the time domain. LMPF can also be computed by just adjusting the electrodes of interest
% 3. TF (aka. time frequency power spectra


% ---> STEPS:
% 0) Load the preprocessed data (EEGLAB preporcessed file)
% 1) iTEP and TEP analsyes
% 2) GMFP analyses
% 3) ERSP time freqeuncy


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



%% 0. Load dataset on interest to analyse.


% Open EEGlab 
 % --------------------------->  Modify subject condition to load
Subject = 'XEB' % XPB or XEB or XTC
Coil = 'coilB65' % coilB35 or coilB65
Target = 'SFG' % SFG or SPL
Orientation = 'APPA' % LM or APPA
Intensity = '90' % 70 or 80 or 90 or 100 
MT_MSO = 'RMT'% RMT or MSO
Paradigm = 'tripulse' % singlepulse or tripulse
pulse = 'no_pulse'; % with_pulse  or  no_pulse
processing_pipeline = 'basic'  % basic   or advanced
% <-----------------------------



% Dataset name structure ----->  [Subject]_[Coil]_[Target]_[Orientation]_[Intensity]_[Paradigm].ext
% Define file to load 
name = [Subject,' ',Coil,' ',Target,' ',Orientation,' ',Intensity,MT_MSO,' ',Paradigm];
name_dataset = [Subject,'_',Coil,'_',Target,'_',Orientation,'_',Intensity,MT_MSO,'_',Paradigm,'_',pulse,'_cleaned_pipeline_',processing_pipeline,'.set']; 
path_dataset = ['/mnt/projects/PARADISE/PARADISE_1/',Subject,'/EEG_clean/'];


% LOAD file with fieltrip
cfg = []; 
     cfg.lpfilter        = 'no';
     cfg.dataset = [path_dataset, name_dataset];
     data = ft_preprocessing(cfg);


%% 1------------------> iTEP and TEP analyses and plots (full frqeuency 0.2 - 2000Hz)
% Set all figures inside matlab, no pop up windows
%  set(0,'DefaultFigureWindowStyle','docked')  % set figure windows inside matlab
%  set(0,'DefaultFigureWindowStyle','remove')  % set figure windows outside matlab
%

cfg = [];
cfg.channel = {'all'};  % alternatively select the channels of interet. Example cf ro channels around left M1:  cfg.channel = {'C1','C3','CP1','CP3'} ; 
cfg.latency  = [-0.2 0.4];
cfg.keeptrials = 'no'; 
TEP = ft_timelockanalysis(cfg,data);

cfg = [];
cfg.baseline = [-0.2 -.05];
cfg.parameter = 'avg';
cfg.baselinetype = 'absolute';
TEP_bc = ft_timelockbaseline(cfg, TEP); %bc = baseline corrected


cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
TEP_bc_gav       = ft_timelockgrandaverage(cfg, TEP_bc);  


cfg = [];
cfg.channel    = 'all';
cfg.bpfilttype = 'firws';
cfg.demean     = 'yes';
cfg.baselinewindow = [-0.2 -0.05] % in seconds, the default is the complete trial (default = 'all')
cfg.bsfreq        = [48 52] %bandstop frequency range, specified as [low high]
TEP_bc_gav = ft_preprocessing(cfg, TEP_bc_gav);


% Inspect TEPs
cfg = [];
cfg.layout = 'biosemi64.lay'; %  !!!!!!!!!! Doble check layout, not sure !!!!!!!!!!!
cfg.showcomment   =  'no';
%cfg.ylim = [-5 5];
cfg.xlim = [-0.05 0.4]
figure; ft_multiplotER(cfg, TEP_bc_gav)

%%
pulse = strtrim(string(pulse));  % normalize to a string scalar
if pulse == "no_pulse" 
        
           % --- pick your original window (same as your code) ---
        idx = 951:1076;                           % sample indices you were plotting
        t_ms  = TEP_bc_gav.time(idx) * 1000;      % time in ms (from FieldTrip: seconds -> ms)
        seg   = TEP_bc_gav.avg(:, idx);           % [channels x N] segment
        
        % --- original and target sampling rates ---
        Fs_old = 1 / median(diff(TEP_bc_gav.time));   % Hz, from your time vector
        Fs_new = 25000;                               % target Hz
        [p,q]  = rat(Fs_new / Fs_old);                % 25000/2000 = 25/2
        
        % --- resample (polyphase interpolation with anti-imaging filter) ---
        seg_up = resample(seg.', p, q).';             % transpose so resample works along time
        
        % --- new time axis (match window endpoints) ---
        t_up_ms = linspace(t_ms(1), t_ms(end), size(seg_up, 2));
        
        % --- plot ---
        figure('Units','normalized','Position',[0.1 0.1 0.4 0.4])
        hold on
        for ch = 1:size(seg_up,1)
            plot(t_up_ms, seg_up(ch,:), 'LineWidth', 2);
        end
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        xlim([-5 15])
        ylim([-15 15])
        title('iTEP EEG Signal 0.2–2000Hz', 'FontSize', 24, 'FontName', 'Arial');
        set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)
        
        % grey patch spanning -0.5 to 6.4 ms
        if Paradigm == 'singlepulse' % singlepulse or tripulse
        
            yl = ylim;
            patch([-0.5 -0.5 2 2], [yl(1) yl(2) yl(2) yl(1)], ...
                  [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
        
        elseif  Paradigm == 'tripulse'
            yl = ylim;
            patch([-0.5 -0.5 6.4 6.4], [yl(1) yl(2) yl(2) yl(1)], ...
                  [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
        end 
        
        hold off
        
        
        
        
        % PLOT FULL TEP
        figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]) 
        % [left bottom width height]
        
        time = TEP_bc_gav.time(251:3001) * 1000; % convert to ms
        
        for i = 1:size(TEP_bc_gav.avg,1)
            plot(time, TEP_bc_gav.avg(i,251:3001),LineWidth=2);
            hold on
        end
        
        % --- add grey rectangle ---
        
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        xlim([-20 400])
        ylim([-15 15])
        
        title('TEP EEG Signal 0.2-2000Hz', 'FontSize', 24, 'FontName', 'Arial');
        set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)
        
        % add grey patch spanning -0.2 to 6.4 ms across whole y-range
        if Paradigm == 'singlepulse' % singlepulse or tripulse
        
            yl = ylim;
            patch([-0.5 -0.5 2 2], [yl(1) yl(2) yl(2) yl(1)], ...
                  [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
        
        elseif  Paradigm == 'tripulse'
            yl = ylim;
            patch([-0.5 -0.5 6.4 6.4], [yl(1) yl(2) yl(2) yl(1)], ...
                  [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
        end 



elseif pulse == "with_pulse" % with_pulse  or  no_pulse

        % PLOT FULL TEP with time in ms on x-axis FOR datasets WITH TMS pulse
        t_ms = TEP_bc_gav.time * 1000;                  % seconds -> ms
        if iscolumn(t_ms), t_ms = t_ms.'; end           % make row to match avg
        
        figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]); hold on
        for i = 1:size(TEP_bc_gav.avg,1)
            plot(t_ms, TEP_bc_gav.avg(i,:), 'LineWidth', 2);
        end
        
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        ylim([-100 100])
        xlim([t_ms(1) t_ms(end)])                       % fit to full time range
        title('TEP EEG Signal 0.2–2000 Hz', 'FontSize', 24, 'FontName', 'Arial');
        set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)

end 

%% The same code but only for high frqeuency data ( aka iTEP time series data plots 400-1000Hz)

% --- BAND-PASS ---
cfg = [];
cfg.channel    = 'all';
cfg.bpfilter   = 'yes';
cfg.bpfreq     = [400 1000];       % passband in Hz
%cfg.bpfilttype = 'firws';
cfg.demean     = 'yes';
data_bp = ft_preprocessing(cfg, data);


cfg = [];
cfg.channel = {'all'};  % alternatively select the channels of interet. Example cf ro channels around left M1:  cfg.channel = {'C1','C3','CP1','CP3'} ; 
cfg.latency  = [-0.2 0.4];
cfg.keeptrials = 'no'; 
TEP = ft_timelockanalysis(cfg,data_bp);

cfg = [];
cfg.baseline = [-0.2 -.05];
cfg.parameter = 'avg';
cfg.baselinetype = 'absolute';
TEP_bc = ft_timelockbaseline(cfg, TEP); %bc = baseline corrected


cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
TEP_bc_gav       = ft_timelockgrandaverage(cfg, TEP_bc);  


cfg = [];
cfg.channel    = 'all';
%cfg.bpfilttype = 'firws';
cfg.demean     = 'yes';
TEP_bc_gav = ft_preprocessing(cfg, TEP_bc_gav);


% Inspect TEPs
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.showcomment   =  'no';
cfg.xlim = [-0.05 0.4]
figure; ft_multiplotER(cfg, TEP_bc_gav)


if pulse == "no_pulse" % with_pulse  or  no_pulse

        % --- select the original window ---
        idx   = 951:1076;                         % samples you were plotting
        t_ms  = TEP_bc_gav.time(idx) * 1000;      % ms
        seg   = TEP_bc_gav.avg(:, idx);           % [channels x N]
        
        % --- original and target sampling rates ---
        Fs_old = 1 / median(diff(TEP_bc_gav.time));  % Hz (from time vector in seconds)
        Fs_new = 25000;                               % target Hz
        [p,q]  = rat(Fs_new / Fs_old);                % 25000/2000 = 25/2
        
        % --- upsample (polyphase, with anti-imaging filter) ---
        seg_up = resample(seg.', p, q).';             % resample along time (transpose in/out)
        
        % --- new time axis (match original window endpoints) ---
        t_up_ms = linspace(t_ms(1), t_ms(end), size(seg_up, 2));
        
        % --- plot iTEP (upsampled) ---
        figure('Units','normalized','Position',[0.1 0.1 0.4 0.4])
        hold on
        for i = 1:size(seg_up,1)
            plot(t_up_ms, seg_up(i,:));               % high-res traces
        end
        
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        xlim([-5 15])
        ylim([-15 15])
        title('iTEP EEG Signal 400–1000Hz', 'FontSize', 24, 'FontName', 'Arial');
        set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)
        
        % grey patch spanning -0.5 to 5.8 ms across whole y-range
        if Paradigm == 'singlepulse' % singlepulse or tripulse
        
            yl = ylim;
            patch([-0.5 -0.5 2 2], [yl(1) yl(2) yl(2) yl(1)], ...
                  [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
        
        elseif  Paradigm == 'tripulse'
            yl = ylim;
            patch([-0.5 -0.5 6.4 6.4], [yl(1) yl(2) yl(2) yl(1)], ...
                  [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
        end 
        
        hold off
        
        
        
        
        
        % PLOT FULL TEP
        figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]) 
        % [left bottom width height]
        
        time = TEP_bc_gav.time(251:3001) * 1000; % convert to ms
        
        for i = 1:size(TEP_bc_gav.avg,1)
            plot(time, TEP_bc_gav.avg(i,251:3001));
            hold on
        end
        
        % --- add grey rectangle ---
        
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        xlim([-20 400])
        ylim([-15 15])
        
        title('TEP EEG Signal 400-1000Hz', 'FontSize', 24, 'FontName', 'Arial');
        set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)
        
        % add grey patch spanning -0.2 to 6.4 ms across whole y-range
        if Paradigm == 'singlepulse' % singlepulse or tripulse
        
            yl = ylim;
            patch([-0.5 -0.5 2 2], [yl(1) yl(2) yl(2) yl(1)], ...
                  [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
        
        elseif  Paradigm == 'tripulse'
            yl = ylim;
            patch([-0.5 -0.5 6.4 6.4], [yl(1) yl(2) yl(2) yl(1)], ...
                  [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
        end 

elseif pulse == "with_pulse" % with_pulse  or  no_pulse



        % PLOT FULL TEP with time in ms on x-axis FOR datasets WITH TMS pulse
        t_ms = TEP_bc_gav.time * 1000;                  % seconds -> ms
        if iscolumn(t_ms), t_ms = t_ms.'; end           % make row to match avg
        
        figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]); hold on
        for i = 1:size(TEP_bc_gav.avg,1)
            plot(t_ms, TEP_bc_gav.avg(i,:), 'LineWidth', 2);
        end
        
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        ylim([-100 100])
        xlim([t_ms(1) t_ms(end)])                       % fit to full time range
        title('TEP EEG Signal 0.2–2000 Hz', 'FontSize', 24, 'FontName', 'Arial');
        set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)




end 



%% 2------------------> GMFP analyses and plots


% --- BAND-PASS ---
cfg = [];
cfg.channel    = 'all';
cfg.bpfilter   = 'yes';
cfg.bpfreq     = [2 2000];       % passband in Hz
cfg.bpfilttype = 'firws';
cfg.demean     = 'yes';
data_gmfp = ft_preprocessing(cfg, data);

cfg = [];
cfg.channel = {'all'};  % alternatively select the channels of interet. Example cf ro channels around left M1:  cfg.channel = {'C1','C3','CP1','CP3'} ; 
cfg.latency  = [-0.2 0.4];
cfg.keeptrials = 'no'; 
TEP = ft_timelockanalysis(cfg,data_gmfp);

cfg = [];
cfg.baseline = [-0.2 -.05];
cfg.parameter = 'avg';
cfg.baselinetype = 'absolute';
TEP_bc = ft_timelockbaseline(cfg, TEP); %bc = baseline corrected


cfg = [];
cfg.method = 'amplitude';
cfg.channel = {'all'};  % Inc ase you want the Local Mean field power, here select the electrodes of interest. Example channels around left M1:  cfg.channel = {'C1','C3','CP1','CP3'} ; 
TEP_bc_gmfp = ft_globalmeanfield(cfg, TEP_bc);

cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
TEP_bc_gmfp_gav       = ft_timelockgrandaverage(cfg, TEP_bc_gmfp); 


cfg = [];
cfg.channel    = 'all';
cfg.bpfilttype = 'firws';
cfg.demean     = 'yes';
cfg.baselinewindow = [-0.2 -0.05] % in seconds, the default is the complete trial (default = 'all')
cfg.bsfreq        = [48 52] %bandstop frequency range, specified as [low high]
TEP_bc_gmfp_gav = ft_preprocessing(cfg, TEP_bc_gmfp_gav);


if pulse == "no_pulse" % with_pulse  or  no_pulse

            
            % --- select the original window ---
            idx   = 951:1076;                              % samples you were plotting
            t_ms  = TEP_bc_gmfp_gav.time(idx) * 1000;      % ms
            seg   = TEP_bc_gmfp_gav.avg(:, idx);           % [channels x N] (GMFP row(s))
            
            % --- original and target sampling rates ---
            Fs_old = 1 / median(diff(TEP_bc_gmfp_gav.time));  % Hz (time is in seconds)
            Fs_new = 25000;                                   % target Hz
            [p,q]  = rat(Fs_new / Fs_old);                    % 25000/2000 = 25/2
            
            % --- upsample (polyphase resampling with anti-imaging filter) ---
            seg_up = resample(seg.', p, q).';                 % resample along time (transpose in/out)
            
            % --- new time axis (match original window endpoints) ---
            t_up_ms = linspace(t_ms(1), t_ms(end), size(seg_up, 2));
            
            % --- plot GMFP (upsampled) ---
            figure('Units','normalized','Position',[0.1 0.1 0.4 0.4])
            hold on
            for i = 1:size(seg_up,1)
                plot(t_up_ms, seg_up(i,:), 'LineWidth', 5);
            end
            
            xlabel('Time (ms)')
            ylabel('Amplitude (µV)')
            xlim([-5 15])
            ylim([-5 5])
            title('GMFP EEG Signal 0.2–2000Hz', 'FontSize', 24, 'FontName', 'Arial');
            set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)
            
            % grey patch spanning -0.5 to 5.8 ms across whole y-range
            if Paradigm == 'singlepulse' % singlepulse or tripulse
            
                yl = ylim;
                patch([-0.5 -0.5 2 2], [yl(1) yl(2) yl(2) yl(1)], ...
                      [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
            
            elseif  Paradigm == 'tripulse'
                yl = ylim;
                patch([-0.5 -0.5 6.4 6.4], [yl(1) yl(2) yl(2) yl(1)], ...
                      [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
            end 
            
            hold off
            
            
            
            %
            % Plot GMFP of FULL TEP
            % --- select the original window ---
            idx   = 251:3001;                               % samples you were plotting
            t_ms  = TEP_bc_gmfp_gav.time(idx) * 1000;       % ms
            seg   = TEP_bc_gmfp_gav.avg(:, idx);            % [rows x N] (GMFP row(s))
            
            % --- original and target sampling rates ---
            Fs_old = 1 / median(diff(TEP_bc_gmfp_gav.time));  % Hz (time is in seconds)
            Fs_new = 25000;                                    % target Hz
            [p,q]  = rat(Fs_new / Fs_old);                     % 25000/2000 = 25/2
            
            % --- upsample (polyphase resampling with anti-imaging filter) ---
            seg_up = resample(seg.', p, q).';                  % resample along time (transpose in/out)
            
            % --- new time axis (match original window endpoints) ---
            t_up_ms = linspace(t_ms(1), t_ms(end), size(seg_up, 2));
            
            % --- plot GMFP (upsampled) ---
            figure('Units','normalized','Position',[0.1 0.1 0.4 0.4])
            hold on
            for i = 1:size(seg_up,1)
                plot(t_up_ms, seg_up(i,:), 'LineWidth', 5);
            end
            
            xlabel('Time (ms)')
            ylabel('Amplitude (µV)')
            xlim([-20 400])
            ylim([-5 5])
            title('GMFP EEG Signal 0.2–2000Hz', 'FontSize', 24, 'FontName', 'Arial');
            set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)
            
            % grey patch spanning -0.5 to 5.8 ms across whole y-range
            if Paradigm == 'singlepulse' % singlepulse or tripulse
            
                yl = ylim;
                patch([-0.5 -0.5 2 2], [yl(1) yl(2) yl(2) yl(1)], ...
                      [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
            
            elseif  Paradigm == 'tripulse'
                yl = ylim;
                patch([-0.5 -0.5 6.4 6.4], [yl(1) yl(2) yl(2) yl(1)], ...
                      [0.8 0.8 0.8], 'FaceAlpha',0.7, 'EdgeColor','none');
            end 
            
            hold off

elseif pulse == "with_pulse" % with_pulse  or  no_pulse

        % PLOT FULL TEP with time in ms on x-axis FOR data WITH PULSE
        t_ms = TEP_bc_gmfp_gav.time * 1000;                  % seconds -> ms
        if iscolumn(t_ms), t_ms = t_ms.'; end           % make row to match avg
        
        figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]); hold on
        for i = 1:size(TEP_bc_gmfp_gav.avg,1)
            plot(t_ms, TEP_bc_gmfp_gav.avg(i,:), 'LineWidth', 2);
        end
        
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        ylim([-100 100])
        xlim([t_ms(1) t_ms(end)])                       % fit to full time range
        title('TEP EEG Signal 0.2–2000 Hz', 'FontSize', 24, 'FontName', 'Arial');
        set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)

end


%% -------------------> FOR TOPOGRAPHIES A TTHE SENSOR LEVEL RUN THIS


% Run this for sensor-level topographies of interest 
% Adjust time points for plots
for i=[0.0058, 0.0064, 0.0075] % timepoints to plot in in sec. For refence 1msec = 0.001s

cfg = [];
cfg.xlim = [i i];
cfg.channel = {'all'};
cfg.zlim = [-5 5];
cfg.layout = 'biosemi64.lay';
cfg.comment   =  'no';
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    cmap = (flipud(brewermap(64,'RdBu')))
    cfg.colormap = cmap;
    cfg.colormap            = cmap;

ft_topoplotER(cfg,TEP_bc_gav); 
c = colorbar
c.FontSize = 15
c.Label.String = 'Amplitude (uV)','FontSize',20;
set(gca, 'FontName', 'Times New Roman');  % Set axes text to Times New Roman
set(gcf, 'DefaultTextFontName', 'Times New Roman');  % Set general text to Times New Roman
title([num2str(i*1000), 'ms'],'FontSize',20)
set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)
end


%% 3. ---------------------------> TIME FRQEUENCY ANALYSES


% --- BAND-PASS ---
cfg = [];
cfg.channel    = 'all';
cfg.bpfilter   = 'yes';
cfg.bpfreq     = [2 1000];       % passband in Hz
%cfg.bpfilttype = 'firws';
cfg.demean     = 'yes';
data_tf = ft_preprocessing(cfg, data);



% ERSP
        cfg = []
        cfg.method= 'mtmconvol';
        cfg.output= 'pow';
        cfg.foi= [2:1000];                          %the frequency range you want to have computed
        cfg.taper = 'hanning';  
        cfg.t_ftimwin= ones(1,numel(cfg.foi))*0.5;      %this is the time window for the moving window; you can %change it between 0.4 and 1; larger number improves frequency resolution
        cfg.toi = [-0.05:0.01:0.35]; % in sec
        cfg.keeptrials = 'yes';
        data_tf = ft_freqanalysis(cfg, data_tf);

        cfg = [];
        cfg.baseline = [-0.05 -0.005]; % in sec
        cfg.baselinetype = 'absolute'; % 'absolute'
        data_tf_b = ft_freqbaseline(cfg, data_tf);
               
        %cfg = [];
        %cfg.baseline = [-6.4 8.4];
        %cfg.baselinetype = 'zscore'; % 'absolute'
        %data_tf_b = ft_freqbaseline(cfg, data_tf_b);
    
        cfg = [];
        data_tf_b_final= ft_freqdescriptives(cfg, data_tf_b);

%% Preinspect TF and topographies
cfg              = [];
cfg.channel      = 'all';
cfg.interactive     = 'yes';
cfg.showoutline     = 'yes';
cfg.layout = 'biosemi64.lay'; 
%cfg.baseline        = [-2.7 -2.5];
%cfg.baselinetype    = 'relchange';
cfg.xlim            = [-0.005  0.05];
cfg.ylim            = [100  900];
%cfg.zlim            = 'maxabs';
cfg.parameter     = 'powspctrm';  % plot the t-value 
cfg.colormap = '*RdBu';
cfg.colorbar     = 'yes';
figure
ft_singleplotTFR(cfg, data_tf_b_final);
    
%% Plot grand average across sensor TF

freq = data_tf_b_final.powspctrm(:,:,:);
tim_interp = linspace(-0.05, 0.15, 512);
freq_interp = linspace(100, 1000, 512);

% We need to make a full time/frequency grid of both the original and
% interpolated coordinates. Matlab's meshgrid() does this for us:

[tim_grid_orig, freq_grid_orig] = meshgrid(data_tf_b_final.time(:),data_tf_b_final.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

meanpow = squeeze(nanmean(freq, 1));

% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow,tim_grid_interp, freq_grid_interp, 'spline');

% plot time-frequency data
figure();
set(gcf, 'Position', [100, 100, 2000, 600]); % Manually adjust figure size [x y width height]
imagesc(tim_interp*1000, freq_interp, pow_interp); % convert x-axis to ms
axis xy;
xlim([-50 150]);  % in ms now
xlabel('Time (ms)', 'FontSize', 20,'FontName', 'Arial');

%xticks([-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8]);
%xticklabels({'-6' '-5' '-4' '-3' '-2' '-1','0', '1', '2' '3', '4' '5' '6' '7' '8'});

ylabel('Frequency (Hz)', 'FontSize', 20,'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 20,'FontName', 'Arial');
title(['ERSP Mean All Sensors']);

ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
cmap = (flipud(brewermap(64,'RdBu')))
colormap(cmap)
colorbar;
h = colorbar('LineWidth',1.5);
%caxis([-5 5])
ylabel(h, 'ERSP Absolute Amplitude (µV)','FontSize', 20, 'FontWeight','bold');
set(gca,'FontName','Arial','fontsize',20,'FontWeight','bold','LineWidth',1.5)

% add grey patch spanning -0.2 to 6.4 ms across whole y-range
yl = ylim;
patch([-0.5 -0.5 6.4 6.4], [yl(1) yl(2) yl(2) yl(1)], ...
      [0.8 0.8 0.8], 'FaceAlpha',1, 'EdgeColor','none');


%%

%END
