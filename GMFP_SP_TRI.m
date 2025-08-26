%% COMPARE GMFP OF TWO DIFFERENT SELECTED CONDITIONS



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

% CONDITION 1
% ---------------------------------------->  Modify subject condition to load
Subject = 'XEB' % XPB or XEB or XTC
Coil = 'coilB65' % coilB35 or coilB65
Target = 'SFG' % SFG or SPL
Orientation = 'LM' % LM or APPA
Intensity = '100' % 70 or 80 or 90 or 100 
MT_MSO = 'RMT'% RMT or MSO
ParadigmSP = 'tripulse' % singlepulse or tripulse
pulseSP = 'no_pulse'; % with_pulse  or  no_pulse
processing_pipeline = 'basic'  % basic   or advanced



% LOAD CONDITION1  
nameSP = [Subject,' ',Coil,' ',Target,' ',Orientation,' ',Intensity,MT_MSO,' ',ParadigmSP];
name_datasetSP = [Subject,'_',Coil,'_',Target,'_',Orientation,'_',Intensity,MT_MSO,'_',ParadigmSP,'_',pulseSP,'_cleaned_pipeline_',processing_pipeline,'.set']; 
path_dataset = ['/mnt/projects/PARADISE/PARADISE_1/',Subject,'/EEG_clean/'];
% LOAD file with fieltrip
cfg = []; 
     cfg.lpfilter        = 'no';
     cfg.dataset = [path_dataset, name_datasetSP];
     dataSP = ft_preprocessing(cfg);
% <--------------------------------------------------------------------
%% LOAD CONDITION 2

% CONDITION 1
% ---------------------------------------->  Modify subject condition to load
Subject = 'XEB' % XPB or XEB or XTC
Coil = 'coilB35' % coilB35 or coilB65
Target = 'SFG' % SFG or SPL
Orientation = 'LM' % LM or APPA
Intensity = '100' % 70 or 80 or 90 or 100 
MT_MSO = 'MSO'% RMT or MSO
ParadigmTRI = 'tripulse' % singlepulse or tripulse
pulseTRI = 'no_pulse'; % with_pulse  or  no_pulse
processing_pipeline = 'basic'  % basic   or advanced


% Define TRI file to load 
nameTRI = [Subject,' ',Coil,' ',Target,' ',Orientation,' ',Intensity,MT_MSO,' ',ParadigmTRI];
name_datasetTRI = [Subject,'_',Coil,'_',Target,'_',Orientation,'_',Intensity,MT_MSO,'_',ParadigmTRI,'_',pulseTRI,'_cleaned_pipeline_',processing_pipeline,'.set']; 
path_dataset = ['/mnt/projects/PARADISE/PARADISE_1/',Subject,'/EEG_clean/'];
% LOAD file with fieltrip
cfg = []; 
     cfg.lpfilter        = 'no';
     cfg.dataset = [path_dataset, name_datasetTRI];
     dataTRI = ft_preprocessing(cfg);

%% 1. DO GMFP



% --- SP  ---
cfg = [];
cfg.channel    = 'all';
cfg.bpfilter   = 'yes';
cfg.bpfreq     = [2 2000];       % passband in Hz
cfg.bpfilttype = 'firws';
cfg.demean     = 'yes';
data_gmfpSP = ft_preprocessing(cfg, dataSP);

cfg = [];
cfg.channel = {'all'};  % alternatively select the channels of interet. Example cf ro channels around left M1:  cfg.channel = {'C1','C3','CP1','CP3'} ; 
cfg.latency  = [-0.2 0.4];
cfg.keeptrials = 'no'; 
TEPSP = ft_timelockanalysis(cfg,data_gmfpSP);

cfg = [];
cfg.baseline = [-0.2 -.05];
cfg.parameter = 'avg';
cfg.baselinetype = 'absolute';
TEP_bcSP = ft_timelockbaseline(cfg, TEPSP); %bc = baseline corrected


cfg = [];
cfg.method = 'amplitude';
cfg.channel = {'all'};  % Inc ase you want the Local Mean field power, here select the electrodes of interest. Example channels around left M1:  cfg.channel = {'C1','C3','CP1','CP3'} ; 
TEP_bc_gmfpSP = ft_globalmeanfield(cfg, TEP_bcSP);

cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
TEP_bc_gmfp_gavSP       = ft_timelockgrandaverage(cfg, TEP_bc_gmfpSP); 


cfg = [];
cfg.channel    = 'all';
cfg.bpfilttype = 'firws';
cfg.demean     = 'yes';
cfg.baselinewindow = [-0.2 -0.05] % in seconds, the default is the complete trial (default = 'all')
cfg.bsfreq        = [48 52] %bandstop frequency range, specified as [low high]
TEP_bc_gmfp_gavSP = ft_preprocessing(cfg, TEP_bc_gmfp_gavSP);




% --- TRI  ---
cfg = [];
cfg.channel    = 'all';
cfg.bpfilter   = 'yes';
cfg.bpfreq     = [2 2000];       % passband in Hz
cfg.bpfilttype = 'firws';
cfg.demean     = 'yes';
data_gmfpTRI = ft_preprocessing(cfg, dataTRI);

cfg = [];
cfg.channel = {'all'};  % alternatively select the channels of interet. Example cf ro channels around left M1:  cfg.channel = {'C1','C3','CP1','CP3'} ; 
cfg.latency  = [-0.2 0.4];
cfg.keeptrials = 'no'; 
TEPTRI = ft_timelockanalysis(cfg,data_gmfpTRI);

cfg = [];
cfg.baseline = [-0.2 -.05];
cfg.parameter = 'avg';
cfg.baselinetype = 'absolute';
TEP_bcTRI = ft_timelockbaseline(cfg, TEPTRI); %bc = baseline corrected


cfg = [];
cfg.method = 'amplitude';
cfg.channel = {'all'};  % Inc ase you want the Local Mean field power, here select the electrodes of interest. Example channels around left M1:  cfg.channel = {'C1','C3','CP1','CP3'} ; 
TEP_bc_gmfpTRI = ft_globalmeanfield(cfg, TEP_bcTRI);

cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
TEP_bc_gmfp_gavTRI       = ft_timelockgrandaverage(cfg, TEP_bc_gmfpTRI); 


cfg = [];
cfg.channel    = 'all';
cfg.bpfilttype = 'firws';
cfg.demean     = 'yes';
cfg.baselinewindow = [-0.2 -0.05] % in seconds, the default is the complete trial (default = 'all')
cfg.bsfreq        = [48 52] %bandstop frequency range, specified as [low high]
TEP_bc_gmfp_gavTRI = ft_preprocessing(cfg, TEP_bc_gmfp_gavTRI);



%% 3. Shift SP data to match the time of the last pulse in the TRI condition

pulseSP = strtrim(string(pulseSP));  % normalize to a string scalar
pulseTRI = strtrim(string(pulseTRI));  % normalize to a string scalar

if pulseSP == "no_pulse" && pulseTRI == "no_pulse"
        
            % PLOT GMFP
            % --- select the original window ---
            idx   = 1:3001;
            t_ms  = TEP_bc_gmfp_gavSP.time(idx) * 1000;     % ms
            seg   = TEP_bc_gmfp_gavSP.avg(:, idx);          % [channels x N]
            
            % --- resample to 25 kHz (as you already do) ---
            Fs_old = 1 / median(diff(TEP_bc_gmfp_gavSP.time));   % Hz
            Fs_new = 25000;
            [p,q]  = rat(Fs_new / Fs_old);                       % 25/2
            seg_up = resample(seg.', p, q).';                    % [channels x N_up]
            t_up_ms = linspace(t_ms(1), t_ms(end), size(seg_up,2));
            
            % --- delay by 4.4 ms with left zero-padding (keep same time axis length) ---
            ParadigmSP = strtrim(string(ParadigmSP));  % normalize to a string scalar
            if ParadigmSP == "singlepulse"
                delay_ms = 4.4;
            elseif ParadigmSP == "tripulse"
                delay_ms = 0;
            end
            dt_ms    = median(diff(t_up_ms));
            n_shift  = max(0, round(delay_ms / dt_ms));          % samples to shift
            
            if n_shift >= size(seg_up,2)
                warning('Shift exceeds signal length; result will be all zeros.');
                seg_shift = zeros(size(seg_up), 'like', seg_up);
            else
                seg_shift = [zeros(size(seg_up,1), n_shift, 'like', seg_up), ...
                             seg_up(:, 1:end-n_shift)];
            end
            % Note: t_up_ms unchanged -> the sample that was at 0 ms now appears at ~4.4 ms
            
            % --- plot ---
            figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]); hold on
            for i = 1:size(seg_shift,1)
                plot(t_up_ms, seg_shift(i,:), 'LineWidth', 5);
            end
            xlabel('Time (ms)'); ylabel('Amplitude (µV)')
            xlim([-5 15]); 
            ylim([-10 10])
            title('GMFP EEG Signal 0.2–2000 Hz', 'FontSize', 24, 'FontName', 'Arial');
            set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)
            
            
            
            
            
            idx   = 1:3001;
            t_ms  = TEP_bc_gmfp_gavTRI.time(idx) * 1000;     % ms
            seg   = TEP_bc_gmfp_gavTRI.avg(:, idx);          % [channels x N]
            
            % --- resample to 25 kHz (as you already do) ---
            Fs_old = 1 / median(diff(TEP_bc_gmfp_gavTRI.time));   % Hz
            Fs_new = 25000;
            [p,q]  = rat(Fs_new / Fs_old);                       % 25/2
            seg_up = resample(seg.', p, q).';                    % [channels x N_up]
            t_up_ms = linspace(t_ms(1), t_ms(end), size(seg_up,2));
            
            % --- delay by 4.4 ms with left zero-padding (keep same time axis length) ---
            
            ParadigmTRI = strtrim(string(ParadigmTRI));  % normalize to a string scalar
            if ParadigmTRI == "singlepulse"
                delay_ms = 4.4;
            elseif ParadigmTRI == "tripulse"
                delay_ms = 0;
            end

            dt_ms    = median(diff(t_up_ms));
            n_shift  = max(0, round(delay_ms / dt_ms));          % samples to shift
            
            if n_shift >= size(seg_up,2)
                warning('Shift exceeds signal length; result will be all zeros.');
                seg_shift = zeros(size(seg_up), 'like', seg_up);
            else
                seg_shift = [zeros(size(seg_up,1), n_shift, 'like', seg_up), ...
                             seg_up(:, 1:end-n_shift)];
            end
            % Note: t_up_ms unchanged -> the sample that was at 0 ms now appears at ~4.4 ms
            
            % --- plot ---
            %figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]); hold on
            for i = 1:size(seg_shift,1)
                plot(t_up_ms, seg_shift(i,:), 'LineWidth', 5);
            end
            xlabel('Time (ms)'); ylabel('Amplitude (µV)')
            xlim([-5 15]); 
            
            name = [Subject,' ',Coil,' ',Target,' ',Orientation,' ',Intensity,MT_MSO,' ',];
            title('GMFP EEG Signal 0.2–2000 Hz', 'FontSize', 24, 'FontName', 'Arial');
            set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)

  
            yl = ylim;
            patch([-0.5 -0.5 6.4 6.4], [yl(1) yl(2) yl(2) yl(1)], ...
                  [0.8 0.8 0.8], 'FaceAlpha',0.9, 'EdgeColor','none');

            legend(nameSP,nameTRI, 'Location','northeast');

elseif pulseSP == "with_pulse" || pulseTRI == "with_pulse"

% --- select the original window ---
            idx   = 1:3001;
            t_ms  = TEP_bc_gmfp_gavSP.time(idx) * 1000;     % ms
            seg   = TEP_bc_gmfp_gavSP.avg(:, idx);          % [channels x N]
            
            % --- resample to 25 kHz (as you already do) ---
            Fs_old = 1 / median(diff(TEP_bc_gmfp_gavSP.time));   % Hz
            Fs_new = 25000;
            [p,q]  = rat(Fs_new / Fs_old);                       % 25/2
            seg_up = resample(seg.', p, q).';                    % [channels x N_up]
            t_up_ms = linspace(t_ms(1), t_ms(end), size(seg_up,2));
            
            % --- delay by 4.4 ms with left zero-padding (keep same time axis length) ---
            ParadigmSP= strtrim(string(ParadigmSP));  % normalize to a string scalar
            if ParadigmSP == "singlepulse"
                delay_ms = 4.4;
            elseif ParadigmSP == "tripulse"
                delay_ms = 0;
            end

            dt_ms    = median(diff(t_up_ms));
            n_shift  = max(0, round(delay_ms / dt_ms));          % samples to shift
            
            if n_shift >= size(seg_up,2)
                warning('Shift exceeds signal length; result will be all zeros.');
                seg_shift = zeros(size(seg_up), 'like', seg_up);
            else
                seg_shift = [zeros(size(seg_up,1), n_shift, 'like', seg_up), ...
                             seg_up(:, 1:end-n_shift)];
            end
            % Note: t_up_ms unchanged -> the sample that was at 0 ms now appears at ~4.4 ms
            
            % --- plot ---
            figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]); hold on
            for i = 1:size(seg_shift,1)
                plot(t_up_ms, seg_shift(i,:), 'LineWidth', 5);
            end
            xlabel('Time (ms)'); ylabel('Amplitude (µV)')
            %xlim([-5 15]); 
            ylim([-10 10])
            title('GMFP EEG Signal 0.2–2000 Hz', 'FontSize', 24, 'FontName', 'Arial');
            set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)

            
            

            idx   = 1:3001;
            t_ms  = TEP_bc_gmfp_gavTRI.time(idx) * 1000;     % ms
            seg   = TEP_bc_gmfp_gavTRI.avg(:, idx);          % [channels x N]
            
            % --- resample to 25 kHz (as you already do) ---
            Fs_old = 1 / median(diff(TEP_bc_gmfp_gavTRI.time));   % Hz
            Fs_new = 25000;
            [p,q]  = rat(Fs_new / Fs_old);                       % 25/2
            seg_up = resample(seg.', p, q).';                    % [channels x N_up]
            t_up_ms = linspace(t_ms(1), t_ms(end), size(seg_up,2));
            
            % --- delay by 4.4 ms with left zero-padding (keep same time axis length) ---
            ParadigmTRI = strtrim(string(ParadigmTRI));  % normalize to a string scalar
            if ParadigmTRI == "singlepulse"
                delay_ms = 4.4;
            elseif ParadigmTRI == "tripulse"
                delay_ms = 0;
            end

            dt_ms    = median(diff(t_up_ms));
            n_shift  = max(0, round(delay_ms / dt_ms));          % samples to shift
            
            if n_shift >= size(seg_up,2)
                warning('Shift exceeds signal length; result will be all zeros.');
                seg_shift = zeros(size(seg_up), 'like', seg_up);
            else
                seg_shift = [zeros(size(seg_up,1), n_shift, 'like', seg_up), ...
                             seg_up(:, 1:end-n_shift)];
            end
            % Note: t_up_ms unchanged -> the sample that was at 0 ms now appears at ~4.4 ms
            
            % --- plot ---
            %figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]); hold on
            for i = 1:size(seg_shift,1)
                plot(t_up_ms, seg_shift(i,:), 'LineWidth', 5);
            end
            xlabel('Time (ms)'); ylabel('Amplitude (µV)')
            %xlim([-5 15]); 
            
            name = [Subject,' ',Coil,' ',Target,' ',Orientation,' ',Intensity,MT_MSO,' ',];
            title(['GMFP ' name], 'FontSize', 24, 'FontName', 'Arial');
            set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',1.5)

  
            yl = ylim;
            patch([-0.5 -0.5 6.4 6.4], [yl(1) yl(2) yl(2) yl(1)], ...
                  [0.8 0.8 0.8], 'FaceAlpha',0.9, 'EdgeColor','none');

            legend(nameSP,nameTRI, 'Location','northeast');


end 


%% END