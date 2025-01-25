close all
clear all
clc

%%%% In this script the ERD/ERS of each subject is computed and visualized.
%%%% The resutls are than saved so that the grand_average.mat file     
%%%% can load them to compute the grand average of ERD/ERS time course 

% loading laplacian filter
load("laplacian16.mat")
load('chanlocs16.mat');

mainFolder     = fullfile(pwd);
allEntries     = dir(mainFolder);                                
subjectFolders = allEntries([allEntries.isdir] & ~startsWith({allEntries.name}, '.'));     

% channel to name correspondance
channelMap = {    1, 'Fz';    2, 'FC3';    3, 'FC1';    4, 'FCz';
    5, 'FC2';    6, 'FC4';    7, 'C3';    8, 'C1';    9, 'Cz';
    10, 'C2';    11, 'C4';    12, 'CP3';    13, 'CP1';    14, 'CPz';
    15, 'CP2';    16, 'CP4'};

for i = 1:length(subjectFolders)      % i changes with the subject

    subjectFolder = fullfile(mainFolder, subjectFolders(i).name);    
    subjectPrefix = extractBefore(subjectFolders(i).name, '_');      
    gdfFiles      = dir(fullfile(subjectFolder, '*.gdf'));                 % all the gdf file corresponding to the current subjects               
    % gdf files selection
    selectedFiles = gdfFiles(contains({gdfFiles.name}, 'offline') & startsWith({gdfFiles.name}, subjectPrefix));
    
    %% file loading and concatenation
    filePath = fullfile(subjectFolder, selectedFiles(1).name); % selection of the first file
    [s, h]   = sload(filePath);
    s(:,17)  = [];                                  % elimination of reference channel      
    signal   = s;
    pos      = h.EVENT.POS;
    typ      = h.EVENT.TYP;
    dur      = h.EVENT.DUR;
    [length_sig, channels] = size(signal);
    Fs       = h.SampleRate;

    for j = 2:length(selectedFiles) 
        filePath = fullfile(subjectFolder, selectedFiles(j).name);
        [s, h]   = sload(filePath);
        s(:,17)  = [];
        
        % updating of signal infos
        signal     = vertcat(signal, s);
        pos        = vertcat(pos, h.EVENT.POS + length_sig); 
        typ        = vertcat(typ, h.EVENT.TYP);
        dur        = vertcat(dur, h.EVENT.DUR);
        length_sig = length(signal);
    end

    %% application of laplacian filter for spacial filtering
    data = signal*lap;                              % spatial filter

    %% filtering the signal in μ and β bands
    n = 4;                                          % filter order

    % MU FILTERING (band-pass 8-12 Hz)
    mu_sig  = data;
    Wn_low  = 8/(Fs/2);                              % Hz normalized in [0 1]
    Wn_high = 12/(Fs/2);                             % Hz normalized in [0 1] 
    [b,a]   = butter(n, [Wn_low Wn_high], 'bandpass');
    mu_sig  = filtfilt(b, a, mu_sig);
            
    % BETA FILTERING (band-pass 18-22 Hz)
    beta_sig = data;
    Wn_low   = 18/(Fs/2);                           % Hz normalized in [0 1]
    Wn_high  = 22/(Fs/2);                           % Hz normalized in [0 1]
    [b,a]    = butter(n, [Wn_low Wn_high], 'bandpass');
    beta_sig = filtfilt(b, a, beta_sig);
    [H_beta, F] = freqz(b,a,512);

    %% squaring of the signal
    sq_mu   = mu_sig.^2;
    sq_beta = beta_sig.^2;

    %% applycation of moving average (1-second window)
    data_beta = movmean(sq_beta,Fs);
    data_mu   = movmean(sq_mu,Fs);

    %% 3. Trial extraction and task-based separation
    % each trial lasts from the fixation cross -786- to the end of the 
    % continuous feedback period -781-
    
    inizio    = pos(typ == 786);
    fine      = pos(typ == 781)+dur(typ == 781);
    inizio_CF = pos(typ == 781) - inizio;
    inizio_CF = max(inizio_CF);
    
    num_trials    = length(inizio);
    length_trials = min(fine - inizio);
    
    trials_beta = zeros(length_trials, channels, num_trials);
    trials_mu   = zeros(length_trials, channels, num_trials);
    
    for k = 1:num_trials
         curr_trial_beta = data_beta(inizio(k):fine(k),:);
         curr_trial_mu   = data_mu(inizio(k):fine(k),:);
         if size(curr_trial_beta,1) > size(trials_beta,1)
             curr_trial_beta(length_trials+1:end, :) = [];
             curr_trial_mu(length_trials+1:end, :)   = [];
         end
    
         trials_beta(:,:,k) = curr_trial_beta;
         trials_mu(:,:,k)   = curr_trial_mu;
    end
    
    types = typ(3:4:end);
    Pk    = zeros(num_trials,1);                       % index to identify trials by task 
    Pk(types==773) = 773;
    Pk(types==771) = 771;

    %% Event related desynchronization and synchronization (ERD/ERS)
    % The ERD/ERS time course is computed by considering the fixation 
    % period as reference period and the continuous feedback as activity period.
    
    % creation of FixData [samples x channels x trials]
    dur_fix = min(dur(2:4:end,1));
    
    FixDataBeta = trials_beta(1:dur_fix,:,:);
    FixDataMu   = trials_mu(1:dur_fix,:,:);
    
    % computation of the ERD/ERS for both band
    ReferenceBeta = repmat(mean(FixDataBeta), [length_trials 1 1]);
    ERD_Beta = 100 * (trials_beta - ReferenceBeta)./ ReferenceBeta;
    
    ReferenceMu = repmat(mean(FixDataMu), [length_trials 1 1]);
    ERD_Mu = 100 * (trials_mu - ReferenceMu)./ ReferenceMu;

    %%  Temporal visualization
    
    % channel selection
    Ch = 7;     % correspond to C3
    
    % averaging of ERD/ERS across trials (μ band)
    ERD_mu_feet = mean(ERD_Mu(:,:,Pk == 771),3);
    ERD_mu_hand = mean(ERD_Mu(:,:,Pk == 773),3);
    
    % SE computation 
    std_ERD_mu_feet = std(ERD_Mu(:,:,Pk == 771),0,3);
    se_ERD_mu_feet  = std_ERD_mu_feet/sqrt(45); 
    std_ERD_mu_hand = std(ERD_Mu(:,:,Pk == 773),0,3);
    se_ERD_mu_hand  = std_ERD_mu_feet/sqrt(45);
    
    % averaging of ERD/ERS across trials (β band)
    ERD_beta_feet = mean(ERD_Beta(:,:,Pk == 771),3);
    ERD_beta_hand = mean(ERD_Beta(:,:,Pk == 773),3);
    
    % SE computation
    std_ERD_beta_feet = std(ERD_Mu(:,:,Pk == 771),0,3);
    se_ERD_beta_feet  = std_ERD_beta_feet/sqrt(45);
    std_ERD_beta_hand = std(ERD_Mu(:,:,Pk == 773),0,3);
    se_ERD_beta_hand  = std_ERD_beta_feet/sqrt(45);
    
    % time vector
    t = 0 : 1/Fs : (length_trials - 1)/Fs;
    
    % Plot of the average and the standard error of the ERD/ERS over time for the two classes
    figure
    subplot 121
    hold on
    plot(t, ERD_mu_hand(:,Ch), 'r', 'DisplayName','Both Hands')
    plot(t, ERD_mu_feet(:,Ch), 'g', 'DisplayName','Both Feet')
    plot(t,ERD_mu_hand(:,Ch) - se_ERD_mu_hand(:,Ch), ':r', t, ERD_mu_hand(:,Ch) + se_ERD_mu_hand(:,Ch), ':r')
    plot(t,ERD_mu_feet(:,Ch) - se_ERD_mu_feet(:,Ch), ':g', t, ERD_mu_feet(:,Ch) + se_ERD_mu_feet(:,Ch), ':g')
    grid on
    title(['ERD/ERS in mu band | Mean +- SE | Channel ', channelMap{Ch,2}])
    xlabel('time [s]')
    xline(dur_fix/Fs, 'DisplayName','end reference')
    xline(inizio_CF/Fs, 'DisplayName', 'start activity')
    legend({'Both hands', 'Both feet'})
    
    subplot 122
    hold on
    plot(t, ERD_beta_hand(:,Ch), 'r', 'DisplayName','Both Hands')
    plot(t, ERD_beta_feet(:,Ch), 'g', 'DisplayName','Both Feet')
    plot(t,ERD_beta_hand(:,Ch) - se_ERD_beta_hand(:,Ch), ':r', t, ERD_beta_hand(:,Ch) + se_ERD_beta_hand(:,Ch), ':r')
    plot(t,ERD_beta_feet(:,Ch) - se_ERD_beta_feet(:,Ch), ':g', t, ERD_beta_feet(:,Ch) + se_ERD_beta_feet(:,Ch), ':g')
    grid on
    title(['ERD in beta band | Mean +- SE | Channel ', channelMap{Ch,2}])
    xlabel('time [s]')
    xline(durFix/Fs, 'DisplayName','end reference')
    xline(CF_start/Fs, 'DisplayName', 'start activity')
    legend({'Both hands', 'Both feet'})

    pause
    
    %% Spatial visualization (only for μ band)

    %average across trials for the two task and in time (μ band)
    ERD_Ref_771_Mu = mean(mean(ERD_Mu(1:dur_fix, :, Pk == 771),3),1);
    ERD_Act_771_Mu = mean(mean(ERD_Mu(inizio_CF:end, :, Pk == 771), 3),1);
    ERD_Ref_773_Mu = mean(mean(ERD_Mu(1:dur_fix, :, Pk == 773),3),1);
    ERD_Act_773_Mu = mean(mean(ERD_Mu(inizio_CF:end, :, Pk == 773), 3),1);

    % individual topoplots
    figure
    subplot 221
    topoplot(squeeze(ERD_Ref_771_Mu), chanlocs16);
    colorbar
    clim([-50, 100])
    title('ERD/ERS during fixation (mu band) - Both Hands')
    subplot 222
    topoplot(squeeze(ERD_Act_771_Mu), chanlocs16);
    colorbar
    clim([-50, 100])
    title('ERD/ERS during activity (mu band) - Both Hands')
    subplot 223
    topoplot(squeeze(ERD_Ref_773_Mu), chanlocs16);
    colorbar
    clim([-50, 100])
    title('ERD/ERS during activity (mu band) - Both Hands')
    subplot 224
    topoplot(squeeze(ERD_Act_773_Mu), chanlocs16);
    colorbar
    clim([-50, 100])
    title('ERD/ERS during activity (mu band) - Both Feet')

    pause

    % saving the results in the same folder with the correspondent gdf
    % file name but .mat format
    [~, fileName, ~] = fileparts(selectedFiles(j).name);
    matFilePath = fullfile(subjectFolder, strcat(fileName, '.mat'));
    save(matFilePath, 'ERD_Beta', 'ERD_Mu', 'dur_fix', 'inizio_CF', 'Pk');

end


