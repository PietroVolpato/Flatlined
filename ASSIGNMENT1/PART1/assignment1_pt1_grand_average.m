close all
clear all
clc

%%%% In this script the grand average of ERD/ERS is computed 
%%%% and visualized                                          

main = fullfile(pwd);
all = dir(main);                                
subjFolders = all([all.isdir] & ~startsWith({all.name}, '.'));     

% useful infos and loading
load('chanlocs16.mat');
Fs = 512;
channelMap = {    1, 'Fz';    2, 'FC3';    3, 'FC1';    4, 'FCz';
    5, 'FC2';    6, 'FC4';    7, 'C3';    8, 'C1';    9, 'Cz';
    10, 'C2';    11, 'C4';    12, 'CP3';    13, 'CP1';    14, 'CPz';
    15, 'CP2';    16, 'CP4'};

Pk_tot = cell(length(subjFolders),1);

for i = 1:length(subjFolders)  % i changes with the subject

    cartellaCorrente = fullfile(main, subjFolders(i).name);    
    prefissoSubj     = extractBefore(subjFolders(i).name, '_');      
    matFiles         = dir(fullfile(cartellaCorrente, '*GA.mat'));           % all the mat file corresponding to the current subjects              
    % mat files selction
    selectedFiles    = matFiles(startsWith({matFiles.name}, prefissoSubj));
    
    % loading
    filePath = fullfile(cartellaCorrente, selectedFiles(1).name);
    load(filePath)

    CF_start(i) = inizio_CF;
    durFix(i)   = dur_fix;
    Pk_tot{i}   = Pk;    
end

clear inizio_CF dur_fix Pk 

CF_start = max(CF_start);
durFix   = min(durFix);
length_trials = size(ERD_Beta,1);

for i = 1:length(subjFolders)
    Pk = Pk_tot{i};

    % computation of ERD/ERS subject level 
    % dimension: [subj x samples x channels]
    ERD_mu_feet(i,:,:) = mean(ERD_Mu(:,:,Pk == 771),3);
    ERD_mu_hand(i,:,:) = mean(ERD_Mu(:,:,Pk == 773),3);

    ERD_beta_feet(i,:,:) = mean(ERD_Beta(:,:,Pk == 771),3);
    ERD_beta_hand(i,:,:) = mean(ERD_Beta(:,:,Pk == 773),3);
end

% Grand average computation
GA_ERD_mu_feet = squeeze(mean(ERD_mu_feet,1));
GA_ERD_mu_hand = squeeze(mean(ERD_mu_hand,1));
GA_ERD_beta_feet = squeeze(mean(ERD_beta_feet,1));
GA_ERD_beta_hand = squeeze(mean(ERD_beta_hand,1));

%% Temporal Visualization
Ch = 7;                                                 % corresponds to C3
% 7 = 'C3' or 9 = 'Cz' or 11 = 'C4'
% time course
t = 0 : 1/Fs : (length_trials - 1)/Fs;

figure
subplot 121
grid on
hold on
plot(t, GA_ERD_mu_hand(:,Ch), 'k', LineWidth=2)
plot(t, ERD_mu_hand(:,:,Ch), '--')
title(['ERD/ERS in μ band| Hands MI | Channel ', channelMap{Ch,2}])
xlabel('time [s]')
legend({'Grand Average'})
xline(durFix/Fs, 'k','DisplayName','end reference', LineWidth=1.2)
xline(CF_start/Fs, 'k','DisplayName', 'start activity', LineWidth=1.2)

subplot 122
grid on
hold on
plot(t, GA_ERD_mu_feet(:,Ch), 'k', LineWidth=2)
plot(t, ERD_mu_feet(:,:,Ch), '--')
title(['ERD/ERS in μ band | Feet MI | Channel ', channelMap{Ch,2}])
xlabel('time [s]')
legend({'Grand Average'})
xline(durFix/Fs, 'k','DisplayName','end reference', LineWidth=1.2)
xline(CF_start/Fs, 'k','DisplayName', 'start activity', LineWidth=1.2)

figure
subplot 121
grid on
hold on
plot(t, GA_ERD_beta_hand(:,Ch), 'k', LineWidth=2)
plot(t, ERD_beta_hand(:,:,Ch), '--')
title(['ERD/ERS in β band | Hands MI | Channel ', channelMap{Ch,2}])
xlabel('time [s]')
legend({'Grand Average'})
xline(durFix/Fs, 'k','DisplayName','end reference', LineWidth=1.2)
xline(CF_start/Fs, 'k','DisplayName', 'start activity', LineWidth=1.2)

subplot 122
grid on
hold on
plot(t, GA_ERD_beta_feet(:,Ch), 'k', LineWidth=2)
plot(t, ERD_beta_feet(:,:,Ch), '--')
title(['ERD/ERS in β band | Feet MI | Channel ', channelMap{Ch,2}])
xlabel('time [s]')
legend({'Grand Average'})
xline(durFix/Fs, 'k','DisplayName','end reference', LineWidth=1.2)
xline(CF_start/Fs, 'k','DisplayName', 'start activity', LineWidth=1.2)


%% Spatial visualization

% average in the time domain (μ band)
ERD_Ref_771_Mu = mean(GA_ERD_mu_hand(1:durFix, :),1);
ERD_Act_771_Mu = mean(GA_ERD_mu_hand(CF_start:end,:),1);
ERD_Ref_773_Mu = mean(GA_ERD_mu_feet(1:durFix, :),1);
ERD_Act_773_Mu = mean(GA_ERD_mu_feet(CF_start:end,:),1);

topoLim = [-50 50];

figure
subplot 221
topoplot(squeeze(ERD_Ref_771_Mu), chanlocs16);
colorbar
clim(topoLim)
title('Fixation | μ band - Both Feet')
subplot 222
topoplot(squeeze(ERD_Act_771_Mu), chanlocs16);
colorbar
clim(topoLim)
title('Activity | μ band - Both Feet')
subplot 223
topoplot(squeeze(ERD_Ref_773_Mu), chanlocs16);
colorbar
clim(topoLim)
title('Fixation | μ band - Both Hands')
subplot 224
topoplot(squeeze(ERD_Act_773_Mu), chanlocs16);
colorbar
clim(topoLim)
title('Activity | μ band - Both Hands')

% average in the time domain (β band)
ERD_Ref_771_beta = mean(GA_ERD_beta_hand(1:durFix, :),1);
ERD_Act_771_beta = mean(GA_ERD_beta_hand(CF_start:end,:),1);
ERD_Ref_773_beta = mean(GA_ERD_beta_feet(1:durFix, :),1);
ERD_Act_773_beta = mean(GA_ERD_beta_feet(CF_start:end,:),1);

topoLim = [-20 20];
figure
subplot 221
topoplot(squeeze(ERD_Ref_771_beta), chanlocs16);
colorbar
clim(topoLim)
title('Fixation | β band - Both Feet')
subplot 222
topoplot(squeeze(ERD_Act_771_beta), chanlocs16);
colorbar
clim(topoLim)
title('Activity | β band - Both Feet')
subplot 223
topoplot(squeeze(ERD_Ref_773_beta), chanlocs16);
colorbar
clim(topoLim)
title('Fixation | β band - Both Hands')
subplot 224
topoplot(squeeze(ERD_Act_773_beta), chanlocs16);
colorbar
clim(topoLim)
title('Activity | β band - Both Hands')
