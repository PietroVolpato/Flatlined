close all
clear all
clc

main = fullfile(pwd);
all = dir(main);                                
subjFolders = all([all.isdir] & ~startsWith({all.name}, '.'));     

for i = 1:length(subjFolders)  
    psd = [];
    % loading of all mat files containing PSD
    cartellaCorrente = fullfile(main, subjFolders(i).name);    
    prefissoSubj = extractBefore(subjFolders(i).name, '_');     
    matFiles = dir(fullfile(cartellaCorrente, '*PSD.mat'));                
    selectedFiles = matFiles(startsWith({matFiles.name}, prefissoSubj));
    
    for j = 1:length(selectedFiles)
        filePath = fullfile(cartellaCorrente, selectedFiles(j).name);
        load(filePath) 
        
        if isempty(psd)                     % for the first iteration
            psd = PSD;                      % PSD: [ windows x freq x channels ]
            pos = events.POS;
            typ = events.TYP;
            dur = events.DUR;
            length_sig = size(PSD,1);
            disp('loading of the first file completed.')
        else                                % updating
            psd = vertcat(psd, PSD);
            pos = vertcat(pos, events.POS + length_sig);
            typ = vertcat(typ, events.TYP);
            dur = vertcat(dur, events.DUR);
            
            length_sig = size(PSD,1) + length_sig;
            disp(['loading of the ', num2str(j), 'th file for the ', num2str(i),'° subject completed.'])
        end
    end
    
    clear events PSD nomefile file

    EVENT.POS = pos;
    EVENT.TYP = typ;
    EVENT.DUR = dur;
    
    channels = size(psd,3);
    f = size(psd,2);
    num_windows = length_sig;
    
    clear pos typ dur length_sig

    %% trial extraction
    % Each trial ranges from the event related to the fixation cross (TYP=786) 
    % to the end of the event related to the continuous feedback (TYP=781)

    inizio = EVENT.POS(EVENT.TYP == 786);
    fine = EVENT.POS(EVENT.TYP == 781)+EVENT.DUR(EVENT.TYP == 781);
    num_trials = length(inizio);
    dur = fine - inizio;
    length_trials = min(dur);

    activity = zeros(length_trials, f, channels, num_trials);

    for k = 1:num_trials
         curr_trial = psd(inizio(k):fine(k),:,:);

         if size(curr_trial,1) > size(activity,1)
             curr_trial(length_trials+1:end, :, :) = [];
         end
         activity(:,:,:,k) = curr_trial;
    end

    types = EVENT.TYP(3:4:end);
    Pk = zeros(num_trials,1);                       % index to distinguish trials by task
    Pk(types==773) = 773;
    Pk(types==771) = 771;

    %% ERD/ERS computation
    % creation of reference [windows x freq x channels x trials]
    fix = EVENT.DUR(2:4:end,1);
    reference = activity(1:min(fix), :, :, :);

    baseline = repmat(mean(reference), [size(activity, 1) 1 1 1]);
    ERD = log(activity./ baseline);                 % ERD: [ windows x freq x channel x trials]
    % or 
    % ERD = 100 * (activity - baseline)./ baseline;
    % ERD: [ windows x freq x channel x trials]

    %% visualization of PSD

    wshift = 0.0625; 
    t = 0 : wshift : wshift*length_trials;
    
    % selection of  meaningful frequencies (from 4 Hz to 48 Hz, step 2 Hz)
    % and channels
    freq = (4:2:48);     
    m_ch = [7, 9, 11];
    sel_ERD = ERD(:, :, m_ch, :);
    
    % averaging across trials
    sel_ERD_hand = mean(sel_ERD(:,:,:,Pk == 773),4);
    sel_ERD_feet = mean(sel_ERD(:,:,:,Pk == 771),4);
    
    limits = [-1,0.5]; % color bar limits
    
    figure(i)
    subplot 231
    imagesc(t, freq, sel_ERD_hand(:,:,1)')
    title('BOTH HANDS, channel C3')
    colorbar
    xlabel('time (s)')
    ylabel('frequency [Hz]')
    clim(limits)
    xline(3)
    xline(4)
    axis xy
    
    subplot 232
    imagesc(t, freq, sel_ERD_hand(:,:,2)')
    title('BOTH HANDS, channel Cz')
    colorbar
    xlabel('time (s)')
    ylabel('frequency [Hz]')
    clim(limits)
    xline(3)
    xline(4)
    axis xy
    
    subplot 233
    imagesc(t, freq, sel_ERD_hand(:,:,3)')
    title('BOTH HANDS, channel C4')
    colorbar
    xlabel('time (s)')
    ylabel('frequency [Hz]')
    clim(limits)
    xline(3)
    xline(4)
    axis xy
    
    subplot 234
    imagesc(t, freq, sel_ERD_feet(:,:,1)')
    title('BOTH FEET, channel C3')
    colorbar
    xlabel('time (s)')
    ylabel('frequency [Hz]')
    clim(limits)
    xline(3)
    xline(4)
    axis xy
    
    subplot 235
    imagesc(t, freq, sel_ERD_feet(:,:,2)')
    title('BOTH FEET, channel Cz')
    colorbar
    xlabel('time (s)')
    ylabel('frequency [Hz]')
    clim(limits)
    xline(3)
    xline(4)
    axis xy
    
    subplot 236
    imagesc(t, freq, sel_ERD_feet(:,:,3)')
    title('BOTH FEET, channel C4')
    colorbar
    xlabel('time (s)')
    ylabel('frequency [Hz]')
    clim(limits)
    xline(3)
    xline(4)
    axis xy
    
    colormap hot
    pause

    disp(['Analisy of the ', num2str(i),'° subject completed.'])
end

fprintf('Analisy of all subjects completed.\n');



