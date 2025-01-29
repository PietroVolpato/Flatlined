close all
clear all
clc

%%%% In this script the PSD is computed and saved in the same subject-
%%%% folder of the processed gdf file.

% loading laplacian filter
load("laplacian16.mat")

mainFolder = fullfile(pwd);
allEntries = dir(mainFolder);                                
subjectFolders = allEntries([allEntries.isdir] & ~startsWith({allEntries.name}, '.'));      

for i = 1:length(subjectFolders)                            

    subjectFolder = fullfile(mainFolder, subjectFolders(i).name);    
    subjectPrefix = extractBefore(subjectFolders(i).name, '_');      
    gdfFiles = dir(fullfile(subjectFolder, '*.gdf'));                

    selectedFiles = gdfFiles(contains({gdfFiles.name}, 'offline') & startsWith({gdfFiles.name}, subjectPrefix));
    
    for j = 1:length(selectedFiles)                             
        % file loading
        filePath = fullfile(subjectFolder, selectedFiles(j).name);
        [s, h] = sload(filePath);

        % Elaborazione dei dati
        data = s(:,1:end-1)*lap;                            % spatial filter

        % PSD parameters
        wlength = 0.5;                          % seconds. Length of the external window
        pshift = 0.25;                          % seconds. Shift of the internal windows
        wshift = 0.0625;                        % seconds. Shift of the external window
        samplerate = h.SampleRate;
        mlength = 1;                            % seconds

        [PSD, f] = proc_spectrogram(data, wlength, wshift, pshift, samplerate, mlength);  % PSD computation
        % PSD : [# of windows, # of freq, # of channels]
        % f = computed frequencies

        % select  meaningful frequencies (e.g., from 4 Hz to 48 Hz, step 2 Hz)
        sel_f = (3:25);
        PSD = PSD(:, sel_f, :);

        % Recompute the h.EVENT.POS and .DUR with respect to the PSD windows 
        winconv = 'backward';
        POS = proc_pos2win(h.EVENT.POS, wshift*samplerate, winconv, wlength*samplerate);
        DUR = floor((h.EVENT.DUR / samplerate)/wshift);

        events.POS = POS;
        events.TYP = h.EVENT.TYP;
        events.DUR = DUR;
        
        % saving the results in the same folder with the correspondent gdf
        % file name but .mat format

        [~, fileName, ~] = fileparts(selectedFiles(j).name);
        matFilePath = fullfile(subjectFolder, strcat(fileName, 'PSD.mat'));

        save(matFilePath, 'PSD', 'events');
        fprintf('File %s elaborato con successo.\n', selectedFiles(j).name);
    end
    
    fprintf('Analisi per il soggetto %s completata.\n', subjectFolders(i).name);
end

fprintf('Elaborazione di tutti i soggetti completata.\n');

