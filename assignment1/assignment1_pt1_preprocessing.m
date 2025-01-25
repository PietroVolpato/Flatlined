close all
clear all
clc

% loading laplacian filter
load("laplacian16.mat")

mainFolder = fullfile(pwd);
allEntries = dir(mainFolder);                                % elenco di tutte i file nella cartella corrente
subjectFolders = allEntries([allEntries.isdir] & ~startsWith({allEntries.name}, '.'));      % estrazione delle cartelle relative ai soggetti

for i = 1:length(subjectFolders)                             % analisi per soggetto

    subjectFolder = fullfile(mainFolder, subjectFolders(i).name);    % nome della cartella del soggetto
    subjectPrefix = extractBefore(subjectFolders(i).name, '_');      % prefisso che identifica il soggetto

    gdfFiles = dir(fullfile(subjectFolder, '*.gdf'));                % tutti i gdf relativi al soggetto corrente

    % da modificare se questa analisi preliminare va fatta su tutti i file

    selectedFiles = gdfFiles(contains({gdfFiles.name}, 'offline') & startsWith({gdfFiles.name}, subjectPrefix));
    % selectedFiles = gdfFiles(startsWith({gdfFiles.name}, subjectPrefix));
    
    for j = 1:length(selectedFiles)                             % analisi a livello del singolo file
        % file loading
        filePath = fullfile(subjectFolder, selectedFiles(j).name);
        [s, h] = sload(filePath);

        % Elaborazione dei dati
        data = s(:,1:end-1)*lap;                            % spatial filter

        % PSD parameters
        wlength = 0.5;          % seconds. Length of the external window
        pshift = 0.25;          % seconds. Shift of the internal windows
        wshift = 0.0625;        % seconds. Shift of the external window
        samplerate = h.SampleRate;
        mlength = 1;            % seconds

        [PSD, f] = proc_spectrogram(data, wlength, wshift, pshift, samplerate, mlength);  % PSD computation
        % PSD : [# of windows, # of freq, # of channels]
        % f = computed frequencies

        % select  meaningful frequencies (e.g., from 4 Hz to 48 Hz, step 2 Hz)
        sel_f = f(3:25);
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
        matFilePath = fullfile(subjectFolder, strcat(fileName, '.mat'));

        save(matFilePath, 'PSD', 'events');

        fprintf('File %s elaborato con successo.\n', selectedFiles(j).name);
    end

    fprintf('Analisi per il soggetto %s completata.\n', subjectFolders(i).name);
end

fprintf('Elaborazione di tutti i soggetti completata.\n');

