close all
clear all
clc

main = fullfile(pwd);
all = dir(main);                               
subjFolders = all([all.isdir] & ~startsWith({all.name}, '.'));      

ch_feat = [];
freq_feat = [];

for i = 1:length(subjFolders)                            
    cartellaCorrente = fullfile(main, subjFolders(i).name);    
    prefissoSubj = extractBefore(subjFolders(i).name, '_');      
    matFiles = dir(fullfile(cartellaCorrente, '*PSD.mat'));               
    selectedFiles = matFiles(startsWith({matFiles.name}, prefissoSubj));

    for j = 1:length(selectedFiles)
        filePath = fullfile(cartellaCorrente, selectedFiles(j).name);
        load(filePath)
        h.EVENT =  events;

        
        [ T,F,A,CF,X ] = label_vector( PSD, h );  
        trials=max(T); 

        % computation of the cue vector
        for i=1:trials
            ck=unique(A(T==i));
            Ck(i)=ck(ck>0);
        end
    
        % masks for the two tasks
        mask773=(A==773) | (CF>0 & sum(T==find(Ck==773),2));
        mask771=(A==771) | (CF>0 & sum(T==find(Ck==771),2));

        Sfeat=reshape(PSD,[size(PSD,1),size(PSD,2)*size(PSD,3)]);
        P1=Sfeat(mask773,:);
        P2=Sfeat(mask771,:);
        
        % Fisher Score computation
        FS1d=abs((mean(P1,1)-mean(P2,1))./sqrt(std(P1,0,1).^2+(std(P2,0,1).^2)));
        FS = reshape(FS1d,[size(PSD,2),size(PSD,3)]);

        % Visualization
        freq_axis = (4:2:48); 

        figure
        imagesc(freq_axis,1:16,FS')
        xlabel('Band [Hz]')
        ylabel('Channel')
        title('Feature analysis for the subject:', prefissoSubj);
        colorbar
        % pause
        
        % Feature extraction
        treshold = 0.8*(max(max(FS)));
        FS_tresh = zeros(size(FS));
        FS_tresh(FS>treshold) = FS(FS>treshold);

        figure
        imagesc(freq_axis, 1:16, FS_tresh')
        xlabel('Band [Hz]')
        ylabel('Channel')
        title('Selected features for the subject:', prefissoSubj);
        colorbar
        % pause

        [ch, freq] = find(FS_tresh' ~= 0);
        ch_feat = vertcat(ch_feat, ch);
        freq_feat = vertcat(freq_feat, freq);

    end
    
end

freq_feat = (freq_feat.*2)+2;    % conversion to columns index to actual frequency

%% histograms computation
C3_freq = freq_feat(ch_feat == 7);
Cz_freq = freq_feat(ch_feat == 9);
C4_freq = freq_feat(ch_feat == 11);

figure
subplot 131
histogram(C3_freq)
xlabel('Frequency component [Hz]')
ylabel('Occurence')
title('C3')
ylim([0 10])
% xlim([10 40])

subplot 132
histogram(Cz_freq)
xlabel('Frequency component [Hz]')
ylabel('Occurence')
title('Cz')
ylim([0 10])
% xlim([10 40])

subplot 133
histogram(C4_freq)
xlabel('Frequency component [Hz]')
ylabel('Occurence')
title('C4')
ylim([0 10])
% xlim([10 40])


figure
hold on
histogram(C3_freq,"BinEdges",freq_axis(4:end-4))
histogram(Cz_freq,"BinEdges",freq_axis(4:end-4))
histogram(C4_freq,"BinEdges",freq_axis(4:end-4))
hold off
legend('C3','Cz','C4')
xlabel('Frequency component [Hz]')
ylabel('Occurence')

