%% THIS SCRIPT LOADS ALL THE .mat OFFLINE FILES COMPUTES THE FEATUREMAP THAT
%% ENABLES THE FEATUTRES SELECTION

clc
close all
clear all

%%
addpath("utils\")

%%
datadir='results\PSD';

SUBJ=dir(datadir);
SUBJ(1:2)=[];

for subject=1:length(SUBJ)

    subjdir=SUBJ(subject).name;
    
    D=dir(fullfile(datadir,subjdir));
    D(1:2)=[];
    
    h=struct('EVENT',[]);
    
    for file=1:length(D)
    
        if D(file).name(22)=='f'   %only offline files  
            filepath=fullfile(pwd,datadir,subjdir,D(file).name);
    
            load( filepath ); 
            [ T,F,A,CF,X ] = label_vector( PSD, h );  

            win=1:size(PSD,1);
    
            trials=max(T); 
            %cration of the vector that tells for every trial wich cue is been
            %given
            for i=1:trials
                ck=unique(A(T==i));
                Ck(i)=ck(ck>0);
            end
        
        
            %% build the masks for the two tasks
        
            mask773=(A==773) | (CF>0 & sum(T==find(Ck==773),2));
            mask771=(A==771) | (CF>0 & sum(T==find(Ck==771),2));
            
            %%   
            Sfeat=reshape(PSD,[size(PSD,1),size(PSD,2)*size(PSD,3)]);
            P1=Sfeat(mask773,:);
            P2=Sfeat(mask771,:);
            
            FScore=abs((mean(P1,1)-mean(P2,1))./sqrt(std(P1,0,1).^2+(std(P2,0,1).^2)));
            
            FSscore2d(:,:,file)=reshape(FScore,[size(PSD,2),size(PSD,3)]);
            
            %%        
        
            clear h
        end
    end

    figure

    x = [2:2:48];
    y = [1:16];

    for i=1:size(FSscore2d,3)

        subplot(1,size(FSscore2d,3),i)
        imagesc(x,y,FSscore2d(:,:,i)',[0 1])
        title(['subject ',subjdir,' file', num2str(i) ])
        xlabel('frequencies [Hz]')
        ylabel('channel')
        
        colorbar
    end

    figure

    subplot(1,2,1)
    imagesc(x,y,prod(FSscore2d,3)',[0 1])
    title(['prod subject ',subjdir])
    xlabel('frequencies [Hz]')
    ylabel('channel')
    colorbar

    subplot(1,2,2 )
    imagesc(x,y,sum(FSscore2d,3)',[0 3])
    title(['sum subject ',subjdir])
    xlabel('frequencies [Hz]')
    ylabel('channel')
    colorbar

    clear FSscore2d
    pause

end

%%

% selected_feat=[137 138 139 nan;
%     23  103 167 183;
%     103 167 nan nan;
%     101 102 165 166;
%     87 102 103 183;
%     102 103 118 119;
%     nan nan nan nan;
%     118 182 246 nan]; %manually

selected_feat=[201 202 203 nan;
    31  151 247 271;
    151 247 nan nan;
    149 150 245 246;
    127 150 151 271;
    150 151 174 175;
    nan nan nan nan;
    174 270 366 nan]; %manually

save(fullfile('results','features'),'selected_feat') 
