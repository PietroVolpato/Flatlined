%% THIS SCRIPT LOADS ALL THE .mat OFFLINE FILES TRAINS THE CLASSIFIER AND
%% SAVES THE RESULTS INTO THE 'results\CLASSIFIER' FOLDER

clc
close all
clear all

%%
addpath("utils\")
mkdir(fullfile('results','MODEL'))
%%
datadir='results\PSD';

SUBJ=dir(datadir);
SUBJ(1:2)=[];
load("results\features.mat")

for subject=1:length(SUBJ)

    subjdir=SUBJ(subject).name;
    
    if strcmp(subjdir , 'aj7')       %we decided not to consider subject aj7
        continue
    end

    D=dir(fullfile(datadir,subjdir));
    D(1:2)=[];
    
    H=struct('EVENT',[]);
    
    for file=1:length(D)

        if strcmp(D(file).name, 'ai7.20180316.102257.offline.mi.mi_bhbf.mat')
            continue
        end
        if strcmp(D(file).name, 'ai7.20180316.104023.offline.mi.mi_bhbf.mat')
            continue
        end

        if D(file).name(22)=='f'   %only offline files  
            filepath=fullfile(pwd,datadir,subjdir,D(file).name);
    
            if isempty(H.EVENT) %only on the first iteration
                load( filepath ); 
                [ T,F,A,CF,X ] = label_vector( PSD, h );  
                H.EVENT.POS=h.EVENT.POS;
                H.EVENT.TYP=h.EVENT.TYP;
                H.EVENT.DUR=h.EVENT.DUR;
                S=PSD;
    
            else    %altre iterazioni
                load( filepath );
                [ Tk,Fk,Ak,CFk,Xk ] = label_vector( PSD, h );  

                % adds the new run datas to the existing structure
                H.EVENT.POS=[H.EVENT.POS;h.EVENT.TYP+length(S(:,1))];
                H.EVENT.TYP=[H.EVENT.TYP;h.EVENT.TYP];
                H.EVENT.DUR=[H.EVENT.DUR;h.EVENT.DUR];
                
                %concates
                S=[S;PSD];   % datas
                T=[T;Tk+max(T)*(Tk>0)]; % trailal lable
                F=[F;Fk];   % fixation
                A=[A;Ak];   % cue
                CF=[CF;CFk];% feedback
                X=[X;Xk];   %hit miss

            end    
        end
    end
    
    %% control plot 
    
    win=1:size(S,1);
    
    %figure used as a test
%     figure
%     hold on
%     plot(win,T,'b+')
%     plot(win,A*0.1)
%     plot(win,CF*0.1,'r')
    
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

    %control plot
%     figure
%     plot(win,mask773)
%     hold on
%     plot(win,mask771)

    %% load subject selected feature

    subj_feat=selected_feat(subject,~isnan(selected_feat(subject,:)));
    
    %%   
    Sfeat=reshape(S,[size(S,1),size(S,2)*size(S,3)]);
    P1=Sfeat(mask773,:);
    P2=Sfeat(mask771,:);
    
    FScore=abs((mean(P1,1)-mean(P2,1))./sqrt(std(P1,0,1).^2+(std(P2,0,1).^2)));
    
    FSscore2d=reshape(FScore,[size(S,2),size(S,3)]);
    
    %%    
    figure
    x = [2:8:48];
    y = {'Fz','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPz','CP2','CP4'};

    subplot(121)
    imagesc(FSscore2d',[0,1])
    title(['sub ',subjdir,'featmap cummulative'])
    xlabel('frequencies [Hz]')
    ylabel('channel')
    xticks(1:4:24)
    xticklabels(x)    
    yticks(1:16)
    yticklabels(y)
    colorbar
    
    subplot(122)
    prova = FSscore2d;
    prova(subj_feat)=100;
    imagesc(prova',[1,100])
    title(['sub ',subjdir,'selected feat'])
    xticks(1:4:24)
    xticklabels(x)    
    yticks(1:16)
    yticklabels(y)
    xlabel('frequencies [Hz]')
    ylabel('channel')

    
    %% train the model
    continuousFeedbackMask=mask771 | mask773;    
    target_values = 771*mask771+773*mask773;
    
    Model = fitcdiscr(Sfeat(continuousFeedbackMask, subj_feat), target_values(continuousFeedbackMask), 'DiscrimType','quadratic');
    
    save(fullfile('results','MODEL',D(file).name(1:3)),'Model')


    %%
    pause
    clear H
end


