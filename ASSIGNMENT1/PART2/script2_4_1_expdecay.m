%% THIS SCRIPT LOADS ALL THE .mat OFFLINE FILES AND THE CLASSIFIER FOR 
%% EVERY SUBJECT THEN USES THE DATA TO PERFORM A GRID SEARCHES FOR THE GLOBAL
%% OPTIMAL ACCUMULATION FRAMEWORK DECAY PARAMETER

clc
close all
clear all

%%
addpath("utils\")

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
    
        if D(file).name(22)=='f'   %ony offline files  
            filepath=fullfile(pwd,datadir,subjdir,D(file).name);
    
            if isempty(H.EVENT) %only on the first iteration
                load( filepath ); 
                [ T,F,A,CF,X ] = label_vector( PSD, h ,0);  
                H.EVENT.POS=h.EVENT.POS;
                H.EVENT.TYP=h.EVENT.TYP;
                H.EVENT.DUR=h.EVENT.DUR;
                S=PSD;
    
            else    %altre iterazioni
                load( filepath );
                [ Tk,Fk,Ak,CFk,Xk ] = label_vector( PSD, h , 0);  

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

    Model=load(fullfile('results','MODEL',D(file).name(1:3)));

    %% build the masks for the two tasks
    win=1:size(S,1);
  
    trials=max(T); 
    for i=1:trials
        ck=unique(A(T==i));
        Ck(i)=ck(ck>0);
    end
    
    % for every window the mask tell if it belongs to a 773 0r 771 type of
    % task (it can also be neither of them
    mask773=(A==773) | (CF>0 & sum(T==find(Ck==773),2));    %problema T è più lungo del vettore delle posizioni
    mask771=(A==771) | (CF>0 & sum(T==find(Ck==771),2));
    
    Sfeat=reshape(S,[size(S,1),size(S,2)*size(S,3)]);
    continuousFeedbackMask=mask771 | mask773;

    %load the subject specific features
    subj_feat=selected_feat(subject,~isnan(selected_feat(subject,:)));
    
    %% test the model
    
    [Gk, pp] = predict(Model.Model, Sfeat(continuousFeedbackMask, subj_feat));

    %% single sample accuracy computation

    groundTruth = 773*mask773 + 771*mask771;
    groundTruth = groundTruth(continuousFeedbackMask);

    [SS_ACCmatrix,classes]=confusionmat(groundTruth,Gk);
%     figure
%     confusionchart(SS_ACCmatrix,'Title',['subj ',subjdir,' OFFLINE single sample acc'],'Normalization','row-normalized')
    
    %% apply exponential smoothing
    

%% apply exponential smoothing
   
    th=0.8;
    Tmasked=T(continuousFeedbackMask);
    
    for iAlpha=1:19
        
        alpha= 0.79+0.01*iAlpha;
        D=0.5*ones(size(pp));
        T0=Tmasked(1);
        trialDecision=zeros(trials,1);   
        timeForDecision=zeros(trials,1);
        time0=0;


        for iW=2:size(pp,1)
            
            trialDecision(T0)=771*(D(iW-1,1)>0.5)+773 *(D(iW-1,1)<0.5);
            time0=time0+1;
    
            if Tmasked(iW-1)~=T0 %all'inizio di ogni trial resetto il valore di D
                D(iW-1,:) = [0.5 0.5];           
                T0=Tmasked(iW-1);
                time0=0;
            else 
                if sum(D(iW-1,:)<th)==2 %control if the decisionboundary has been crossed
                    D(iW,:) = D(iW-1,:) * alpha + pp(iW,:) * (1-alpha); 
                else 
                    D(iW,:) = D(iW-1,:);
                    if timeForDecision(T0)==0
                        timeForDecision(T0)=0.0625*time0;
                    end 
                end
            end
        end
        T_ACCmatrix(:,:,iAlpha)=confusionmat(Ck(timeForDecision>0),trialDecision(timeForDecision>0));
        AVGtime(iAlpha)=mean(timeForDecision(timeForDecision>0));
        NUMdecisions(iAlpha)=sum(timeForDecision>0);
    end
    
    %%

    figure
    alpha=0.8:0.01:0.98;
    
    chosen_alpha=0.94;

    subplot(131)
    plot(alpha,AVGtime,'LineWidth',2)
    xline(chosen_alpha,'k','LineWidth',2)
    ylim([0,7])
    xlabel('decay rate []')
    ylabel('time [s]')
    title('OFFLINE avg decision time')
    grid on

    subplot(132)
    hold on
    plot(alpha,NUMdecisions/T0,'B','LineWidth',2)
    plot(alpha,squeeze(T_ACCmatrix(1,1,:)+T_ACCmatrix(2,2,:))/T0,'r','LineWidth',2)
    plot(alpha,squeeze(T_ACCmatrix(1,1,:)+T_ACCmatrix(2,2,:))./squeeze(sum(T_ACCmatrix,[1,2])),'g','LineWidth',2)
    yline((SS_ACCmatrix(1,1)+SS_ACCmatrix(2,2))/sum(SS_ACCmatrix,'all'),'.-g','LineWidth',2)
    xline(chosen_alpha,'k','LineWidth',2)
    ylim([0,1])
    legend({'fract of decisions','fract of correct dec','model accuracy','sing samp acc','optimal \alpha'},'Location','south','FontSize',9 )
    title(['subj ',subjdir,' OFFLINE'])
    xlabel('decay rate []')
    grid on
    
    subplot(133)
    hold on
    plot(alpha,squeeze(T_ACCmatrix(1,1,:))./squeeze(T_ACCmatrix(1,1,:)+T_ACCmatrix(1,2,:)),'r','LineWidth',2)
    plot(alpha,squeeze(T_ACCmatrix(2,2,:))./squeeze(T_ACCmatrix(2,2,:)+T_ACCmatrix(2,1,:)),'b','LineWidth',2)
    yline(SS_ACCmatrix(1,1)/(SS_ACCmatrix(1,1)+SS_ACCmatrix(1,2)),'-.r','LineWidth',2)
    yline(SS_ACCmatrix(2,2)/(SS_ACCmatrix(2,2)+SS_ACCmatrix(2,1)),'-.b','LineWidth',2)
    xline(chosen_alpha,'k','LineWidth',2)
    ylim([0,1])
    legend({['class' num2str(classes(1)) ' acc'],['class' num2str(classes(2)) ' acc'],['sing samp ' num2str(classes(1)) ' acc'],['sing  samp' num2str(classes(2)) ' acc'],'optimal \alpha'},'Location','south','FontSize',9)
    title(['subj ',subjdir,' OFFLINE'])
    xlabel('decay rate []') 
    grid on


    %%
    pause

end