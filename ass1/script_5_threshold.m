%% THIS SCRIPT LOADS ALL THE .mat OFFLINE FILES AND THE CLASSIFIER FOR 
%% EVERY SUBJECT THEN USES THE DATA TO EVALUATE THE PERFORMANCE OF THE 
%% CLASSIFIER ON TRAINING DATA, FURTHERMORE IT SEARCHES FOR THE OPTIMAL 
%% ACCUMULATION FRAMEWORK THRESHOLD PARAMETERS (SUBJECT-SPECIFIC)

clc
close all
clear all

%%
addpath("utils\")

mkdir(fullfile('results','ACC'))
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

    SS_ACC_OFFLINEmatrix=confusionmat(groundTruth,Gk);
    
    %% apply exponential smoothing
    
    
    alpha=0.94;
    Tmasked=T(continuousFeedbackMask);

    th_vec=0.6:0.05:0.95;

    for iTH=1:length(th_vec)
        
        th=th_vec(iTH);
        D=0.5*ones(size(pp));
        T0=Tmasked(1);
        trialDecision=zeros(trials,1);   
        timeForDecision=zeros(trials,1,length(th_vec));
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
                    if timeForDecision(T0,iTH)==0
                        timeForDecision(T0,iTH)=0.0625*time0;
                    end 
                end
            end
        end
        Trial_ACCmatrix(:,:,iTH)=confusionmat(Ck(timeForDecision(:,iTH)>0),trialDecision(timeForDecision(:,iTH)>0));
        AVGtime(iTH)=mean(timeForDecision(timeForDecision(:,iTH)>0,iTH));
        NUMdecisions(iTH)=sum(timeForDecision(:,iTH)>0);
    end
    

    %%
    SUBJ_TH=[0.8 0.8 0.75 0.75 0.75 0.8 nan 0.75];

    figure

    subplot(131)
    plot(th_vec,AVGtime)
    ylim([0,7])
    xlabel('threshold value []')
    ylabel('time [s]')
    xline(SUBJ_TH(subject),'k')
    title('OFFLINE avg decision time')
    grid on

    subplot(132)
    hold on
    plot(th_vec,NUMdecisions/T0,'B')
    plot(th_vec,squeeze(Trial_ACCmatrix(1,1,:)+Trial_ACCmatrix(2,2,:))/T0,'r')
    plot(th_vec,squeeze(Trial_ACCmatrix(1,1,:)+Trial_ACCmatrix(2,2,:))./squeeze(sum(Trial_ACCmatrix,[1,2])),'g')
    ylim([0,1])
    yline((SS_ACC_OFFLINEmatrix(1,1)+SS_ACC_OFFLINEmatrix(2,2))/sum(SS_ACC_OFFLINEmatrix,'all'),'-g')
    xline(SUBJ_TH(subject),'k')
    legend({'fract of decisions','fract of correct dec','model accuracy','sing samp acc','chosen_threshold'},'Location','south')
    title(['subj ',subjdir,' OFFLINE'])
    xlabel('threshold value []')
    grid on
    
    subplot(133)
    hold on
    plot(th_vec,squeeze(Trial_ACCmatrix(1,1,:))./squeeze(Trial_ACCmatrix(1,1,:)+Trial_ACCmatrix(1,2,:)),'r')
    plot(th_vec,squeeze(Trial_ACCmatrix(2,2,:))./squeeze(Trial_ACCmatrix(2,2,:)+Trial_ACCmatrix(2,1,:)),'b')
    yline(SS_ACC_OFFLINEmatrix(1,1)/(SS_ACC_OFFLINEmatrix(1,1)+SS_ACC_OFFLINEmatrix(1,2)),'-.r')
    yline(SS_ACC_OFFLINEmatrix(2,2)/(SS_ACC_OFFLINEmatrix(2,2)+SS_ACC_OFFLINEmatrix(2,1)),'-.b')
    xline(SUBJ_TH(subject),'k')
    ylim([0,1])
    legend({'class0 acc','class1 acc','sing samp 0 acc','sing samp 1 acc','chosen threshold'},'Location','south')
    title(['subj ',subjdir,' OFFLINE'])
    xlabel('threshold value []')
    grid on


    %%
 
    figure
    subplot(121)
    confusionchart(SS_ACC_OFFLINEmatrix,'Title',['subj ',subjdir,' OFFLINE single sample acc'],'Normalization','row-normalized')

    subjITH=abs(th_vec-SUBJ_TH(subject))<eps;
    T_ACC_OFFLINEmatrix=Trial_ACCmatrix(:,:,subjITH);
    timeForDecisionOFFLINE=timeForDecision(:,subjITH);
    subplot(122)
    confusionchart(T_ACC_OFFLINEmatrix,'Title',['subj ',subjdir,' OFFLINE trial acc'],'Normalization','row-normalized')

    mkdir(fullfile('results','ACC',subjdir))
    save(fullfile('results','ACC',subjdir, 'OFFLINE_ACC.mat'),'T_ACC_OFFLINEmatrix','SS_ACC_OFFLINEmatrix','timeForDecisionOFFLINE')
    pause
end

save(fullfile('results','thresholds'),'SUBJ_TH')
