%% THIS SCRIPT LOADS ALL THE .mat ONLINE FILES AND THE CLASSIFIER FOR 
%% EVERY SUBJECT THEN USES THE DATA TO EVALUATE THE PERFORMANCE OF THE 
%% CLASSIFIER

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
    
        if D(file).name(22)=='n'   %ony online files  
            filepath=fullfile(pwd,datadir,subjdir,D(file).name);
    
            if isempty(H.EVENT) %only on the first iteration
                load( filepath ); 
                [ T,F,A,CF,X ] = label_vector( PSD, h ,1);  
                H.EVENT.POS=h.EVENT.POS;
                H.EVENT.TYP=h.EVENT.TYP;
                H.EVENT.DUR=h.EVENT.DUR;
                S=PSD;
    
            else    %altre iterazioni
                load( filepath );
                [ Tk,Fk,Ak,CFk,Xk ] = label_vector( PSD, h , 1);  

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

    matrix=confusionmat(groundTruth,Gk);
    confusionchart(matrix,'Title',['subj ',subjdir,' single sample acc'],'Normalization','row-normalized')
    
    %% apply exponential smoothing
    
    D=0.5*ones(size(pp));
    alpha=0.95;
    Tmasked=T(continuousFeedbackMask);
    T0=Tmasked(1);
    th=0.7;
    trialDecision=zeros(trials,1);

    
    for iW=2:size(pp,1)
        
        trialDecision(T0)=771*(D(iW-1,1)>0.5)+773 *(D(iW-1,1)<0.5);
        
        if Tmasked(iW)~=T0 %all'inizio di ogni trial resetto il valore di D
            D(iW,:) = [0.5 0.5];           
            T0=Tmasked(iW);
            
        else 
            if sum(D(iW-1,:)<th)==2 %control if the decisionboundary has been crossed
            D(iW,:) = D(iW-1,:) * alpha + pp(iW,:) * (1-alpha);    
            else 
                D(iW,:) = D(iW-1,:);
            end
        end
    end
    
    matrix=confusionmat(Ck,trialDecision);
    confusionchart(matrix,'Title',['subj ',subjdir,' trial acc'],'Normalization','row-normalized')

    %%
    figure
    winmasked=win(continuousFeedbackMask);
    %plot(winmasked(Tmasked<2),D(Tmasked<2,1))
    plot(winmasked,D(:,1))
    %xline(winmasked(diff(Tmasked)>0))

    %%
    pause

end