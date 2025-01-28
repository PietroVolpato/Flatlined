clc
close all
clear all

%%
addpath("utils\")

%%
datadir='results\ACC';

SUBJ=dir(datadir);
SUBJ(1:2)=[];

for subject=1:length(SUBJ)

    subjdir=SUBJ(subject).name;
    load(fullfile(pwd,datadir,subjdir,'ONLINE_ACC.mat'));
    load(fullfile(pwd,datadir,subjdir,'OFFLINE_ACC.mat'));

    T_accuracy_offline(subject)=(T_ACC_OFFLINEmatrix(1,1)+T_ACC_OFFLINEmatrix(2,2))/sum(T_ACC_OFFLINEmatrix,'all');
    T_accuracy_online(subject)=(T_ACC_ONLINEmatrix(1,1)+T_ACC_ONLINEmatrix(2,2))/sum(T_ACC_ONLINEmatrix,'all');
    SS_accuracy_offline(subject)=(SS_ACC_OFFLINEmatrix(1,1)+SS_ACC_OFFLINEmatrix(2,2))/sum(SS_ACC_OFFLINEmatrix,'all');
    SS_accuracy_online(subject)=(SS_ACC_ONLINEmatrix(1,1)+SS_ACC_ONLINEmatrix(2,2))/sum(SS_ACC_ONLINEmatrix,'all');
    mean_time_offline(subject)=mean(timeForDecisionOFFLINE(timeForDecisionOFFLINE>0));
    mean_time_online(subject)=mean(timeForDecisionONLINE(timeForDecisionONLINE>0));
    takendecisions_offline(subject)=sum(timeForDecisionOFFLINE>0)/length(timeForDecisionOFFLINE);
    takendecisions_online(subject)=sum(timeForDecisionONLINE>0)/length(timeForDecisionONLINE);
end

figure
subplot(221)
hold on
[H,P]=ttest(SS_accuracy_offline,SS_accuracy_online);
plot([1 2],[0.55 0.55],'+-k')
boxplot([SS_accuracy_offline' SS_accuracy_online'],'Labels',{'offline','online'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.5,1])
title('Trial ACCURACY')

subplot(222)
hold on
[H,P]=ttest(T_accuracy_offline,T_accuracy_online);
plot([1 2],[0.55 0.55],'+-k')
boxplot([T_accuracy_offline' T_accuracy_online'],'Labels',{'offline','online'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.5,1])
title('Trial ACCURACY')

subplot(223)
hold on
[H,P]=ttest(mean_time_offline,mean_time_online);
plot([1 2],[0.55 0.55],'+-k')
boxplot([mean_time_offline' mean_time_online'],'Labels',{'offline','online'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0,5])
title('Mean Time For Decision')

subplot(224)
hold on
[H,P]=ttest(takendecisions_offline,takendecisions_online);
plot([1 2],[0.8 0.8],'+-k')
boxplot([takendecisions_offline',takendecisions_online'],'Labels',{'offline','online'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.75,1])
title('percentage of taken decisions')




figure
subplot(121)
hold on
[H,P]=ttest(SS_accuracy_offline,T_accuracy_offline);
plot([1 2],[0.55 0.55],'+-k')
boxplot([SS_accuracy_offline' T_accuracy_offline'],'Labels',{'single sample','trial based'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.5,1])
title('offline ACCURACY')

subplot(122)
hold on
[H,P]=ttest(SS_accuracy_online,T_accuracy_online);
plot([1 2],[0.55 0.55],'+-k')
boxplot([SS_accuracy_online' T_accuracy_online'],'Labels',{'single sample','trial based'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.5,1])
title('online ACCURACY')





