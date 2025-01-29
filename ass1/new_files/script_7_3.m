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

    load(fullfile(pwd,datadir,subjdir,'ONLINE_ACC_expsys.mat'));
    load(fullfile(pwd,datadir,subjdir,'OFFLINE_EXPACC.mat'));
    load(fullfile(pwd,datadir,subjdir,'ONLINE_ACC_dynsys.mat'));
    load(fullfile(pwd,datadir,subjdir,'OFFLINE_DYNACC.mat'));


% T_ACC_EXPOFFLINEmatrix=1
% T_ACC_EXPONLINEmatrix=1
% timeForDecisionEXPOFFLINE=1
% timeForDecisionEXPONLINE=1
% SS_ACC_EXPOFFLINEmatrix=1
% SS_ACC_EXPONLINEmatrix=0

    SS_accuracy_offline(subject)=(SS_ACC_EXPOFFLINEmatrix(1,1)+SS_ACC_EXPOFFLINEmatrix(2,2))/sum(SS_ACC_EXPOFFLINEmatrix,'all');
    SS_accuracy_online(subject)=(SS_ACC_EXPONLINEmatrix(1,1)+SS_ACC_EXPONLINEmatrix(2,2))/sum(SS_ACC_EXPONLINEmatrix,'all');

    T_accuracy_dynoffline(subject)=(T_ACC_DYNOFFLINEmatrix(1,1)+T_ACC_DYNOFFLINEmatrix(2,2))/sum(T_ACC_DYNOFFLINEmatrix,'all');
    T_accuracy_dynonline(subject)=(T_ACC_DYNONLINEmatrix(1,1)+T_ACC_DYNONLINEmatrix(2,2))/sum(T_ACC_DYNONLINEmatrix,'all');
    mean_time_dynoffline(subject)=mean(timeForDecisionDYNOFFLINE(timeForDecisionDYNOFFLINE>0));
    mean_time_dynonline(subject)=mean(timeForDecisionDYNONLINE(timeForDecisionDYNONLINE>0));
    takendecisions_dynoffline(subject)=sum(timeForDecisionDYNOFFLINE>0)/length(timeForDecisionDYNOFFLINE);
    takendecisions_dynonline(subject)=sum(timeForDecisionDYNONLINE>0)/length(timeForDecisionDYNONLINE);

    T_accuracy_expoffline(subject)=(T_ACC_EXPOFFLINEmatrix(1,1)+T_ACC_EXPOFFLINEmatrix(2,2))/sum(T_ACC_EXPOFFLINEmatrix,'all');
    T_accuracy_exponline(subject)=(T_ACC_EXPONLINEmatrix(1,1)+T_ACC_EXPONLINEmatrix(2,2))/sum(T_ACC_EXPONLINEmatrix,'all');
    mean_time_expoffline(subject)=mean(timeForDecisionEXPOFFLINE(timeForDecisionEXPOFFLINE>0));
    mean_time_exponline(subject)=mean(timeForDecisionEXPONLINE(timeForDecisionEXPONLINE>0));
    takendecisions_expoffline(subject)=sum(timeForDecisionEXPOFFLINE>0)/length(timeForDecisionEXPOFFLINE);
    takendecisions_oexpnline(subject)=sum(timeForDecisionEXPONLINE>0)/length(timeForDecisionEXPONLINE);
end

figure

subplot(131)
hold on
[H,P]=ttest(SS_accuracy_offline,SS_accuracy_online);
plot([1 2],[0.55 0.55],'+-k')
boxplot([SS_accuracy_offline' SS_accuracy_online'],'Labels',{'offline','online'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.5,1])
title('single sample ACCURACY')

subplot(132)
hold on
[H,P]=ttest(T_accuracy_expoffline,T_accuracy_exponline);
plot([1 2],[0.55 0.55],'+-k')
boxplot([T_accuracy_expoffline' T_accuracy_exponline'],'Labels',{'offline','online'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.5,1])
title('Trial ACC EXP')

subplot(133)
hold on
[H,P]=ttest(T_accuracy_dynoffline,T_accuracy_dynonline);
plot([1 2],[0.55 0.55],'+-k')
boxplot([T_accuracy_dynoffline' T_accuracy_dynonline'],'Labels',{'offline','online'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.5,1])
title('Trial ACC DYN')


figure

%prima riga exp seconda dyn
subplot(221)
hold on
[H,P]=ttest(mean_time_dynoffline,mean_time_dynonline);
plot([1 2],[0.55 0.55],'+-k')
boxplot([mean_time_dynoffline' mean_time_dynonline'],'Labels',{'offline','online'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0,5])
title('Mean Time For Decision')

subplot(222)
hold on
[H,P]=ttest(takendecisions_dynoffline,takendecisions_dynonline);
plot([1 2],[0.8 0.8],'+-k')
boxplot([takendecisions_dynoffline',takendecisions_dynonline'],'Labels',{'offline','online'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.75,1])
title('percentage of taken decisions')

% subplot(223)
% hold on
% [H,P]=ttest(mean_time_dynoffline,mean_time_dynonline);
% plot([1 2],[0.55 0.55],'+-k')
% boxplot([mean_time_dynoffline' mean_time_dynonline'],'Labels',{'offline','online'})
% legend({['p=' num2str(P)]},"Location","south")
% ylim([0,5])
% title('Mean Time For Decision')

% subplot(224)
% hold on
% [H,P]=ttest(takendecisions_dynoffline,takendecisions_dynonline);
% plot([1 2],[0.8 0.8],'+-k')
% boxplot([takendecisions_dynoffline',takendecisions_dynonline'],'Labels',{'offline','online'})
% legend({['p=' num2str(P)]},"Location","south")
% ylim([0.75,1])
% title('percentage of taken decisions')



figure
subplot(121)
hold on
[H,P]=ttest(SS_accuracy_offline,T_accuracy_dynoffline);
plot([1 2],[0.55 0.55],'+-k')
boxplot([SS_accuracy_offline' T_accuracy_expoffline' T_accuracy_dynoffline'],'Labels',{'single sample','trial exp','trial dyn'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.5,1])
title('offline ACCURACY')

subplot(122)
hold on
[H,P]=ttest(SS_accuracy_online,T_accuracy_dynonline);
plot([1 2],[0.55 0.55],'+-k')
boxplot([SS_accuracy_online'  T_accuracy_exponline' T_accuracy_dynonline'],'Labels',{'single sample','trial exp','trial dyn'})
legend({['p=' num2str(P)]},"Location","south")
ylim([0.5,1])
title('online ACCURACY')


% figure
% subplot(221)
% hold on
% [H,P]=ttest(SS_accuracy_offline,SS_accuracy_online);
% plot([1 2],[0.55 0.55],'+-k')
% boxplot([SS_accuracy_offline' SS_accuracy_online'],'Labels',{'offline','online'})
% legend({['p=' num2str(P)]},"Location","south")
% ylim([0.5,1])
% title('Trial ACCURACY')
% 
% subplot(222)
% hold on
% [H,P]=ttest(T_accuracy_dynoffline,T_accuracy_dynonline);
% plot([1 2],[0.55 0.55],'+-k')
% boxplot([T_accuracy_dynoffline' T_accuracy_dynonline'],'Labels',{'offline','online'})
% legend({['p=' num2str(P)]},"Location","south")
% ylim([0.5,1])
% title('Trial ACCURACY')
% 
% subplot(223)
% hold on
% [H,P]=ttest(mean_time_dynoffline,mean_time_dynonline);
% plot([1 2],[0.55 0.55],'+-k')
% boxplot([mean_time_dynoffline' mean_time_dynonline'],'Labels',{'offline','online'})
% legend({['p=' num2str(P)]},"Location","south")
% ylim([0,5])
% title('Mean Time For Decision')
% 
% subplot(224)
% hold on
% [H,P]=ttest(takendecisions_dynoffline,takendecisions_dynonline);
% plot([1 2],[0.8 0.8],'+-k')
% boxplot([takendecisions_dynoffline',takendecisions_dynonline'],'Labels',{'offline','online'})
% legend({['p=' num2str(P)]},"Location","south")
% ylim([0.75,1])
% title('percentage of taken decisions')
% 
% 
% 
% 
% figure
% subplot(121)
% hold on
% [H,P]=ttest(SS_accuracy_offline,T_accuracy_dynoffline);
% plot([1 2],[0.55 0.55],'+-k')
% boxplot([SS_accuracy_offline' T_accuracy_dynoffline'],'Labels',{'single sample','trial based'})
% legend({['p=' num2str(P)]},"Location","south")
% ylim([0.5,1])
% title('offline ACCURACY')
% 
% subplot(122)
% hold on
% [H,P]=ttest(SS_accuracy_online,T_accuracy_dynonline);
% plot([1 2],[0.55 0.55],'+-k')
% boxplot([SS_accuracy_online' T_accuracy_dynonline'],'Labels',{'single sample','trial based'})
% legend({['p=' num2str(P)]},"Location","south")
% ylim([0.5,1])
% title('online ACCURACY')
% 
% 
