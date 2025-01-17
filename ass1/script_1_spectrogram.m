%% THIS SCRIPT LOADS ALL THE .gdf FILES COMPUTES THE SPECTROGRAM AND LOADS 
%% THE RESULTS INTO THE 'results\PSD' FOLDER
clc
close all
clear all

addpath("data\")
addpath("utils\")

load("laplacian16.mat")
%%

datadir='data';

SUBJ=dir(datadir);
SUBJ(1:2)=[];

wlength=0.5;    %externalwindow length
pshift=0.25;    %internalwindow shift
wshift= 0.0625; %externalwindow shift
mlength=1;      %seconds

freqmin=1;      %Hz
freqmax=49;     %Hz

wincov='backward';

for subject=1:length(SUBJ)

    subjdir=SUBJ(subject).name;
    
    D=dir(fullfile(datadir,subjdir));
    D(1:2)=[];

    for file=1:length(D)

        if length(D(file).name)>10 
        
            filepath=fullfile(pwd,datadir,subjdir,D(file).name);
    
            [ S,h0 ] = sload( filepath );        %file extraction    
            S=S(:,1:end-1)*lap;                 %application of the laplacian mask
    
            [PSD,f] = proc_spectrogram(S,wlength,wshift,pshift,h0.SampleRate,mlength);
            
            %select only meaningfull frequencies
            PSD=PSD(:,(f>freqmin & f<freqmax),:);
            f=f(f>freqmin & f<freqmax);
    
            %recompute h.EVENT.POS & .DUR with respect to PSD window indexes
            h.EVENT.POS = proc_pos2win(h0.EVENT.POS, wshift*h0.SampleRate, wincov, wlength*h0.SampleRate);
            h.EVENT.DUR = floor(h0.EVENT.DUR/h0.SampleRate/wshift);
            h.EVENT.TYP = h0.EVENT.TYP;
    
            mkdir(fullfile('results','PSD',D(file).name(1:3)))
            save(fullfile('results','PSD',D(file).name(1:3),[D(file).name(1:end-4) '.mat']),'PSD','f','h')
       
    
        end
    end
end
