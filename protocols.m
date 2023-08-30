%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for  "Electrophysiological effects of stretch-activated ion channels: 
% a systematic computational characterization" 
% by Melania Buonocunto, Aurore Lyon, Tammo Delhaas, Jordi Heijman, and Joost Lumens
% DOI: 10.1113/JP284439, The Journal of Physiology 2023
% contact: work: m.buonocunto@maastrichtuniversity.nl 
%          personal: melania.buonoc@gmail.com           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Simulation parameters
par=[0.2203 10.3203 12.8208 1.4110 0.8742]; 
ISAC_factor=[4.47 0.001248]; %change to [0 0] if you want the diseased condition

%% Inputs
CL=1000; %cycle length
prepace=0; 
beats=3; %how many beats display
nb=[2,3]; %which beat to stretch
Lambda_seq=[1.1] %stretch amplitude
tstretch=[600]; %stretch time of application
interval=10; %stretch duration

%% Simulate AP
for Lambda_seq=Lambda_seq
for tstretch=tstretch
    trelax=tstretch+interval; 
    run_ToRORd_SAC(CL,prepace,beats,Lambda_seq,par,ISAC_factor,tstretch,trelax,nb);
end
end
