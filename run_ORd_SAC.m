function [output1 output2 output3]=run_ORd_SAC(CL,prepace,beats,Lambda_seq,par,ISAC_factor,tstretch,trelax,nb)


count_lambda=1;

for Lambda=Lambda_seq 
    count=1;
    load total_prepaced500_ToRORdSAC_35species_09ikbfactor
    X0=total(end,2:44);
    options=[];

%%
%run these lnes without loading the data.mat (fill X0 vector)
%% %X0 is the vector for initial sconditions for state variables
% v=-88.8691566357934;
% nai=12.0996647655188;
% nass=12.1000028563765;
% ki=142.412524737626;
% kss=142.412481425842;
% cai=7.45541572746214e-05;
% cass=6.50418928341426e-05;
% cansr=1.53037019085812;
% cajsr=1.52803094224238;
% m=0.000787657400526199;
% hp=0.674096901201792;
% h=0.830658198588696;
% j=0.830466744399495;
% jp=0.830093612199637;
% mL=0.000159670117055769;
% hL=0.528261721740178;
% hLp=0.288775833197764;
% a=0.000944249645410894;
% iF=0.999616956857814;
% iS=0.593680589620082;
% ap=0.000481107253796778;
% iFp=0.999616964658062;
% iSp=0.654092074678260;
% d=8.86091322819384e-29;
% ff=0.999999992783113;
% fs=0.938965241412012;
% fcaf=0.999999992783179;
% fcas=0.999900458262832;
% jca=0.999977476316330;
% nca=0.000492094765239740;
% nca_i=0.000833711885764158;
% ffp=0.999999992566681;
% fcafp=0.999999992766279;
% xs1=0.247156543918935;
% xs2=0.000175017075236424;
% Jrel_np=3.90843796133124e-24;
% CaMKt=0.0110752904836162;
% t_ikr_c0=0.998073652444028;
% t_ikr_c1=0.000844745297078649;
% t_ikr_c2=0.000698171876592920;
% t_ikr_o=0.000370404872169913;
% t_ikr_i=1.30239063420973e-05;
% Jrel_p=-1.88428892080206e-22;


%X0=[v nai nass ki kss cai cass cansr cajsr m hp h j jp mL hL hLp a iF iS ap iFp iSp,...
%        d ff fs fcaf fcas jca nca nca_i ffp fcafp xs1 xs2 Jrel_np CaMKt,...
%        t_ikr_c0 t_ikr_c1 t_ikr_c2 t_ikr_o t_ikr_i Jrel_p]';
  
   
% create memory for some variables
Xsim=[];
timeCL=[0];
timesim=[0];
Xsim2=[];
timeCL2=[0];
timesim2=[0];
Xprepace=[];
timeprepace=[0];

%% prepacing for multiple beats

for n=[1:prepace]
    [time X]=ode15s(@model_SAC,[0 CL],X0,options,1,1,par,ISAC_factor,0,0,0);
    X0=X(size(X,1),:);
    Xprepace=[Xprepace;X];
    timeprepace=[timeprepace;(time+timeprepace(end))];
    if n==1
        timeprepace=timeprepace(2:end);
    end
    n;
end
%%

for n=[1:beats]
    if ismember(n,nb) %apply stretch only during the stretched beat
        [time1 X1]=ode15s(@model_SAC,[0 tstretch],X0,options,1,Lambda,par,ISAC_factor,tstretch,trelax,0);
        X0=X1(size(X1,1),:);
        [time2 X2]=ode15s(@model_SAC,[tstretch trelax],X0,options,1,Lambda,par,ISAC_factor,tstretch,trelax,0);
        X0=X2(size(X2,1),:);
        [time3 X3]=ode15s(@model_SAC,[trelax CL],X0,options,1,Lambda,par,ISAC_factor,tstretch,trelax,0);
        
    else
        [time X]=ode15s(@model_SAC,[0 CL],X0,options,1,1,par,ISAC_factor,0,0,0);
        
    end
    
    if ismember(n,nb)&&tstretch~=0 
        X=[X1;X2;X3];
        time=[time1;time2;time3];
    end
    %
    X0=X(size(X,1),:);
    Xsim=[Xsim;X];
    timeCL=[timeCL;time];
    timesim=[timesim;(time+timesim(end))];
    if n==1
        timeCL=timeCL(2:end);
        timesim=timesim(2:end);%remove first 0 in timesim array
    end
    n; %output beat number to the screen to monitor runtime progress
end


%%%%%%%%%%%%%% dependent variables -> prepacing
Yprepace=[];
for i=[1:size(Xprepace,1)];
    IsJs=model_SAC(timeprepace(i),Xprepace(i,:),0,1,par,ISAC_factor,0,0,0); 
    Yprepace=[Yprepace;IsJs];
end
prepace=horzcat(timeprepace,Xprepace,Yprepace);

%%
%%%%%%%%%%%%%% dependent variables -> SIMULATION
Ysim=[];
cc=0; %count every beginning of the beat
for i=[1:size(Xsim,1)];
    if timeCL(i)==0
        cc=cc+1;
    end
    if ismember(cc,nb) %apply stretch at the stretched beat
        IsJs=model_SAC(timeCL(i),Xsim(i,:),0,Lambda,par,ISAC_factor,tstretch,trelax,0);
    else
        IsJs=model_SAC(timeCL(i),Xsim(i,:),0,1,par,ISAC_factor,0,0,0);
    end
    Ysim=[Ysim;IsJs];
end
sim=horzcat(timesim,Xsim,Ysim);
%total=vertcat(prepace,sim);
 total=sim;

%%%%%%%%%%%%%% name state variables
time=total(:,1);
v=total(:,2);
nai=total(:,3);
nass=total(:,4);
ki=total(:,5);
kss=total(:,6);
cai=total(:,7);
cass=total(:,8);
cansr=total(:,9);
cajsr=total(:,10);
m=total(:,11);
hp=total(:,12);
h=total(:,13);
j=total(:,14);
jp=total(:,15);
mL=total(:,16);
hL=total(:,17);
hLp=total(:,18);
a=total(:,19);
iF=total(:,20);
iS=total(:,21);
ap=total(:,22);
iFp=total(:,23);
iSp=total(:,24);
d=total(:,25);
ff=total(:,26);
fs=total(:,27);
fcaf=total(:,28);
fcas=total(:,29);
jca=total(:,30);
nca=total(:,31);
nca_i=total(:,32);
ffp=total(:,33);
fcafp=total(:,34);
xs1=total(:,35);
xs2=total(:,36);
Jrel_np=total(:,37);
CaMKt=total(:,38);
t_ikr_c0=total(:,39);
t_ikr_c1=total(:,40);
t_ikr_c2=total(:,41);
t_ikr_o=total(:,42);
dt_ikr_i=total(:,43);
dJrel_p=total(:,44);
%
%       %%%%%%%%%%%%%% name dependent variables
INa=total(:,45);
INaL=total(:,46);
Ito=total(:,47);
ICaL=total(:,48);
IKr=total(:,49);
IKs=total(:,50);
IK1=total(:,51);
INaCa_i=total(:,52);
INaCa_ss=total(:,53);
INaK=total(:,54);
IKb=total(:,55);
INab=total(:,56);
ICab=total(:,57);
IpCa=total(:,58);
Jdiff=total(:,59);
JdiffNa=total(:,60);
JdiffK=total(:,61);
Jup=total(:,62);
Jleak=total(:,63);
Jtr=total(:,64);
Jrel=total(:,65);
CaMKa=total(:,66);
Istim=total(:,67);


INsNa=total(:,68);
INsK=total(:,69);
INs=total(:,70);
IKo=total(:,71);
strain=total(:,72);

fINap=total(:,73);
fINaLp=total(:,74);
fICaLp=total(:,75);
fJrelp=total(:,76);
fJupp=total(:,77);
cajsr=total(:,78);
cansr=total(:,79);
PhiCaL_ss=total(:,80);
ICaL_i=total(:,81);
I_ClCa=total(:,82);
I_Clbk=total(:,83);
ICaL_tot=total(:,84);

output1=v;
output2=time;
%

%%Plot

figure(1)
subplot(12,1,[1,2,3,4,5,6,7])
plot(time,v,'LineWidth',2.5), hold on, ylabel({'Membrane';'potential (mV)'}),xlabel('Time (ms)')
set(gca,'linewidth',2)
%set(gca,'XTickLabel',[]);
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize=10;
box off
figure(1)
subplot(12,1,[10,11,12])
plot((time),strain*100,'LineWidth',2.5,'Color','k'),  hold on, ylabel({'Stretch (%)'}),xlabel('Time (ms)')
ylim([0 50])
set(gca,'linewidth',2)
box off
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize=10;


clear INa INaL Ito ICaL IKr IKs IK1 INaCa_i INaCa_ss INaK IKb INab ICab IpCa Jdiff JdiffNa JdiffK Jup Jleak Jtr Jrel CaMKa Istim

clear INsNa INsK INs IKo


end
count_lambda = count_lambda+1;
end

