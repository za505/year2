%Growth Curve annd growth rate code with standard deviation

clear, close all

filename=('01242022_growthCurve.xlsx');
fullpath=('/Users/zarina/Downloads/NYU/Year3_2022_Spring/growthCurves');
xlRange='B50:DH103';
nwells=54;
T=111;

%Conditions
cond1=[7:12]; % Control
cond2=[13:15]; % vanco1
cond3=[16:18]; % vanco2
cond4=[19:21]; % amp1
cond5=[22:24]; % amp2 
cond6=[25:27]; % cef1
cond7=[28:30]; % cef2 
cond8=[31:33]; % fos1
cond9=[34:36]; % fos2
cond10=[37:39]; % carb1
cond11=[40:42]; % carb2
cond12=[43:45]; % tun1
cond13=[46:48]; % tun2
cond14=[49:51]; % moe1
cond15=[52:54]; % moe2
blank=[1:6];

condlab={'Control', 'vancomycin 1', 'vancomycin 2', 'ampicillin 1', 'ampicillin 2', 'cefsulodin 1', 'cefsulodin 2', 'phosphomycin 1', 'phosphomycin 2', 'carbenicillin 1', 'carbenicillin 2', 'tunicamycin 1', 'tunicamycin 2', 'moenomycin 1', 'moenomycin 2'};

%upload file
cd(fullpath)
GCtable=xlsread(filename,xlRange);
OD=GCtable;

%parse data
WellInd={cond1; %first column, cond wells, second column blanks
  cond2;
  cond3;
  cond4;
  cond5;
  cond6;
  cond7;
  cond8;
  cond9;
  cond10;
  cond11;
  cond12; 
  cond13;
  cond14;
  cond15;   
  };

[Ncond ~]=size(WellInd); % # rows = # conditions
wellODi=cell(Ncond,1);
%blankODi=cell(Ncond,1);
OD_av=zeros(Ncond,T);
%OD_av_bl=zeros(Ncond,T);
OD_norm=zeros(Ncond,T);
std_OD=zeros(Ncond,T);

%get mean OD for blank
blankODi=OD(blank, :);
OD_av_bl=mean(blankODi,1);

for i=1:1:Ncond
    wellODi{i,:}=OD(WellInd{i,1},:); %get the OD from the # wells in the first column of the cell array
    %blankODi{i,:}=OD(WellInd{i,2},:); %index values in the second column of the cell array
end

for i=1:Ncond
    OD_av(i,:)=mean(wellODi{i},1); %take the temporal average
    %OD_av_bl(i,:)=mean(blankODi{i},1);
    OD_norm=(OD_av-OD_av_bl);
    std_OD(i,:)=std(OD_norm,1);
end

OD_norm_smooth=movingaverage(OD_norm,2); %keep smooth number low
OD_norm_smooth(OD_norm_smooth<=0.001)=0.001;

%calculate growthrate for i=1:Ncond
time=([0:T-1]*600);
dOD=(OD_norm_smooth(:,2:1:end))-(OD_norm_smooth(:,1:1:end-1));
eOD=dOD./(OD_norm_smooth(:,1:1:end-1))/600;
eOD_smooth=movingaverage(eOD,5); 
eOD_smooth(eOD_smooth<=0)=0; %any negative growth rates become 0


tmid=(time(2:2:end)+time(1:2:end-1))/2;
tim=tmid./3600;
time2=([0:T-2]*600)/3600;
time3=([0:T-1]*600)/3600;


% figure
plot(time2,eOD_smooth*3600)
title('Growth Rate')
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
legend(condlab)
fig2pretty

figure
plot(time3,OD_norm_smooth)
title('Growth Curve')
xlabel('Time (h)')
ylabel('OD(AU)')
legend(condlab)
fig2pretty
% 
% 
% figure
% %ciplot((OD_norm_smooth-std_OD),(OD_norm_smooth+std_OD),time3,[0.75 0.75 1])
% xlabel('Time (h)')
% ylabel('OD(AU)')
% fig2pretty
% hold on
% plot(time3,OD_norm_smooth,'LineWidth',2)
% title('Growth Curve')
% xlabel('Time (h)')
% ylabel('OD(AU)')
% legend(condlab)

%% Prediction
% [~, idx] = max(OD_norm_smooth(1,:));
% y=OD_norm_smooth(1, 1:idx);
% x=time3(1, 1:idx);
% f=fit(x',y','exp1')
% 
% figure(1)
% plot(f, x,y)
% 
% [~, idx] = max(OD_norm_smooth(2,:));
% y=OD_norm_smooth(2, 1:idx);
% x=time3(1, 1:idx);
% f=fit(x',y','exp1')
% 
% figure(2)
% plot(f, x,y)
%% Function

% function [dilution_factor]=dilPredict(current_OD, incubation_time, OD_curve, time)
% % 
% [~, idx] = max(OD_curve);
% y=OD_curve(1, 1:idx);
% x=time(1, 1:idx);
% f=fit(x',y','exp1');


% end