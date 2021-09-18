%Growth Curve annd growth rate code with standard deviation

clear, close all

filename=('09162021_growthCurve.xlsx');
fullpath=('/Users/zarina/Downloads/NYU/Year3_2021_Fall/growthCurves');
xlRange='B51:EL80';
nwells=96;
T=141;

%Conditions
cond1=[1:5:11]; % ER466 LB
blank1=[16:5:26]; % LB blank 
cond2=[2:5:12]; % ER466 BHI
blank2=[17:5:27]; % BHI blank
cond3=[3:5:13]; % ER466 TSB
blank3=[18:5:28]; % TSB blank
cond4=[4:5:14]; % ER451
blank4=[19:5:29]; %JH9 blank
cond5=[5:5:15]; % ER451
blank5=[20:5:30]; %LB blank

condlab={'ER466 LB','ER466 BHI', 'ER466 TSB', 'ER451 JH9', 'ER451 LB'};

%upload file
cd(fullpath)
GCtable=xlsread(filename,xlRange);
OD=GCtable;
tscale=GCtable(1,:)/600;

%parse data
WellInd={cond1,blank1; %first column, cond wells, second column blanks
  cond2,blank2;
  cond3,blank3;
  cond4,blank4;
  cond5,blank5;
  };

[Ncond ~]=size(WellInd); % # rows = # conditions
wellODi=cell(Ncond,1);
blankODi=cell(Ncond,1);
OD_av=zeros(Ncond,T);
OD_av_bl=zeros(Ncond,T);
OD_norm=zeros(Ncond,T);
std_OD=zeros(Ncond,T);

for i=1:1:Ncond
    wellODi{i,:}=OD(WellInd{i,1},:); %get the OD from the # wells in the first column of the cell array
    if i==2
        blankODi{i,:}=repelem(OD(WellInd{i,2},1),T); %index values in the second column of the cell array
    else
        blankODi{i,:}=OD(WellInd{i,2},:); %index values in the second column of the cell array
    end
end

for i=1:Ncond
    OD_av(i,:)=mean(wellODi{i},1); %take the temporal average
    OD_av_bl(i,:)=mean(blankODi{i},1);
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


figure
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


figure
%ciplot((OD_norm_smooth-std_OD),(OD_norm_smooth+std_OD),time3,[0.75 0.75 1])
xlabel('Time (h)')
ylabel('OD(AU)')
fig2pretty
hold on
plot(time3,OD_norm_smooth,'LineWidth',2)
title('Growth Curve')
xlabel('Time (h)')
ylabel('OD(AU)')
legend(condlab)
