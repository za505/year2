%GrowthCurveAnalysis.m
%Rico Rojas, 7/3/19
%Paola Bardetti, 1/6/19
%Zarina Akbary, 8/2/2020
%Calculates and plots growth growth curves from all wells from Tecan M200.
%Manipulated TecanGrowthCurve.m and InfiniteOriz.m to generate
%GrowthCurveAnalysis.m
%INPUT
%dirname: full path to excel file with plate reader data 
%OUTPUT
%strain:strain rate of cells.

clear, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT
filename=('12242020_growthCurve');
dirname=('/Users/zarina/Downloads/NYU/Year2_2020_Fall/Lab/growthCurves/12242020/');
filepath=[dirname filename '.xlsx'];

xlRange='B46:Y93'; %what is the range of values for the OD data
timeRange='B44:Y44'; %what is the range of time values
wellRange='B46:B93'; %wells are we pulling data from

OD=xlsread(filepath, xlRange); %OD data
t=xlsread(filepath, timeRange); %time data
[data, wells]=xlsread(filepath, wellRange); %why does this need to be formatted differently?

T=length(t);
tscale=T/600; %why divide by 600? probably bc that's the # second in 10 min
nwells=size(wells);

%Now, let's define conditions
condlab={'ampicillin', 'tunicamycin', 'fosfomycin', 'vancomycin', 'moenomycin', 'lysozyme', 'pronase', 'proteinase', 'LB control'};

amp1=[1:3];
amp2=[4:6];
tun1=[7:9];
tun2=[10:12];
fos1=[13:15];
fos2=[16:18];
van1=[19:21];
van2=[22:24];
moe1=[25:27];
moe2=[28:30];
lys1=[31:33];
lys2=[34:36];
prE1=[37:39];
prE2=[40:42];
prK1=[43:45];
prK2=[46:48];
con1=[49:51];
con2=[52:54];

blank=[55, 56, 60];


WellInd={amp1,blank;
    amp2,blank;
    tun1,blank;
    tun2,blank;
    fos1,blank;
    fos2,blank;
    van1,blank;
    van2,blank;
    moe1,blank;
    moe2,blank;
    lys1,blank;
    lys2,blank;
    prE1,blank;
    prE2,blank;
    prK1kl;,blank;
    prK2,blank;
    con1,blank;
    con2,blank;
    };

%Initialize arrays
[Ncond ~]=size(WellInd); 
wellODi=cell(Ncond,1);
blankODi=cell(Ncond,1);
OD_av=zeros(Ncond,T);
OD_av_bl=zeros(Ncond,T);
OD_norm=zeros(Ncond,T);

%Fill in those empty arrays with data from sheet
for i=1:1:Ncond
    wellODi{i,:}=OD(WellInd{i,1},:);
    blankODi{i,:}=OD(WellInd{i,2},:);
    OD_av(i,:)=mean(wellODi{i},1);
    OD_av_bl(i,:)=mean(blankODi{i},1);
    OD_norm=(OD_av-OD_av_bl);
    OD_st=std(OD_norm);
end 
   
OD_norm1=OD_norm(1:3:end, :); %the OD for 'a' conditions
OD_norm2=OD_norm(2:3:end, :); %the OD for 'b' conditions
OD_norm3=OD_norm(3:3:end, :); %the OD for 'c' conditions

%calculate growthrate  for i=1:Ncond

time=([0:T-1]*600); %why multiply by 600?
dOD=(OD_norm(:,2:2:end))-(OD_norm(:,1:2:end-1));
eOD=dOD./(OD_norm(:,1:2:end-1))/600;
eODsmooth=movingaverage(eOD,10); %what is moving average? the derivative?
%eODsmooth= eODsmooth .* 3600;

tmid=(time(2:2:end)+time(1:2:end-1))/2;
tim=tmid./3600;
time2=([0:T-1]*600)/3600;

eODsmooth1=eODsmooth(1:3:end,:); %the OD for 'a' conditions
eODsmooth2=eODsmooth(2:3:end,:); %the OD for 'b' conditions
%eODsmooth3=eODsmooth(3:3:end,:); %the OD for 'c' conditions

%a=find(eODsmooth(1, :) == max(eODsmooth(1,:)));
%b=find(eODsmooth(2, :) == max(eODsmooth(2,:)));

%ta = tim(a); 
%tb = tim(b);
[nrow, ncol] = size(eODsmooth1);

% for i=1:nrow
%     a(i)=find(eODsmooth1(i, :) == max(eODsmooth1(i,:))); %find index of max growth rate
%     b(i)=find(eODsmooth2(i, :) == max(eODsmooth2(i,:)));
%     c(i)=find(eODsmooth3(i, :) == max(eODsmooth3(i,:)));
%     ta(i) = tim(a(i)); %match it with the time index
%     tb(i) = tim(b(i));
%     tc(i) = tim(c(i));
% end

figure
plot(tim,eODsmooth1*3600)
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
hold on
for i=1:length(ta)
    xline(ta(i), '--')
end
fig2pretty
legend(condlab)
cd(dirname)
saveas(gcf, [filename '_gr'])
saveas(gcf, [filename '_gr.png'])

figure 
plot(time2,OD_norm1)
xlabel('Time (h)')
ylabel('OD(AU)')
fig2pretty
hold on
for i=1:length(ta)
    xline(ta(i), '--')
end
legend(condlab, 'Location', 'southeast')
saveas(gcf, [filename '_OD'])
saveas(gcf, [filename '_OD.png'])

%cd(dirname)
save([filename '_GCA.mat'])