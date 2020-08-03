%TecanGrowthCurve.m
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
filename=('07292020_growthCurve');
dirname=('/Users/zarina/Downloads/NYU/Lab_2020_Summer/growthCurves/');
filepath=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/growthCurves/' filename '.xlsx'];

xlRange='B50:EJ145';
timeRange='B48:EJ48';
wellRange='A50:A145';

OD=xlsread(filepath, xlRange);
t=xlsread(filepath, timeRange);
[data, wells]=xlsread(filepath, wellRange); %why does this need to be formatted differently?

T=length(t);
tscale=T/600;
nwells=size(wells);
condlab={'WT', 'dSigM'};

%xlRange='B50:EL145';
%nwells=96
%T=141
array1=[1:6];
cond1=[1:6];

array2=[7:12];
cond2=[7:12];

blank=[73:77,79:96];

for i=1:5
    array1New=array1 + (12*i)
    array2New=array2 + (12*i)
    cond1=[cond1, array1New]
    cond2=[cond2, array2New]
   
end

%condition1=od(cond1, :);
%condition2=od(cond2, :);
%blank=od(control, :);

% for i=1:T
%     blankAvg(i)=mean(blank(:, i));
%     condition1=condition1(:, i) - blankaAvg(i);
%     condition2=condition2(:, i) - blankaAvg(i);
%     cond1Avg(i)=mean(condition1(:, i));
%     cond2Avg(i)=mean(condition2(:, i));   
% end

WellInd={cond1,blank;
  cond2,blank};

[Ncond ~]=size(WellInd); 
wellODi=cell(Ncond,1);
blankODi=cell(Ncond,1);
OD_av=zeros(Ncond,T);
OD_av_bl=zeros(Ncond,T);
OD_norm=zeros(Ncond,T);

for i=1:1:Ncond
    wellODi{i,:}=OD(WellInd{i,1},:)
    blankODi{i,:}=OD(WellInd{i,2},:)
end 

for i=1:Ncond
    OD_av(i,:)=mean(wellODi{i},1)
    OD_av_bl(i,:)=mean(blankODi{i},1)
    OD_norm=(OD_av-OD_av_bl)
    %OD_st=std(OD_norm)
end
   
%calculate growthrate  for i=1:Ncond

time=([0:T-1]*600)
dOD=(OD_norm(:,2:2:end))-(OD_norm(:,1:2:end-1))
eOD=dOD./(OD_norm(:,1:2:end-1))/600
eODsmooth=movingaverage(eOD,10)
tmid=(time(2:2:end)+time(1:2:end-1))/2
tim=tmid./3600
time2=([0:T-1]*600)/3600

% [sod,~]=size(od);
% figure
% hold on
% cmap=jet;
% for i=1:sod
%     c=round((64/sod)*i);
%     col=cmap(c,:);
%     plot(t/3600,od(i,:),'Color',col)
% end
% xlabel('time (h)')
% ylabel('absorbance (A.U.)')
% fig2pretty
% 
% save([filename '_GC.mat'])

figure
plot(tim,eODsmooth*3600)
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
legend(condlab)
fig2pretty
cd(dirname)
saveas(gcf, [filename '_gr'])
 
figure 
plot(time2,OD_norm)
xlabel('Time (h)')
ylabel('OD(AU)')
legend(condlab)
fig2pretty
saveas(gcf, [filename '_OD'])