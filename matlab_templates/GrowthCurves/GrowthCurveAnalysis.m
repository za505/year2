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
filename=('07312020_growthCurve');
dirname=('/Users/zarina/Downloads/NYU/Lab_2020_Summer/growthCurves/');
filepath=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/growthCurves/' filename '.xlsx'];

xlRange='B51:EM110'; %what is the range of values for the OD data
timeRange='B49:EM49'; %what is the range of time values
wellRange='A51:A110'; %wells are we pulling data from

OD=xlsread(filepath, xlRange); %OD data
t=xlsread(filepath, timeRange); %time data
[data, wells]=xlsread(filepath, wellRange); %why does this need to be formatted differently?

T=length(t);
tscale=T/600; %why divide by 600?
nwells=size(wells);
condlab={'WT', 'dSigM'};

%xlRange='B50:EL145';
%nwells=96
%T=141
array1=[1:5];
cond1=[1:5]; 

array2=[6:10];
cond2=[6:10];

blank=[41:50,51:60]; %index for controls

for i=1:3
    array1New=array1 + (10*i) 
    array2New=array2 + (10*i)
    cond1=[cond1, array1New] %calculate index for cond1
    cond2=[cond2, array2New] %calculate index for cond2
   
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

%Initialize arrays
[Ncond ~]=size(WellInd); 
wellODi=cell(Ncond,1);
blankODi=cell(Ncond,1);
OD_av=zeros(Ncond,T);
OD_av_bl=zeros(Ncond,T);
OD_norm=zeros(Ncond,T);

%Fill in those empty arrays with data from sheet
for i=1:1:Ncond
    wellODi{i,:}=OD(WellInd{i,1},:)
    blankODi{i,:}=OD(WellInd{i,2},:)
    OD_av(i,:)=mean(wellODi{i},1)
    OD_av_bl(i,:)=mean(blankODi{i},1)
    OD_norm=(OD_av-OD_av_bl)
    OD_st=std(OD_norm)
end 
   
%calculate growthrate  for i=1:Ncond

time=([0:T-1]*600) %why multiply by 600?
dOD=(OD_norm(:,2:2:end))-(OD_norm(:,1:2:end-1))
eOD=dOD./(OD_norm(:,1:2:end-1))/600
eODsmooth=movingaverage(eOD,10); %what is moving average? the derivative?
%eODsmooth= eODsmooth .* 3600;

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

%{
dim=size(eODsmooth);

for j=1:dim(1)
    for i=1:dim(2)-1
        dOD(j,i) = eODsmooth(j, i+1) - eODsmooth(j, i);
        dTim(i)  =  tim(i + 1) - tim(i);
        dOT(j,i) = dOD(j,i)/dTim(i); 
    end
end

figure
plot(tim, dOT)
xlabel('Time (h)')
ylabel('Derivative')
legend(condlab)
fig2pretty
%}

a=find(eODsmooth(1, :) == max(eODsmooth(1,:)));
b=find(eODsmooth(2, :) == max(eODsmooth(2,:)));

ta = tim(a); 
tb = tim(b);

figure
plot(tim,eODsmooth*3600)
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
xline(ta, '--b')
xline(tb, '--', 'Color', [0.8500 0.3250 0.0980])
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

cd(dirname)
save([filename '_GCA.mat'])