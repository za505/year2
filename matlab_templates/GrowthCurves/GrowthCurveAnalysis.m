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
filename=('08062020_growthCurve');
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
condlab={'WTsgRNA', 'yhdL', 'JP4', 'JP5', 'JP6', 'JP27', 'JP28', 'JP29'};

%xlRange='B50:EL145';
%nwells=96
%T=141
cond1a=[1:10:30]; 
cond1b=[31:10:60]; 
cond2a=[2:10:30];
cond2b=[32:10:60];
cond3a=[3:10:30];
cond3b=[33:10:60];
cond4a=[4:10:30];
cond4b=[34:10:60];
cond5a=[5:10:30];
cond5b=[35:10:60];
cond6a=[6:10:30];
cond6b=[36:10:60];
cond7a=[7:10:30];
cond7b=[37:10:60];
cond8a=[8:10:30];
cond8b=[38:10:60];

blank1=[9:10:30;39:10:60]; %index for controls
blank2=[10:10:30;40:10:60];

%{
for i=1:3
    array1New=array1 + (10*i) 
    array2New=array2 + (10*i)
    cond1=[cond1, array1New] %calculate index for cond1
    cond2=[cond2, array2New] %calculate index for cond2
   
end
%}

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

WellInd={cond1a,blank1;
    cond1b,blank2;
    cond2a,blank1;
    cond2b,blank2;
    cond3a,blank1;
    cond3b,blank2;
    cond4a,blank1;
    cond4b,blank2;
    cond5a,blank1;
    cond5b,blank2;
    cond6a,blank1;
    cond6b,blank2;
    cond7a,blank1;
    cond7b,blank2;
    cond8a,blank1;
    cond8b,blank2};

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
   
OD_norm1=OD_norm(1:2:end, :);
OD_norm2=OD_norm(2:2:end, :);

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

eODsmooth1=eODsmooth(1:2:16,:);
eODsmooth2=eODsmooth(2:2:16,:);

%a=find(eODsmooth(1, :) == max(eODsmooth(1,:)));
%b=find(eODsmooth(2, :) == max(eODsmooth(2,:)));

%ta = tim(a); 
%tb = tim(b);
[nrow, ncol] = size(eODsmooth1);

for i=1:nrow
    a(i)=find(eODsmooth1(i, :) == max(eODsmooth1(i,:)));
    b(i)=find(eODsmooth2(i, :) == max(eODsmooth2(i,:)));
    ta(i) = tim(a(i)); 
    tb(i) = tim(b(i));
end

figure
plot(tim,eODsmooth1*3600)
xlabel('Time (h)')
ylabel('Growth Rate in LB (h^{-1})')
%xline(ta, '--b')
%xline(tb, '--', 'Color', [0.8500 0.3250 0.0980])
legend(condlab)
fig2pretty
cd(dirname)
saveas(gcf, [filename '_gr1'])

figure
plot(tim,eODsmooth2*3600)
xlabel('Time (h)')
ylabel('Growth Rate in LB + 1M sorbitol (h^{-1})')
%xline(ta, '--b')
%xline(tb, '--', 'Color', [0.8500 0.3250 0.0980])
legend(condlab)
fig2pretty
cd(dirname)
saveas(gcf, [filename '_gr2'])

figure 
plot(time2,OD_norm1)
xlabel('Time (h)')
ylabel('OD(AU) in LB')
legend(condlab)
fig2pretty
saveas(gcf, [filename '_OD1'])

figure 
plot(time2,OD_norm2)
xlabel('Time (h)')
ylabel('OD(AU) in LB + 1M sorbitol')
legend(condlab)
fig2pretty
saveas(gcf, [filename '_OD2'])

cd(dirname)
save([filename '_GCA.mat'])
