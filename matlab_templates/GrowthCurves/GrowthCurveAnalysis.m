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
filename='07292020_growthCurve';%B. subtilis with nucleic acid stains
dirname=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/growthCurves/' filename '.xlsx'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlRange='B50:EJ145';
timeRange='B48:EJ48';
wellRange='A50:A145';

od=xlsread(dirname, xlRange);
t=xlsread(dirname, timeRange);
[data, wells]=xlsread(dirname, wellRange); %why does this need to be formatted differently?

T=length(t);
nwells=96;
condlab={'WT', 'dSigM', 'Control'};

array1=[1:6];
cond1=[1:6];

array2=[7:12];
cond2=[7:12];

control=[73:77,79:96];

for i=1:5
    array1New=array1 + (12*i)
    array2New=array2 + (12*i)
    cond1=[cond1, array1New]
    cond2=[cond2, array2New]
   
end

wellInd={cond1, control;
    cond2, control};

[sod,~]=size(od);
figure
hold on
cmap=jet;
for i=1:sod
    c=round((64/sod)*i);
    col=cmap(c,:);
    plot(t/3600,od(i,:),'Color',col)
end
xlabel('time (h)')
ylabel('absorbance (A.U.)')
fig2pretty

save([filename '_GC.mat'])