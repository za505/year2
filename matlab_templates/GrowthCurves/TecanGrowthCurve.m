%TecanGrowthCurve.m
%Rico Rojas, 7/3/19
%Calculates and plots growth growth curves from all wells from Tecan M200.
%
%INPUT
%dirname: full path to excel file with plate reader data 
%OUTPUT
%strain:strain rate of cells.

clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT
filename='growthcurve070319';%B. subtilis with nucleic acid stains
dirname=['/Rico/Bacteria Morphogenesis/Data/Growth Curves/' filename '.xlsx'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[data,text,~]=xlsread(dirname);

od=data(40:end,:);
t=data(38,:);
T=length(t);
wells=text(42,4:end);

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