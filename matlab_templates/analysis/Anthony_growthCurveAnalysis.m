
%Growth Curve annd growth rate code with standard deviation

clear, close all

dirsave = ['/Users/zarina/Downloads/NYU/Year3_2022_Spring/growthCurves'];

okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

cd(dirsave)

basename = ['04052022'];
strains = {'ER005'};
conditions = {'LB', 'RDM + glucose', 'RDM + sucrose', 'RDM + sorbitol', 'RDM + glycerol'};

filename = [basename '_growthCurve.xlsx'];

xlRange='B60:CV119';
nwells=60;
T=99;
contamination = [];

OD = xlsread(filename,xlRange);

WellInd = {[10:12], [1:3];
    [22:24], [13:15];
    [34:36], [25:27];
    [46:48], [37:39];
    [58:60], [49:51]};

[Ncond ~]=size(WellInd); % # rows = # conditions
wellODi=cell(Ncond,1);
blankODi=cell(Ncond,1);
OD_av=zeros(Ncond,T);
OD_av_bl=zeros(Ncond,T);
OD_norm=zeros(Ncond,T);
std_OD=zeros(Ncond,T);

%get mean OD for blank
%blankODi=OD(blank, :);
%OD_av_bl=mean(blankODi,1);

for i=1:1:Ncond
    wellODi{i,:}=OD(WellInd{i,1},:); %get the OD from the # wells in the first column of the cell array
    blankODi{i,:}=OD(WellInd{i,2},:); %index values in the second column of the cell array
end

for i=1:Ncond
    OD_av(i,:)=mean(wellODi{i},1); %take the average
    OD_av_bl(i,:)=mean(blankODi{i},1);
    OD_norm=(OD_av-OD_av_bl);
    std_OD(i,:)=std(OD_norm,1);
end

OD_norm_smooth=movingaverage(OD_norm,2); %keep smooth number low
OD_norm_smooth(OD_norm_smooth<=0.001)=0.001;

%calculate growth rate for i=1:Ncond
time=([0:T-1]*600); %this is the time in seconds
dOD=(OD_norm_smooth(:,2:1:end))-(OD_norm_smooth(:,1:1:end-1));
eOD=dOD./(OD_norm_smooth(:,1:1:end-1))/600;
eOD_smooth=movingaverage(eOD,2); 
eOD_smooth(eOD_smooth<=0)=0; %any negative growth rates become 0

%calculate the logOD
logOD = log10(OD_norm_smooth);

tmid=(time(2:2:end)+time(1:2:end-1))/2;
tim=tmid./3600;
time2=([0:T-2]*600)/3600;
time3=([0:T-1]*600)/3600; %this is the time in hours

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
subplot(1, 2, 1)
for i=1:length(conditions)
    plot(time2,eOD_smooth(i, :)*3600, 'Color', okabeIto{i}, 'LineWidth', 1.5), hold on
end
title('Growth Rate')
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
legend(conditions)

subplot(1, 2, 2)
for i=1:length(conditions)
    plot(time3,OD_norm_smooth(i, :), 'Color', okabeIto{i}, 'LineWidth', 1.5), hold on
end
title('Growth Curve')
xlabel('Time (h)')
ylabel('OD (AU)')
saveas(gcf, [basename '_growthCurve_ER005.fig'])
saveas(gcf, [basename '_growthCurve_ER005.png'])
save([basename '_growthCurve_ER005.mat'])
