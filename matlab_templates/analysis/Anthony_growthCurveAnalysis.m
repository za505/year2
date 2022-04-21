%Growth Curve annd growth rate code with standard deviation

clear, close all

fullpath=('/Users/zarina/Downloads/NYU/Year3_2022_Spring/growthCurves');
dirsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/figures/';

okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

cd(fullpath)

ER300 = [2, 5, 8, 11, 14];
labels = {'LB', 'glucose', 'sorbitol', 'sucrose', 'glycerol'};

load('04052022_growthCurve.mat')

cd(dirsave)

v = 1;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
subplot(1, 2, 1)
for i=ER300
    plot(time2,eOD_smooth(i, :)*3600, 'Color', okabeIto{v}), hold on
    v = v + 1;
end
legend(labels)
title('Growth Rate')
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
%xline(time2(16), '--k')
fig2pretty
% pause
% saveas(gcf, 'figure25.png')
% saveas(gcf, 'figure25.fig')

v = 1;
%figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
subplot(1, 2, 2)
for i=ER300
    plot(time3,OD_norm_smooth(i, :), 'Color', okabeIto{v}), hold on
    v = v + 1;
end
title('Growth Curve')
xlabel('Time (h)')
ylabel('OD(AU)')
legend(labels)
%xline(time2(16), '--k')
fig2pretty
pause
saveas(gcf, 'figure25.png')
saveas(gcf, 'figure25.fig')

cd(fullpath)

ER300 = [2, 5, 8, 11, 14];
labels = {'LB', 'glucose', 'sorbitol', 'sucrose', 'glycerol'};

load('04062022_growthCurve.mat')

cd(dirsave)

v = 1;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
subplot(1, 2, 1)
for i=ER300
    plot(time2,eOD_smooth(i, :)*3600, 'Color', okabeIto{v}), hold on
    v = v + 1;
end
title('Growth Rate')
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
%xline(time2(16), '--k')
legend(labels)
fig2pretty
% pause
% saveas(gcf, 'figure27.png')
% saveas(gcf, 'figure27.fig')

v = 1;
%figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
subplot(1, 2, 2)
for i=ER300
    plot(time3,OD_norm_smooth(i, :), 'Color', okabeIto{v}), hold on
    v = v + 1;
end
title('Growth Curve')
xlabel('Time (h)')
ylabel('OD(AU)')
legend(labels)
%xline(time2(16), '--k')
fig2pretty
pause
saveas(gcf, 'figure26.png')
saveas(gcf, 'figure26.fig')

cd(fullpath)
load('04072022_growthCurve_01.mat')

cd(dirsave)

v = 1;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
subplot(1, 2, 1)
for i=ER300
    plot(time2,eOD_smooth(i, :)*3600, 'Color', okabeIto{v}), hold on
    v = v + 1;
end
title('Growth Rate')
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
%xline(time2(16), '--k')
legend(labels)
fig2pretty
% pause
% saveas(gcf, 'figure29.png')
% saveas(gcf, 'figure29.fig')

v = 1;
%figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
subplot(1, 2, 2)
for i=ER300
    plot(time3,OD_norm_smooth(i, :), 'Color', okabeIto{v}), hold on
    v = v + 1;
end
title('Growth Curve')
xlabel('Time (h)')
ylabel('OD(AU)')
legend(labels)
%xline(time2(16), '--k')
fig2pretty
pause
saveas(gcf, 'figure27.png')
saveas(gcf, 'figure27.fig')
%% original
fullpath=('/Users/zarina/Downloads/NYU/Year3_2022_Spring/growthCurves');
filename=('04072022_growthCurve_01.xlsx');

xlRange='B62:EL121';
nwells=60;
T=141;
contamination = [];

%Conditions
media = {'LB', 'RDM glucose', 'RDM sorbitol', 'RDM sucrose', 'RDM glycerol'};
strains = {'ER002', 'ER300', 'ER005'};

layout = reshape(1:60, 3, 20)';
bidx = [1:4:height(layout)];
cidx = setdiff(1:20, bidx);
conditions = mat2cell(layout(cidx, :), ones(1, 15));
blank = mat2cell(layout(bidx, :), ones(1, 5));

condlab = {'LB 1', 'LB 2', 'LB 3', 'RDM glucose 1', 'RDM glucose 2', 'RDM glucose 3', 'RDM sorbitol 1', 'RDM sorbitol 2', 'RDM sorbitol 3', 'RDM sucrose 1', 'RDM sucrose 2', 'RDM sucrose 3', 'RDM glycerol 1', 'RDM glycerol 2', 'RDM glycerol 3'};

%upload file
cd(fullpath)
GCtable=xlsread(filename,xlRange);
OD=GCtable;

%parse data
count = 1;
for i=1:height(conditions)
    
    if isempty(contamination)==0
        blank{count,1} = setdiff(blank{count,1}, contamination);
    end
    
    conditions(i, 2) = blank(count,1);
    
    if mod(i,3)==0
        count = count+1;
    end
               
end

% WellInd={cond1; %first column, cond wells, second column blanks
%   cond2;
%   cond3;
%   cond4;
%   cond5;
%   cond6;  
%   };

WellInd = conditions;

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


% figure
plot(time2,eOD_smooth*3600)
title('Growth Rate')
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
xline(time2(16), '--k')
legend(condlab)
fig2pretty


figure
plot(time3,OD_norm_smooth)
title('Growth Curve')
xlabel('Time (h)')
ylabel('OD(AU)')
legend(condlab)
xline(time2(16), '--k')
fig2pretty

maxGR = nan(5, 2);
maxTime = nan(5, 2);

[maxGR(:, 1), idx] = max(eOD_smooth([1, 4, 7, 10, 13], :)*3600, [], 2);
maxTime(:, 1) = time2(idx') * 60;
[maxGR(:, 2), idx] = max(eOD_smooth([2, 5, 8, 11, 14], :)*3600, [], 2);
maxTime(:, 2) = time2(idx') * 60;

savename = filename(1:length(filename)-5);
save(savename)

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