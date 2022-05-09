%Growth Curve annd growth rate code with standard deviation

clear, close all

okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

%% original
fullpath=('/Users/zarina/Downloads/NYU/Year3_2022_Spring/growthCurves');
filename=('04272022_growthCurve.xlsx');

xlRange='B51:EM110';
nwells=60;
T=142;
contamination = [33];

%Conditions
media = {'LB', '2.5 hr spent', '3.5 hr spent', '4.5 hr spent', '5.5 hr spent'};
strains = {'ER300'};

layout = reshape(1:60, 6, 10)';
bidx = [1:2:height(layout)];
cidx = setdiff(1:10, bidx);
if isempty(contamination)
    conditions = mat2cell(layout(cidx, :), ones(1, 5));
    blank = mat2cell(layout(bidx, :), ones(1, 5));
else
    conditions = cell(length(cidx), 1);
    blank = cell(length(bidx), 1);
    for i=1:length(cidx)
        j = cidx(i);
        wval = setdiff(layout(j, :), contamination);
        conditions{i} = wval;
    end
    
    for i=1:length(bidx)
        j = bidx(i);
        wval = setdiff(layout(j, :), contamination);
        blank{i} = wval;
    end
end

conditions = [conditions, blank];

condlab = {'LB', '2.5 hr spent', '3.5 hr spent', '4.5 hr spent', '5.5 hr spent'};

%upload file
cd(fullpath)
GCtable=xlsread(filename,xlRange);
OD=GCtable;

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
% plot(time2,eOD_smooth*3600)
% title('Growth Rate')
% xlabel('Time (h)')
% ylabel('Growth Rate (h^{-1})')
% xline(time2(16), '--k')
% legend(condlab)
% fig2pretty
% savename = filename(1:length(filename)-5);
% saveas(gcf, [savename '_growthRate.fig'])
% 
% figure
% plot(time3,OD_norm_smooth)
% title('Growth Curve')
% xlabel('Time (h)')
% ylabel('OD(AU)')
% legend(condlab)
% xline(time2(16), '--k')
% fig2pretty
% savename = filename(1:length(filename)-5);
% saveas(gcf, [savename '_growthCurve.fig'])

% maxGR = nan(5, 2);
% maxTime = nan(5, 2);
% 
% [maxGR(:, 1), idx] = max(eOD_smooth([1, 4, 7, 10, 13], :)*3600, [], 2);
% maxTime(:, 1) = time2(idx') * 60;
% [maxGR(:, 2), idx] = max(eOD_smooth([2, 5, 8, 11, 14], :)*3600, [], 2);
% maxTime(:, 2) = time2(idx') * 60;
% figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
% subplot(1, 2, 1)
% for i=1:length(condlab)
%     plot(time2,eOD_smooth(i, :)*3600, 'Color', okabeIto{i}), hold on
% end
% title('Growth Rate')
% xlabel('Time (h)')
% ylabel('Growth Rate (h^{-1})')
% %xline(time2(16), '--k')
% legend(condlab)
% fig2pretty
% % pause
% % saveas(gcf, 'figure29.png')
% % saveas(gcf, 'figure29.fig')
% 
% %figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
% subplot(1, 2, 2)
% for i=1:length(condlab)
%     plot(time3,OD_norm_smooth(i, :), 'Color', okabeIto{i}), hold on
% end
% title('Growth Curve')
% xlabel('Time (h)')
% ylabel('OD(AU)')
% legend(condlab)
% %xline(time2(16), '--k')
% fig2pretty

%% for the thesis update
idx = [1,5];
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
subplot(1, 2, 1)
for i=idx
    plot(time2,eOD_smooth(i, :)*3600, 'Color', okabeIto{i}), hold on
end
title('Growth Rate')
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
%xline(time2(16), '--k')
legend({'LB', 'Spent Media'})
fig2pretty
% pause
% saveas(gcf, 'figure29.png')
% saveas(gcf, 'figure29.fig')

%figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
subplot(1, 2, 2)
for i=idx
    plot(time3,OD_norm_smooth(i, :), 'Color', okabeIto{i}), hold on
end
title('Growth Curve')
xlabel('Time (h)')
ylabel('OD(AU)')
legend({'LB', 'Spent Media'})
%xline(time2(16), '--k')
fig2pretty

% savename = filename(1:length(filename)-5);
% save(savename)

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