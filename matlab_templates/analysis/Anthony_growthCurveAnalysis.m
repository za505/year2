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
% saveas(gcf, 'figure25.png')
% saveas(gcf, 'figure25.fig')

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
% saveas(gcf, 'figure26.png')
% saveas(gcf, 'figure26.fig')

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
% saveas(gcf, 'figure27.png')
% saveas(gcf, 'figure27.fig')
