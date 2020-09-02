%BacTrackAvg.m
%Loads data from BacTrack and generates plots for average 
%growth rate


%INSTRUCTIONS FOR USE:
%INPUT:
%basename: name of the image stack.
%dirname:the full pathname of the directory where you saved the image
%        stack.
%OUTPUT:
%lscale: microscope calibration in microns per pixels.
%sm: width of the Gaussian filter used in edge finder equals sm*sqrt(2).
%minL: minimum length of cells;
%minW: minimum width of cells;
%maxW: maximum width of cells;
%T: number of time points.
%time: vector of length T with time points.
%tmid: vector of length T-1 with interstitial time points.
%ncells: number of individual cells tracked.
%lcell: ncells x T matrix of cell lengths.
%wcell: ncells x T matrix of cell widths.
%acell: ncells x T matrix of cell areas
%ew: ncells x T matrix of circumferential strains.
%acell: ncells x T matrix of cell areas.
%v: ncells x T-1 matrix of cell strain rates.
%B: ncells x T cell array with cell contours.
%mlines: ncells x T cell array with cell midlines
%wav: vector of length T with average cell widths.
%wstd: vector of length T with standard deviations of cell widths.
%wste: vector of length T with standard error of cell widths.
%vav: vector of length T-1 with average cell strain rate.
%vstd: vector of length T-1 with standard deviation of strain rates.
%vste: vector of length T-1 with standard error of strain rates.
%avav: vector of length T-1 with average cell areal strain rate.
%avstd: vector of length T-1 with standard deviation of areal strain rates.
%avste: vector of length T-1 with standard error of areal strain rates.
%ndp: vecotr of lenth T-1 with number of data points averaged over.

%Calls on the following m-files:
%norm16bit.m
%polefinder.m
%cellcurvature.m
%metadata.m
%extrema.m
%EffectiveLength.m
%fig2pretty.m
%movingaverage.m

clear
close all

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maindir=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/08112020_hypershock/08112020_WT_hypershock/'];

cd(maindir);
subdir=dir('*colony*');
subdir=struct2cell(subdir);
[~, numFile]=size(subdir);
%condlab=[1, numFile];

cmap=colormap(gray);
%cmap2=colormap(parula);
%cmap3=colormap(hsv);

for i=1:numFile
    basename=char(subdir(1,i));
    savedir=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/08112020_hypershock/08112020_WT_hypershock/' basename '/' basename '_figureErase'];
    cd(savedir);
    
    load([basename '_BTphase.mat'], 'vav', 'tmid')
    
    data{1, i}=vav;
    data{2, i}=tmid;
end

len = cellfun(@length, data)

if all(len==len(1))
    display("all arrays are of equal lengths")
    MAX = max(len(:));
else
    MAX = max(len(:));
    
    for i=1:numFile
        measure = length(data{1,i}(1, :));
        if  measure < MAX
            col= MAX - measure;
            data{1,i}(1, MAX-col+1:MAX)=NaN(1,col);
            data{2,i}(1, MAX-col+1:MAX)=NaN(1,col);
        else
            display("array of equal lengths")
            tmid=data{2,i}(1, :);
        end 
    end
end

for i=1:numFile
    growthRate(i, :)=data{1,i}(1, :).*3600;
    preShock(i, :)=growthRate(i,1:(MAX/2));
    postShock(i,:)=growthRate(i,(MAX/2+1):MAX);   
end

preMean=mean(preShock);
postMean=mean(postShock);

h1=figure, hold on
for i=1:numFile
    condlab=['Colony ' num2str(i)];
    plot(tmid, growthRate(i, :), 'LineWidth', 0.5,'Color', cmap(256-(i*30), :), 'Display', condlab)
    %title('Growth Rate vs. Time')
    %xlabel('Time (s)')
    %ylabel(strcat('Growth Rate (s^{-1})')) %Shouldn't this be h^-1?
    %fig2pretty
    %savefig([basename,'_growthRate.fig'])
    
   %{
    figure, hold on
    for i=1:numFile
    condlab=['Colony ' num2str(i)]
    plot(tmid(1, 1:(numMid/2)),preShock(i, :), 'LineWidth', 0.5,'Color', cmap2(256-(i*30), :), 'Display', condlab)
    plot(tmid(1,(numMid/2+1):numMid),postShock(i, :), 'LineWidth', 0.5, 'Color', cmap3(256-(i*10), :), 'Display', condlab)
    end
    %}
%dpts=size
end
plot(tmid(1, 1:(MAX/2)),preMean,'LineWidth', 3, 'Color', 'blue', 'Display', 'Pre-shock Avg')
plot(tmid(1,(MAX/2+1):MAX),postMean,'LineWidth', 3, 'Color', 'red', 'Display', 'Post-shock Avg')
%condlab(numFile + 1)= 'Average Pre-shock';
%condlab(numFile + 2)= 'Average Post-shock';
title('Average Growth Rate vs. Time')
xlabel('Time (s)')
ylabel('Growth Rate (s^{-1})')
fig2pretty
legend('Location', 'southeast')
%legend(condlab)
%legend(condlab, Average)
hold off
figure(h1)
cd(maindir);
savefig('growthRate.fig')
saveas(gcf, 'growthRate.png')
%}