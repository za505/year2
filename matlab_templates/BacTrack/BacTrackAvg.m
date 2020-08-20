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
maindir=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/08062020_hyposhock/08062020_dSigM_hyposhock/'];

cd(maindir);
subdir=dir('*colony*');
subdir=struct2cell(subdir);
[~, numFile]=size(subdir);
%condlab=[1, numFile];

cmap=colormap(lines);
cmap2=colormap(parula);
cmap3=colormap(hsv);

for i=1:numFile
    basename=char(subdir(1,i));
    savedir=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/08062020_hyposhock/08062020_dSigM_hyposhock/' basename '/' basename '_figureErase'];
    cd(savedir);
    
    load([basename '_BTphase.mat'], 'vav', 'tmid')
    
    vavMat(i, :)=vav;
    tmidMat(i, :)=tmid;
    numMid=length(vav);
    
    growthRate(i, :)=vav.*3600;
    preShock(i, :)=growthRate(i,1:(numMid/2));
    postShock(i,:)=growthRate(i,(numMid/2+1):numMid);   
end

figure, hold on
for i=1:numFile
    title('Growth Rate vs. Time')
    condlab=['Colony ' num2str(i)]
    plot(tmid, growthRate(i, :), 'Color', cmap(i, :), 'Display', condlab)
    xlabel('Time (s)')
    ylabel(strcat('Growth Rate (s^{-1})')) %Shouldn't this be h^-1?
    fig2pretty
    %savefig([basename,'_growthRate.fig'])
end
hold off

figure, hold on
for i=1:numFile
    condlab=['Colony ' num2str(i)]
    plot(tmid(1, 1:(numMid/2)),preShock(i, :), 'Color', cmap2(i*30, :), 'Display', condlab)
    plot(tmid(1,(numMid/2+1):numMid),postShock(i, :), 'Color', cmap3(i*10, :), 'Display', condlab)
end
title('Average Growth Rate vs. Time')
%dpts=size
preMean=mean(preShock);
postMean=mean(postShock);
plot(tmid(1, 1:(numMid/2)),preMean,'LineWidth', 1.5, 'Color', 'blue')
plot(tmid(1,(numMid/2+1):numMid),postMean,'LineWidth', 1.5, 'Color', 'red')
%condlab(numFile + 1)= 'Average Pre-shock';
%condlab(numFile + 2)= 'Average Post-shock';
xlabel('Time (s)')
ylabel('Growth Rate (s^{-1})')
fig2pretty
%legend(condlab)
%legend(condlab, Average)
%saveas(gcf, [basename,'_eTraces.png'])
%}
