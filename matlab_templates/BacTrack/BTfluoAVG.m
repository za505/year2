%BTfluoAVG.m Zarina Akbary, 03/25/21 based on BTfluo by Rico Rojas, updated
%1/21/19 Calculates the average cytoplasmic fluorescence intensity and
%intensity/background ratio from cell tracked with BacTrack.m and BTfluo.m

clear, close all

%INSTRUCTIONS FOR USE:
%Align phase and fluorescent stacks with imagealign.m, remove pillars from
%the phase image with eraseimage.m, and then analyze using BacTrack.m.
%Calculate fluorescent intensities using BTfluo.m. That data from that
%script will be load here.

%INPUT
%basename: name to obtain data from.
%channels: list of directories containing fluorescent image stacks.

%OUTPUT:



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename=["03172021_Exp1_colony1", "03172021_Exp1_colony2"];%Name of the image stack, used to save file.
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03172021_analysis/03172021_Exp1'];
channel=['_FSS'];
recrunch=0;
B=length(basename);%number of main directories to analyze
mgGradient=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==0
    
%load the time variables
base=char(basename(1))
cd([filename '/' base '/' base channel '/' base '_figures'])
filelist{1}=dir([base '_BTfluo.mat']);
T = cell2mat(struct2cell(load([filelist{1}.name],'T')));
time = cell2mat(struct2cell(load([filelist{1}.name],'time')));
frameAuto = cell2mat(struct2cell(load([filelist{1}.name],'frameAuto')));

if mgGradient==1
    mgRange = cell2mat(struct2cell(load([filelist{1}.name],'mgRange')));
end

%pre-allocate matricies
%there we just want to have in one place to see
cellAdj=zeros(B,1);
bgAdj=zeros(B,1);
bgIntensity=zeros(B,T);
cellIntensity=[];
ratio=[];
Iout=zeros(B, T);
Iin=[];

%these are the ones we'll be re-calculating/plotting
avgIin=zeros(1, T);
stdIin=zeros(1, T);
avgRatio=zeros(1, T);
stdRatio=zeros(1, T);

for b=1:B
    
    base=char(basename(b))
    cd([filename '/' base '/' base channel '/' base '_figures'])
    filelist{b}=dir([base '_BTfluo.mat']);
    
    bgAdj(b,1) = cell2mat(struct2cell(load([filelist{b}.name],'bgAdj')));
    cellAdj(b,1) = cell2mat(struct2cell(load([filelist{b}.name],'cellAdj')));
    bgIntensity(b, :) = cell2mat(struct2cell(load([filelist{b}.name],'bgIntensity')));
    Iout(b, :) = cell2mat(struct2cell(load([filelist{b}.name],'Iout')));
    temp=struct2cell(load([filelist{b}.name],'cellIntensity'));
    ratio_temp=struct2cell(load([filelist{b}.name],'ratio'));
    Iin_temp=struct2cell(load([filelist{b}.name],'Iin'));
    
    %let's get all the raw intensities, adj intensities, and ratios in one place
    temp=temp{1,1};
    cellIntensity=[cellIntensity; temp];
    
    ratio_temp=ratio_temp{1,1};
    ratio=[ratio; ratio_temp];
    
    Iin_temp=Iin_temp{1,1};
    Iin=[Iin; Iin_temp];
  
end
    
%let's calulcate the population average intensity based on all the cells
avgIin = mean(Iin, 'omitnan');
stdIin = std(Iin, 'omitnan');

%let's calulcate the population average background based on all the cells
avgIout = mean(Iout, 'omitnan');
stdIout = std(Iout, 'omitnan');

%remove ratios from frames without dye
ratio(:,1:frameAuto)=NaN;

%let's calulcate the population average based on all the cells
avgRatio = mean(ratio, 'omitnan');
stdRatio = std(ratio, 'omitnan');

%now, calculate the Iin/Iout ratio as a function of Mg2+ concentration
avgMg = [];
stdMg = [];

avgMg(:,1) = mean(ratio(:, 19:30), 2, 'omitnan');
stdMg(1) = std(avgMg(:, 1), 'omitnan');

avgMg(:,2) = mean(ratio(:, 31:42), 2, 'omitnan');
stdMg(2) = std(avgMg(:, 2), 'omitnan');

avgMg(:,3) = mean(ratio(:, 43:54), 2, 'omitnan');
stdMg(3) = std(avgMg(:, 3), 'omitnan');

avgMg(:,4) = mean(ratio(:, 55:61), 2, 'omitnan');
stdMg(4) = std(avgMg(:, 4), 'omitnan');

avgMg = mean(avgMg, 'omitnan');

%change directory
cd(filename);
save(['BTfluoAVG' channel])

else
    %change directory
    cd(filename);
    load(['BTfluoAVG' channel])
end

%let's plot the population avereage intensity
figure, hold on
title('Population Average Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
ylim([-1 Inf])
ciplot(avgIin - stdIin, avgIin + stdIin, time, [0.75 0.75 1])
%plot(time, exp_popavg)
plot(time, avgIin, '-r')
xline(60, '--k', '*PBS + 5% detergent') %frame 1-18
xline(190, '--k', '*PBS + 647 + FSS') %frame 19-30
xline(310, '--k', '*PBS + 647 + FSS + 6.66 mM Mg2+') %frame 31-42
xline(430, '--k', '*PBS + 647 + FSS + 12.33 mM Mg2+') %frame 43-54
xline(550, '--k', '*PBS + 647 + FSS + 20 mM Mg2+') %frame 55-61
saveas(gcf, ['avgIin' channel '.png'])

%let's plot the population avereage intensity
figure, hold on
title('Average Background Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
ylim([-1 Inf])
ciplot(avgIout - stdIout, avgIout + stdIout, time, [0.75 0.75 1])
%plot(time, exp_popavg)
plot(time, avgIout, '-r')
xline(60, '--k', '*PBS + 5% detergent') %frame 1-18
xline(190, '--k', '*PBS + 647 + FSS') %frame 19-30
xline(310, '--k', '*PBS + 647 + FSS + 6.66 mM Mg2+') %frame 31-42
xline(430, '--k', '*PBS + 647 + FSS + 12.33 mM Mg2+') %frame 43-54
xline(550, '--k', '*PBS + 647 + FSS + 20 mM Mg2+') %frame 55-61
saveas(gcf, ['avgIout' channel '.png'])

%let's plot the average ratio too
figure, hold on
title('Intensity/Background vs Time')
xlabel('Time (s)')
ylabel('Intensity/Background')
ylim([-1 Inf])
fig2pretty
ciplot(avgRatio - stdRatio, avgRatio + stdRatio, time, [0.75 0.75 1])
plot(time, avgRatio, '-r') 
xline(60, '--k', '*PBS + 5% detergent') %frame 1-18
xline(190, '--k', '*PBS + 647 + FSS') %frame 19-30
xline(310, '--k', '*PBS + 647 + FSS + 6.66 mM Mg2+') %frame 31-42
xline(430, '--k', '*PBS + 647 + FSS + 12.33 mM Mg2+') %frame 43-54
xline(550, '--k', '*PBS + 647 + FSS + 20 mM Mg2+') %frame 55-61
saveas(gcf, ['avgRatio' channel '.png'])

if mgGradient==1
    %let's plot the average ratio vs Mg2+ 
    figure, hold on
    title('Average Intensity/Background vs Mg^{2+} Concentration')
    errorbar(mgRange,avgMg,stdMg, 'both', 'o')
    xlabel('Mg^{2+} concentration (mM)')
    ylabel('Intensity/Background')
    yline(1, '--k')
    xlim([-2 22])
    ylim([0 Inf])
    xticks(mgRange)
    fig2pretty 
    saveas(gcf, ['avgMg' channel '.png'])
end 
