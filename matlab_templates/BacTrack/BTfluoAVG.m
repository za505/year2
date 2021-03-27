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
basename=["03252021_Exp2_colony1", "03252021_Exp2_colony2", "03252021_Exp2_colony3", "03252021_Exp2_colony4"];%Name of the image stack, used to save file.
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03252021_analysis/03252021_Exp2'];
channel=['_647'];
%recrunch=1;
B=length(basename);%number of main directories to analyze
mgGradient=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the time variables
base=char(basename(1))
cd([filename '/' base '/' base channel '/' base '_figures'])
filelist{1}=dir([base '_BTfluo.mat']);
T = cell2mat(struct2cell(load([filelist{1}.name],'T')));
time = cell2mat(struct2cell(load([filelist{1}.name],'time')));
frameAuto = cell2mat(struct2cell(load([filelist{1}.name],'frameAuto')));
mgConc = cell2mat(struct2cell(load([filelist{1}.name],'mgConc')));
mgRange = cell2mat(struct2cell(load([filelist{1}.name],'mgRange')));

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

%remove ratios from frames without dye
ratio(:,1:frameAuto)=NaN;

%let's calulcate the population average based on all the cells
avgRatio = mean(ratio, 'omitnan');
stdRatio = std(ratio, 'omitnan');

%change directory
cd(filename);
save(['BTfluoAVG' channel])

%let's plot the population avereage intensity
figure, hold on
title('Population Average Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
ylim([0 Inf])
ciplot(avgIin - stdIin, avgIin + stdIin, time, [0.75 0.75 1])
%plot(time, exp_popavg)
plot(time, avgIin, '-r')
xline(200, '--k', 'PBS + 647') %frame 20-32
xline(330, '--k', 'PBS + 647 + 12 mM Mg') %frame 33-45
xline(460, '--k', 'PBS + 647 + 15 mM Mg') %frame 46-58
xline(590, '--k', 'PBS + 647 + 20 mM Mg') %frame 59-66
saveas(gcf, ['avgIin' channel '.png'])

%let's plot the average ratio too
figure, hold on
title('Intensity/Background vs Time')
xlabel('Time (s)')
ylabel('Intensity/Background')
ylim([0 Inf])
fig2pretty
ciplot(avgRatio - stdRatio, avgRatio + stdRatio, time, [0.75 0.75 1])
plot(time, avgRatio, '-r') 
xline(200, '--k', 'PBS + 647') %frame 20-32
xline(330, '--k', 'PBS + 647 + 12 mM Mg') %frame 33-45
xline(460, '--k', 'PBS + 647 + 15 mM Mg') %frame 46-58
xline(590, '--k', 'PBS + 647 + 20 mM Mg') %frame 59-66
saveas(gcf, ['avgRatio' channel '.png'])

%let's plot the average ratio vs Mg2+ 
figure, hold on
title('Average Intensity/Background vs Mg^{2+} Concentration')
xlabel('Mg^{2+} (mM)')
ylabel('Intensity/Background')
yline(1, '--b')
xlim([-2 22])
xticks(mgRange)
fig2pretty
plot(mgConc, avgRatio, '-r') 
saveas(gcf, ['avgRatio_mgConc' channel '.png'])
