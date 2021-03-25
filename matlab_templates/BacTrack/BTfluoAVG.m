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
basename=["03172021_Exp2_colony1", "03172021_Exp2_colony2","03172021_Exp2_colony3",];%Name of the image stack, used to save file.
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03172021_analysis/03172021_Exp2'];
channel=['_647'];
recrunch=0;
B=length(basename); %number of main directories to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the time variables
base=char(basename(1))
cd([filename '/' base '/' base channel '/' base '_figures'])
filelist{1}=dir([base '_BTfluo.mat']);
T = cell2mat(struct2cell(load([filelist{1}.name],'T')));
time = cell2mat(struct2cell(load([filelist{1}.name],'time')));

%pre-allocate matricies
Iout=zeros(B, T);
ratio=[];
%intensity_avg=zeros(B,T);
intensity=[];
bgAdj=[];
frameAuto=[];


for b=1:B
    
    base=char(basename(b))
    cd([filename '/' base '/' base channel '/' base '_figures'])
    filelist{b}=dir([base '_BTfluo.mat']);
    
    bgAdj(b) = cell2mat(struct2cell(load([filelist{b}.name],'bgAdj')));
    frameAuto(b) = cell2mat(struct2cell(load([filelist{b}.name],'frameAuto')));
    Iout(b, :) = cell2mat(struct2cell(load([filelist{b}.name],'Iout')));
    ratio(b, :) = cell2mat(struct2cell(load([filelist{b}.name],'ratio')));
    temp=struct2cell(load([filelist{b}.name],'icell_temp'));
    %temp2=struct2cell(load([filelist{b}.name],'icellAvg_temp'));
    
    %subtract background from the intensity
    temp=cell2mat(temp{1,1});
    temp=temp-bgAdj(b);
    
    %now let's calulate the ratio
    ratio_temp=temp./Iout(b, :)
    %ratio_temp(:,1:frameAuto(b))=NaN;
    
    intensity=[intensity; temp];
    ratio=[ratio; ratio_temp];
    %intensity_avg(b,:)=cell2mat(temp2{1,1});
    
end
    
%what intensity do you expect the population average to be?
%exp_popavg = nanmean(intensity_avg);

%remove ratios from frames without dye
ratio(:,1:frameAuto(b))=NaN;

%let's calulcate the population average intensity based on all the cells
intensity_popavg = nanmean(intensity);
intensity_popstd = std(intensity, 'omitnan');

%let's calulcate the population average based on all the cells
ratio_popavg = nanmean(ratio);
ratio_popstd = std(ratio, 'omitnan');


%let's plot the population avereage intensity
figure, hold on
title('Population Average Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
ylim([0 Inf])
ciplot(intensity_popavg - intensity_popstd, intensity_popavg + intensity_popstd, time, [0.75 0.75 1])
%plot(time, exp_popavg)
plot(time, intensity_popavg, '-r')
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647') 
xline(330, '--k', '*PBS + 647 + FSS') 
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')

%let's plot the average ratio too
figure, hold on
title('Intensity/Background vs Time')
xlabel('Time (s)')
ylabel('Intensity/Background')
fig2pretty
ciplot(ratio_popavg - ratio_popstd, ratio_popavg + ratio_popstd, time, [0.75 0.75 1])
plot(time, ratio_popavg, '-r')
xline(210, '--k', '*PBS + 647') 
xline(330, '--k', '*PBS + 647 + FSS') 
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')