%BTfluo.m
%Rico Rojas, updated 1/21/19
%Calculates the average cytoplasmic fluorescence intensity from cell
%tracked with BacTrack.m.  

clear, close all

%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save fluorescent image stacks
%directories by themselves. 

%INPUT
%basename: name to save the results to.
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%cellIntensity=vector with raw intensity values
%bgIntensity=vector with raw background intensity values
%Iin=cellular intensity values adjusted for autofluorescence
%Iout=background autointensity values adjusted for background

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='03262021_Exp1_colony3';
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03262021_analysis/' basename];
channel=[filename '/' basename '_FSS/' basename '_aligned'];
midSwitch=0; %0=all the frames before frameAuto have no dye
frameInitial=1; %this is the FIRST frame that has no dye
frameAuto=17; %this is the LAST frame that has no dye
frameBg=23; %this is the frame that you'll pick the background area from
mg=1; %is there a Mg gradient?
mgRange=[0 6 9 12]; %which concentrations?
mgConc=[repelem(0,frameAuto+12), repelem(6,13), repelem(9,13), repelem(12,9)]; 
recrunch=0; %recrunch=1, just redo the plots, don't recalculate any values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==0

%load fluorescent images
curdir=cd;
cd(channel); 
fluo_directory=dir('*.tif');

%load BacTrack data
load([filename '/' basename '_phase/' basename '_figures/' basename '_BTphase'], 'T', 'time', 'pixels')

%preallocate cells
cellIntensity=[]; %cellular intensities
avgIntensity=zeros(1,T); %population average cellular intensity
stdIntensity=zeros(1,T); %std of population cellular intensity

bgIntensity=zeros(1,T); %background intensity
cellAuto=zeros(height(pixels),frameAuto-frameInitial); %intensity of cells prior to dye perfusion, aka autofluorescence
bgAuto=zeros(1,frameAuto-frameInitial); %intensity of background prior to dye perfusion, aka autofluorescence

%determine region where you'll measure background intensity
imagename=fluo_directory(frameBg).name;
im=imread(imagename);    
[p1, p2]=getBackground(imagename);

%now let's track intensity over time
for t=1:T

    t

    %load the image
    imagename=fluo_directory(t).name;
    im=imread(imagename);
    
    %measure background level
    bglevel = measureBackground(imagename, p1, p2);
    bgIntensity(t)=bglevel; 
    
    %calculate autofluorescence
    if t <= frameAuto & midSwitch==0
        bgAuto(1, t)=measureBackground(imagename, p1, p2);
    elseif (t <= frameAuto | t >= frameInitial) & midSwitch==1
        bgAuto(1, t)=measureBackground(imagename, p1, p2);
    elseif t > frameAuto 
        bgAuto(1, t)=NaN;
    end
    
    %measure cellular intensity
    for j=1:height(pixels)
         
        %calculate intensity
        cellIntensity(j,t)=mean(im(pixels{j,t}));
        
        %calculate autofluorescence
        if t <= frameAuto & midSwitch==0
            cellAuto(j, t)=mean(im(pixels{j,t}));
        elseif (t <= frameAuto | t >= frameInitial) & midSwitch==1
            cellAuto(j, t)=mean(im(pixels{j,t}));
        elseif t > frameAuto 
            cellAuto(j, t)=NaN;
        end
        
    end

end
    
%cellIntensity(cellIntensity==0)=NaN;
%cellAuto(cellAuto==0)=NaN;
        
%take the population average of the autofluorescence 
cellAuto=mean(cellAuto, 'omitnan');

%now take the temporal average of the background and cellular
%autofluorescence
cellAdj=mean(cellAuto, 'all', 'omitnan');
bgAdj = mean(bgAuto, 'all', 'omitnan');

%now subtract the background autofluorescence from background intensity
Iout = bgIntensity - bgAdj;

%now subtract the background autofluorescence from background intensity
Iin = cellIntensity - cellAdj;

%calculate average intensity and standard deviation
avgIntensity = mean(cellIntensity, 'omitnan');
stdIntensity = std(cellIntensity, 'omitnan');

%now, calculate the Iin/Iout ratio
ratio = Iin ./ Iout;
ratio(frameInitial:frameAuto)=NaN;

avgRatio = mean(ratio, 'omitnan');
stdRatio = std(ratio, 'omitnan');

elseif recrunch==1
    load ([filename '/' basename '_FSS/' basename '_figures/' basename '_BTfluo'])
end

%let's change folders to save the plots and variables
cd([filename '/' basename '_FSS/' basename '_figures'])
    
%save the variables
save([basename '_BTfluo'])

%Plot data
%Let's plot cellular intensity first
figure, hold on, title('Intensity vs Time')
for i=1:height(pixels)
    plot(time,Iin(i,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
ylim([0 Inf])
xline(170, '--k', 'PBS + FSS') %frame 18-29
xline(300, '--k', 'PBS + FSS + 6 mM Mg') %frame 30-42
xline(430, '--k', 'PBS + FSS + 9 mM Mg') %frame 43-55
xline(560, '--k', 'PBS + FSS + 12 mM Mg') %frame 56-64
saveas(gcf, [basename '_intensity.png'])

%now population average cell intensity
figure
ciplot(avgIntensity - stdIntensity, avgIntensity + stdIntensity, time, [0.75 0.75 1])
plot(time, avgIntensity,'-r')
title('Average Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
ylim([0 Inf])
xline(170, '--k', 'PBS + FSS') %frame 18-29
xline(300, '--k', 'PBS + FSS + 6 mM Mg') %frame 30-42
xline(430, '--k', 'PBS + FSS + 9 mM Mg') %frame 43-55
xline(560, '--k', 'PBS + FSS + 12 mM Mg') %frame 56-64
saveas(gcf, [basename,'_intensityAvg.png'])

%Plot average background fluorescence
figure, hold on 
title('Background Intensity vs Time')
plot(time,Iout)
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
ylim([0 Inf])
xline(170, '--k', 'PBS + FSS') %frame 18-29
xline(300, '--k', 'PBS + FSS + 6 mM Mg') %frame 30-42
xline(430, '--k', 'PBS + FSS + 9 mM Mg') %frame 43-55
xline(560, '--k', 'PBS + FSS + 12 mM Mg') %frame 56-64
hold off
saveas(gcf, [basename,'_background.png'])

%Plot the Iin/Iout ratio over time
figure, hold on
title('Intensity/Background vs Time')
ciplot(avgRatio - stdRatio, avgRatio + stdRatio, time, [0.75 0.75 1])
plot(time,avgRatio)
xlabel('Time (s)')
ylabel('Intensity/Background')
fig2pretty 
ylim([0 Inf])
xline(170, '--k', 'PBS + FSS') %frame 18-29
xline(300, '--k', 'PBS + FSS + 6 mM Mg') %frame 30-42
xline(430, '--k', 'PBS + FSS + 9 mM Mg') %frame 43-55
xline(560, '--k', 'PBS + FSS + 12 mM Mg') %frame 56-64
hold off
saveas(gcf, [basename,'_ratioTime.png'])

%Plot the Iin/Iout ratio over [Mg2+]
figure, hold on
title('Intensity/Background vs Mg^{2+}')
plot(mgConc,avgRatio)
xlabel('Mg^{2+} concentration (mM)')
ylabel('Intensity/Background')
yline(1, '--k')
xlim([-2 13])
ylim([0 Inf])
xticks(mgRange)
fig2pretty 
hold off
saveas(gcf, [basename,'_ratioMg.png'])

%%%%%Functions
function [p1, p2]=getBackground(imagename)
        
        %Load last image
        %imagename=fluo_directory{i}(t).name;
        im2=imread(imagename);

        %Determine Background
        figure,imshow(im2,[]), hold on, title('Select Background')
        k=waitforbuttonpress;
        set(gcf,'Pointer')
        hold on
        axis manual
        point1=get(gca,'CurrentPoint');
        finalRect=rbbox;
        point2=get(gca,'CurrentPoint');
        point1=point1(1,1:2);
        point2=point2(1,1:2);
        point1(point1<1)=1;
        point2(point2<1)=1;
        p1=min(point1,point2);%Calculate locations
        p2=max(point1,point2);
        offset = abs(point1-point2);%And dimensions
        x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
        y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
        plot(x,y)
        p1=round(p1);
        p2=round(p2);  
end 

function bglevel = measureBackground(imagename, p1, p2)
        
        %Load last image
        %imagename=fluo_directory{i}(t).name;
        im2=imread(imagename);
        
        %Determine background
        backim=im2(p1(2):p2(2),p1(1):p2(1));
        [counts,bins]=imhist(backim);
        [~,binnum]=max(counts);
        maxpos=bins(binnum);
        bglevel=mean(mean(backim));
        
end 