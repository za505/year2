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
%icell: Cell array with length equal to the number of fluorescent
        %channels.  Each entry is a matrix (ncellxT) with the fluorescent intensities of each
        %cell, where rows are the cells and columns are time points.
%icell_av:  Cell array with length equal to the number of fluorescent
        %channels.  Each entry is a vector containing the population-
        %average of the single-cell fluorescent intensities.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='03172021_Exp2_colony1';
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03172021_analysis/' basename];
channels={[filename '/' basename '_647/' basename '_full']};
frameAuto=20; %this is the frame that you'll pick the autofluorescence from
frameBack=40; %this is the frame that you'll pick the background area from
recrunch=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==0

curdir=cd;
for i=1:length(channels)
    cd(channels{i}); 
    fluo_directory{i}=dir('*.tif');
end

load([filename '/' basename '_phase/' basename '_figures/' basename '_BTphase'])

%preallocate cells
icell_temp=cell(length(channels),1);
icellAvg_temp=cell(length(channels),1);

icell_adj=cell(length(channels),1);
icellAvg_adj=cell(length(channels),1);

icell_auto=cell(length(channels),1);
icellAvg_auto=cell(length(channels),1);

icell_ratio=cell(length(channels),1);
icellAvg_ratio=cell(length(channels),1);

for i=1:length(channels)
    
    cd(channels{i}); 
    intensity_temp=zeros(size(lcell)); %just a measure of mean intensity
    intensity_adj=zeros(size(lcell)); %mean intensity - bglevel
    intensity_auto=zeros(size(lcell)); %mean intensity - bglevel - autofluorescence
    intensity_ratio=zeros(size(lcell)); %ratio intensity_auto/bglevel
    
    %determine region where you'll measure background intensity
    imagename=fluo_directory{i}(frameBack).name;
    im=imread(imagename);
    
    [p1, p2]=getBackground(imagename);
    
    %measure autointensity
    imagename=fluo_directory{i}(frameAuto).name;
    im=imread(imagename);
    
    for j=1:ncells
        autoFluo=mean(im(pixels{j,frameAuto}));
    end
    
    %now let's track intensity over time
    for t=1:T
        
        t
        
        %load the image
        imagename=fluo_directory{i}(t).name;
        im=imread(imagename);
       
        %measure background level
        bglevel = measureBackground(imagename, p1, p2);
        
        for j=1:ncells
     
            %calculate intensity
            intensity_temp(j,t)=mean(im(pixels{j,t}));
            %subtract the background intensity
            intensity_adj(j,t)=intensity_temp(j,t) - bglevel;
            %subtract the autofluorescence
            intensity_auto(j,t)=intensity_adj(j,t) - autoFluo;
            %calculate the ratio
            intensity_ratio(j,t)= intensity_auto(j,t)/bglevel;
            
        end
        
    end
    
    intensity_temp(intensity_temp==0)=NaN;
    icell_temp{i}=intensity_temp;
    
    intensity_adj(intensity_adj==0)=NaN;
    icell_adj{i}=intensity_adj;
    
    intensity_auto(intensity_auto==0)=NaN;
    icell_auto{i}=intensity_auto;
    
    intensity_ratio(intensity_ratio==0)=NaN;
    icell_ratio{i}=intensity_ratio;
    
end

for i=1:length(channels)
    icellAvg_temp{i}=nanmean(icell_temp{i});
    icellAvg_adj{i}=nanmean(icell_adj{i});
    icellAvg_auto{i}=nanmean(icell_auto{i});
    icellAvg_ratio{i}=nanmean(icell_ratio{i});
    
end

%cd([filename])
%save([filename '/' basename '_647/' basename '_BTFSS'])

elseif recrunch==1
    load ([filename '/' basename '_647/' basename '_BTFSS'])
end

%Plot data
%Let's just measure intensity data first
figure, hold on, title('Intensity vs Time')
for i=1:ncells
    plot(time,icell_temp{1}(i,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')
saveas(gcf, [filename '/' basename '_647/' basename,'_647itemp.png'])

figure, title('Avgerage Intensity vs Time')
plot(time,icellAvg_temp{1},'-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')
saveas(gcf, [filename '/' basename '_647/' basename,'_647itempAvg.png'])

%Now let's measure adj intensity data and see if it's better
figure, hold on
for i=1:ncells
    plot(time,icell_adj{1}(i,:))
end
title('Intensity (adjusted for background)')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')
saveas(gcf, [filename '/' basename '_647/' basename,'_647iadj.png'])

figure, title('Avgerage Intensity (adjusted for background)')
plot(time,icellAvg_adj{1},'-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')
saveas(gcf, [filename '/' basename '_647/' basename,'_647iadjAvg.png'])

%Now let's measure adj intensity data - autoFluo and see if it's better
figure, hold on 
for i=1:ncells
    plot(time,icell_auto{1}(i,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
title('Intensity (adjusted for autofluorescence and background)')
fig2pretty
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')
saveas(gcf, [filename '/' basename '_647/' basename,'_647auto.png'])

figure, title('Avgerage Intensity (adjusted for autofluorescence and background)')
plot(time,icellAvg_auto{1},'-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')
saveas(gcf, [filename '/' basename '_647/' basename,'_647iautoAvg.png'])

%Finally, let's plot the ratio
figure, hold on
for i=1:ncells
    plot(time,icell_ratio{1}(i,:))
end
title('Intensity/Background Ratio')
xlabel('Time (s)')
ylabel('Intensity/Background (A.U.)')
fig2pretty
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')
saveas(gcf, [filename '/' basename '_647/' basename,'_647ratio.png'])

figure, title('Avgerage Intensity/Background Ratio')
plot(time,icellAvg_ratio{1},'-r')
xlabel('Time (s)')
ylabel('Intensity/Background (A.U.)')
fig2pretty
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')
saveas(gcf, [filename '/' basename '_647/' basename,'_647iratioAvg.png'])

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