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
basename='03172021_Exp3_colony1';
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03172021_analysis/' basename];
channels={[filename '/' basename '_647/' basename '_aligned']};
midSwitch=1; %not all the frames before a certain point have no dye
frameMid=17; %this is the FIRST frame that has no dye, only used if midSwitch=1
frameAuto=28; %this is the LAST frame that has no dye
frameBack=33; %this is the frame that you'll pick the background area from
recrunch=0;
auto=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==0

curdir=cd;
for i=1:length(channels)
    cd(channels{i}); 
    fluo_directory{i}=dir('*.tif');
end

load([filename '/' basename '_phase/' basename '_figures/' basename '_BTphase'])

%preallocate cells
icell_temp=cell(length(channels),1); %cellular intensities
icellAvg_temp=cell(length(channels),1); %population average cellular intensity
bgIntensity=zeros(1,T); %background intensity

if midSwitch==1
     autoFluo=zeros(height(pxls),frameAuto-frameMid); %intensity of cells prior to dye perfusion
else
    autoFluo=zeros(height(pxls),frameAuto);
end

for i=1:length(channels)
    
    cd(channels{i}); 
    intensity_temp=zeros(height(pxls)); %just a measure of cellular intensity
    
    %determine region where you'll measure background intensity
    imagename=fluo_directory{i}(frameBack).name;
    im=imread(imagename);
    
    [p1, p2]=getBackground(imagename);
    
    %measure autointensity
    imagename=fluo_directory{i}(frameAuto).name;
    im=imread(imagename);
    
    %now let's track intensity over time
    for t=1:T
        
        t
        
        %load the image
        imagename=fluo_directory{i}(t).name;
        im=imread(imagename);
       
        %measure background level
        bglevel = measureBackground(imagename, p1, p2);
        bgIntensity(t)=bglevel; 
        
        for j=1:height(pxls)
            
            if t <= frameAuto & midSwitch==0
                %calculate autofluorescence
                autoFluo(j, t)=mean(im(pxls{j,t}));
            elseif t <= frameAuto t >= frameMid
                %calculate autofluorescence
                autoFluo(j, t)=mean(im(pxls{j,t}));
            end
   
            %calculate intensity
            intensity_temp(j,t)=mean(im(pxls{j,t}));
            
%            figure
%            imshow(im)
%            hold on
%            
%            if isempty(boun{j,t})~=1
%            plot(boun{j,t}(:,1),boun{j,t}(:,2),'-w')
%            pause(0.5)
%            close
%            else
%                continue
%            end
        end
        
    end
    
    intensity_temp(intensity_temp==0)=NaN;
    autoFluo(autoFluo==0)=NaN;
    icell_temp{i}=intensity_temp;
    
end

for i=1:length(channels)
    icellAvg_temp{i}=nanmean(icell_temp{i});    
end

    
%take the population average of the autofluorescence 
autoFluo=mean(autoFluo);

%now take the temporal average and the background
%without dye
autoFluo=mean(autoFluo, 'all');
bgAdj = mean(bgIntensity(1, frameMid:frameAuto), 'all');

%now subtract the background without dye from the rest of the background
Iout = bgIntensity - bgAdj;

if auto==1
    
    %subtract the autofluorescence from the intensity
    Iin = icellAvg_temp{1} - autoFluo; 
    
else
    
    %just substract the background without dye
    Iin = icellAvg_temp{1} - bgAdj;
    
end 

%now, calculate the Iin/Iout ratio
ratio = Iin ./ Iout;
%remove values in the gap where there is no dye (weird things happen)
ratio(1,frameMid:frameAuto)=NaN

%cd([filename])
%save([filename '/' basename '_647/' basename '_BTFSS'])

elseif recrunch==1
    load ([filename '/' basename '_647/' basename '_figures/' basename '_BTfluo'])
end

%let's make a new directory to save plots and data
mkdir([filename '/' basename '_647/' basename '_figures'])
cd([filename '/' basename '_647/' basename '_figures'])
    
%Save the data
save([basename '_BTfluo'])

%Plot data
%Let's just measure intensity data first
figure, hold on, title('Intensity vs Time')
for i=1:height(pxls)
    plot(time,icell_temp{1}(i,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(0, '--k', 'LB + 647') %frame 1-16
xline(170, '--k', '*PBS + 5% detergent') %frame 17-28
xline(290, '--k', '*PBS + 647') %frame 29-41
xline(420, '--k', '*PBS + 647 + 20 mM NaCl') %frame 42+
saveas(gcf, [basename '_intensity.png'])

figure
plot(time,icellAvg_temp{1},'-r')
title('Avgerage Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(0, '--k', 'LB + 647') %frame 1-16
xline(170, '--k', '*PBS + 5% detergent') %frame 17-28
xline(290, '--k', '*PBS + 647') %frame 29-41
xline(420, '--k', '*PBS + 647 + 20 mM NaCl') %frame 42+
saveas(gcf, [basename,'_intensityAvg.png'])

%Now let's plot average background fluorescence
figure, hold on, title('Average Background Intensity vs Time')
plot(time,bgIntensity)
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(0, '--k', 'LB + 647') %frame 1-16
xline(170, '--k', '*PBS + 5% detergent') %frame 17-28
xline(290, '--k', '*PBS + 647') %frame 29-41
xline(420, '--k', '*PBS + 647 + 20 mM NaCl') %frame 42+
hold off
saveas(gcf, [basename,'_backgroundAvg.png'])

%Now let's plot adjusted background fluorescence
figure, hold on, title('Average Background Intensity vs Time')
plot(time,Iout)
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(0, '--k', 'LB + 647') %frame 1-16
xline(170, '--k', '*PBS + 5% detergent') %frame 17-28
xline(290, '--k', '*PBS + 647') %frame 29-41
xline(420, '--k', '*PBS + 647 + 20 mM NaCl') %frame 42+
hold off
saveas(gcf, [basename,'_backgroundAdj.png'])


%Now let's plot adjusted fluorescence
figure, hold on, title('Intensity vs Time')
plot(time,Iin)
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(0, '--k', 'LB + 647') %frame 1-16
xline(170, '--k', '*PBS + 5% detergent') %frame 17-28
xline(290, '--k', '*PBS + 647') %frame 29-41
xline(420, '--k', '*PBS + 647 + 20 mM NaCl') %frame 42+
hold off
saveas(gcf, [basename,'_IntensityAdj.png'])

%Plot the Iin/Iout ratio over time
figure, hold on
plot(time(1,:),ratio(1,:))
xlabel('Time (s)')
ylabel('Intensity/Background')
fig2pretty
xline(0, '--k', 'LB + 647') %frame 1-16
xline(170, '--k', '*PBS + 5% detergent') %frame 17-28
xline(290, '--k', '*PBS + 647') %frame 29-41
xline(420, '--k', '*PBS + 647 + 20 mM NaCl') %frame 42+
hold off
saveas(gcf, [basename,'_ratio.png'])


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