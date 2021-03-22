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
basename='03172021_Exp1_colony1';
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03172021_analysis/' basename];
channels={[filename '/' basename '_FSS/' basename '_full']};
%channels={[filename '/' basename '_FSS/' basename '_full']};
switch1=20; %frame during switch 1
switch2=31;
switch3=44;
switch4=57;
recrunch=0;
removeData=0; %change to 0 if you don't want to remove data points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==0;

curdir=cd;
for i=1:length(channels)
    cd(channels{i}); 
    fluo_directory{i}=dir('*.tif');
end

load([filename '/' basename '_phase/' basename '_figures/' basename '_BTphase'])
%load([basename '_BTlab'])

intensities=cell(length(channels),1);


for i=1:length(channels)
    cd(channels{i}); 
    intensities_temp=zeros(size(lcell));
    intensity_ratio=zeros(size(lcell));
    
    for t=1:T
        t
        imagename=fluo_directory{i}(t).name;
        im=imread(imagename);
        
        %let's normalize the image to see if it helps
        %ppix=0.5;
        %im=norm16bit(im,ppix);
        
        %let's check some of the images just to be sure
%         if t==10 | t==30 | t==50 | t==60
%             imtool(im)
%             pause
%         end

        if t==1
            
            [p1, p2]=getBackground(imagename);
            
        elseif t==switch1
            
            [p1, p2]=getBackground(imagename);
        
        elseif t==switch2
            
            [p1, p2]=getBackground(imagename);
            
        elseif t==switch3
            
            [p1, p2]=getBackground(imagename);
            
        elseif t==switch4
            
            [p1, p2]=getBackground(imagename);
            
        end 
        
        %measure background level
        bglevel = measureBackground(imagename, p1, p2);
        
        for j=1:ncells
            
            %remove weird data points
            if removeData~=0 & j == removeData
                intensities_temp(j,t)=NaN;
            else
                %calculate intensity and intensity ratio
                intensities_temp(j,t)=mean(im(pixels{j,t}));
                intensity_ratio(j,t)= intensities_temp(j,t)/bglevel;
            end
            
        end
    end
    
    intensities_temp(intensities_temp==0)=NaN;
    icell{i}=intensities_temp;
end

icell_av=cell(length(channels),1);
for i=1:length(channels)
    icell_av{i}=nanmean(icell{i});
end
    
cd([filename])
%save([filename '/' basename '_FSS/' basename '_BTFSS'])

elseif recrunch==1
    load ([filename '/' basename '_FSS/' basename '_BTFSS'])
end

%Plot data
figure, hold on, 
for i=1:ncells
    plot(time,icell{1}(i,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(60, '--k', '*PBS + 5% detergent')
xline(114, '--k', '*PBS + FSS + FSS') %frame 19-30
xline(234, '--k', '*PBS + FSS + FSS + 6.66 mM Mg2+') %frame 31-43
xline(354, '--k', '*PBS + FSS + FSS + 12.33 mM Mg2+') %frame 44-56
xline(474, '--k', '*PBS + FSS + FSS + 20 mM Mg2+') %frame 57-69
saveas(gcf, [filename '/' basename '_FSS/' basename,'_FSStrace.png'])

figure
%plot(time(tstart:end),icell_av{1}(tstart:end),'-r')
plot(time,icell_av{1},'-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(60, '--k', '*PBS + 5% detergent')
xline(114, '--k', '*PBS + FSS + FSS') %frame 19-30
xline(234, '--k', '*PBS + FSS + FSS + 6.66 mM Mg2+') %frame 31-43
xline(354, '--k', '*PBS + FSS + FSS + 12.33 mM Mg2+') %frame 44-56
xline(474, '--k', '*PBS + FSS + FSS + 20 mM Mg2+') %frame 57-69
saveas(gcf, [filename '/' basename '_FSS/' basename,'_FSSavg.png'])

figure, hold on, 
for i=1:ncells
    plot(time, intensity_ratio(i,:))
end
xlabel('Time (s)')
ylabel('Intensity Ratio (cell intensity/background)')
fig2pretty
xline(60, '--k', '*PBS + 5% detergent')
xline(114, '--k', '*PBS + FSS + FSS') %frame 19-30
xline(234, '--k', '*PBS + FSS + FSS + 6.66 mM Mg2+') %frame 31-43
xline(354, '--k', '*PBS + FSS + FSS + 12.33 mM Mg2+') %frame 44-56
xline(474, '--k', '*PBS + FSS + FSS + 20 mM Mg2+') %frame 57-69
saveas(gcf, [filename '/' basename '_FSS/' basename,'_FSSratio.png'])

save([filename '/' basename '_FSS/' basename '_BTFSS'])

function [p1, p2]=getBackground(imagename)
        
        %Load last image
        %imagename=fluo_directory{i}(t).name;
        im2=imread(imagename);

        %Determine Background
        figure,imshow(im2,[]), hold on, title('Select Background')
        k=waitforbuttonpress;
        set(gcf,'Pointer','fullcross')
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