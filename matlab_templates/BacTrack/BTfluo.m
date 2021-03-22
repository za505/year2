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
basename='03172021_Exp2_colony2';
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03172021_analysis/' basename];
channels={[filename '/' basename '_FSS/' basename '_full']};
%channels={[filename '/' basename '_FSS/' basename '_full']};
switch1=40; %frame during switch 1
switch2=50;
switch3=60;
recrunch=0;
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
            
            bglevel=getBackground(imagename);
            
        elseif t==switch1
            
            bglevel=getBackground(imagename);
        
        elseif t==switch2
            
            bglevel=getBackground(imagename);
            
        elseif t==switch3
            
            bglevel=getBackground(imagename);
        
        end 
        
        for j=1:ncells
            intensities_temp(j,t)=mean(im(pixels{j,t}));
            intensity_ratio(j,t)= bglevel/intensities_temp(j,t);
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
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')
%saveas(gcf, [filename '/' basename '_FSS/' basename,'_FSStrace.png'])

figure
%plot(time(tstart:end),icell_av{1}(tstart:end),'-r')
plot(time,icell_av{1},'-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')
%saveas(gcf, [filename '/' basename '_FSS/' basename,'_FSSavg.png'])

figure, hold on, 
for i=1:ncells
    plot(time, intensity_ratio(i,:))
end
xlabel('Time (s)')
ylabel('Intensity Ratio (background/cell insensity)')
fig2pretty
xline(90, '--k', '*PBS + 5% detergent')
xline(210, '--k', '*PBS + 647')
xline(330, '--k', '*PBS + 647 + FSS')
xline(450, '--k', '*PBS + 647 + CF')
xline(570, '--k', '*PBS + 647 + AF')

%save([filename '/' basename '_FSS/' basename '_BTFSS'])

function bglevel = getBackground(imagename)
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

        %Determine background
        backim=im2(p1(2):p2(2),p1(1):p2(1));
        [counts,bins]=imhist(backim);
        [~,binnum]=max(counts);
        maxpos=bins(binnum);
        bglevel=mean(mean(backim));
end 