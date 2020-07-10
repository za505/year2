%FluoColonyTrack.m
%
%Rico Rojas
%11/19/19
%
%Tracks colony size from epifluorescent image stacks in which the media is
%labeled with a tracer dye.  Save registered image stack to a folder by itself.
%The program will ask you to manually select an area that is background and
%then also select the area(s) in which cells are to determine where posts
%are.
%
%INPUT:
%
%dirname: full path of directory in which image stack is saved
%thresh: intensity cutoff below which cells are identified
%dr: Radius of disk with which to dilate posts

clear, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT 
basename='06162020_yhdL_nutrientShift_colony1';
folder=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/06162020_nutrientShiftAssay/' basename];
dirname=['/Users
    

/zarina/Downloads/NYU/Lab_2020_Summer/06162020_nutrientShiftAssay/' basename '/' basename '_2_a'];
thresh=0.9;%Parameter used to threshold image, to estimate this, use checkhist.
dr=6;%Radius to dilate posts by

tscale=60;%Frame rate in seconds
lscale=0.098;%Length scale in microns/pixel
smooth=20;%Number of frames over which to smooth growth rate profiles

vis=1;%Visualize tracking %vis=1 is tracking, vis=0 no tracking
checkhist=0;%Look at image histogram
checkposts=0;%Check microfluidic posts
recrunch=0;%Load previously cruched data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==1
    load([basename '_FCT'])
else
    
%Determine number of frames in movie
cd(dirname);
directory=dir('*.tif');
T=length(directory);

%Load first image
imagename=directory(1).name;
im=imread(imagename);

%Load last image
imagename=directory(T).name;
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

%Circle Initial Cell(s)
figure,imshow(im,[]), hold on, title('Select Initial Cell(s) and Press Enter')
k=0;
count=0;
while k~=1
count=count+1;
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
p1_2=min(point1,point2);%Calculate locations
p2_2=max(point1,point2);
offset = abs(point1-point2);%And dimensions
x = [p1_2(1) p1_2(1)+offset(1) p1_2(1)+offset(1) p1_2(1) p1_2(1)];
y = [p1_2(2) p1_2(2) p1_2(2)+offset(2) p1_2(2)+offset(2) p1_2(2)];
plot(x,y)
p1_2=round(p1_2);
p2_2=round(p2_2); 
rp1(count,:)=round(p1_2);
rp2(count,:)=round(p2_2);
end

    %Normalize images
    %im=norm16bit(im,ppix);
    
    %Determine background
    backim=im(p1(2):p2(2),p1(1):p2(1));
    bglevel=mean(mean(backim));

    %Normalize image
    [vals,bins]=imhist(im);
    [~,maxbin]=max(vals);
    nonzvals=find(vals);
    minind=bins(nonzvals(1));
  
    im=imadjust(im,[minind bglevel]/2^16,[0 1]);
    bw=im2bw(im,thresh);

    %Delete Posts
    bw=~bw;
    for n=1:count
        bw(rp1(n,2):rp2(n,2),rp1(n,1):rp2(n,1))=0;
    end
    posts=bw;
    posts=imdilate(posts,strel('disk',dr,4));
    for n=1:count
        posts(rp1(n,2):rp2(n,2),rp1(n,1):rp2(n,1))=0;
    end
    
    if checkposts==1
        figure,imshow(posts),pause
    end
    
for t=1:T
    t
    
    %Load images
    imagename=directory(t).name;
    im=imread(imagename);
  
    %Determine background
    backim=im(p1(2):p2(2),p1(1):p2(1));
    [counts,bins]=imhist(backim);
    [~,binnum]=max(counts);
    maxpos=bins(binnum);
    bglevel=mean(mean(backim));
    
    %Normalize image
    [vals,bins]=imhist(im);
    nonzvals=find(vals);
    minind=bins(nonzvals(1));


    
    %Delete Posts
    im2=im;
    im(posts)=bglevel;
     
    im=imadjust(im,[minind bglevel]/(2^16),[0 1]);
    im2=imadjust(im2,[minind bglevel]/(2^16),[0 1]);
    
    if checkhist==1
        figure
        imhist(im2)
        pause
    end

    %De-speckle image
    im=medfilt2(im);
   
    %Threshold image
    bw=im2bw(im,thresh);
    bw=~bw;
    
    %Clean image
    bw=bwmorph(bw,'clean');
    
    colonies_dil=imdilate(bw,strel('disk',1));
    col_edges=colonies_dil-bw;
    
    if vis==1 & t < 40
        figure
        imshow(imoverlay(im,col_edges,[1 0 0]))
        pause
        close all
    end
    
    %Count nonzero pixels
    bw_nz=nonzeros(bw);
    fluo_pix(t)=length(bw_nz);
    
end

%Dimensionalize variables
ColonyArea=fluo_pix*lscale^2;
Ncells=ColonyArea/20;
dA=ColonyArea(2:end)-ColonyArea(1:end-1);
eA=dA./ColonyArea(1:end-1)/(tscale);
eAsmooth=movingaverage(eA,smooth);

time=[0:T-1]*tscale;
tmid=(time(1:end-1)+time(2:end))/2;

end

%Plot data
figure
plot(time,Ncells/(Ncells(1)))
xlabel('Time (h)')
ylabel(strcat('Colony Size (# of cells) threshold = ', num2str(thresh)))
fig2pretty
cd(folder)
savefig([basename,'_colonySize.fig'])
%saveas(strcat(basename,'_colonySize.png'))

figure
plot(tmid/3600,eAsmooth*3600)
xlabel('Time (h)')
ylabel(strcat('Growth Rate (s^{-1}) threshold = ', num2str(thresh)))
fig2pretty
cd(folder)
savefig([basename,'_growthRate.fig'])
%saveas([basename,'_growthRateSize.png'])

cd(folder)
save([basename '_FCT'])


