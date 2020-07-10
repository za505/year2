%BacTrack.m
%Tracks bacterial growth from phase image stacks.  
%Customized for B. subtilis.
%Bactrack2.m is customized for B. subtilis filaments.


%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save phase image stack in a directory
%by itself.  Also save the micromanager metadata file as 'basename.txt' in
%the matlab path.
%
%INPUT:
%basename: name of the image stack.
%dirname:the full pathname of the directory where you saved the image
%        stack.
%metaname(optional):full or relative pathname of micromanager metadata file from
%which to extract time points.  If it is relative path name, the
%directory in which it is saved must be on the matlab path.
%lscale: microscope calibration in microns per pixels.
%sm: width of the Gaussian filter used in edge finder equals sm*sqrt(2).
%minL: minimum length of cells;
%minW: minimum width of cells;
%maxW: maximum width of cells;
%recrunch:0 or 1.  if you've already tracked the data set and just want to
%         re-plot the data enter 1.
%
%OUTPUT:
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
%%User Input
basename='06162020_yhdL_nutrientShift_colony1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/06162020_nutrientShiftAssay/' basename '/' basename '_2_b'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/06162020_nutrientShiftAssay/' basename];%Directory to save the output .mat file to.
%metaname=['/Users/Rico/Documents/MATLAB/Matlab Ready/' basename '/metadata.txt'];%Name of metadata file.  Will only work if images were taken with micromanager.
lscale=0.08;%%Microns per pixel.
tscale=60;%Frame rate.
thresh=0;%For default, enter zero.
IntThresh=20000;%Threshold used to enhance contrast. Default:35000
dr=1;%Radius of dilation before watershed 
sm=2;%Parameter used in edge detection
minL=2;%Minimum cell length
minW=0.2;%Minimum cell width
maxW=1.5;%Maximum cell width
minA=50;%Minimum cell area. default 50
cellLink=4;%Number of frames to ignore missing cells when tracking frame to frame
recrunch=0;%Display data from previously crunched data? 0=No, 1=Yes.
vis=1;%Display cell tracking? 0=No, 1=Yes.
checkhist=0;%Display image histogram? 0=No, 1=Yes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==1
    load([basename '_BT'])
else

%Determine number of frames
curdir=cd;
cd(dirname);
directory=dir('*.tif');
T=length(directory);

cd(curdir);
path(dirname,path)

nc=zeros(1,T);
allcentroids=[];
cellnum=[];
tstamp=[];

%Pre-allocate matrices
ewav=zeros(1,T);
ewstd=zeros(1,T);
ewste=zeros(1,T);
ewndp=zeros(1,T);
wav=zeros(1,T);
wstd=zeros(1,T);
wste=zeros(1,T);
vav=zeros(1,T-1);
vstd=zeros(1,T-1);
vste=zeros(1,T-1);
ndp=zeros(1,T-1);
avav=zeros(1,T-1);
avstd=zeros(1,T-1);
avste=zeros(1,T-1);
a=zeros(1,T);
w=zeros(1,T);
l=zeros(1,T);
DS=zeros(1,T);
boun=cell(1,T);
pole=zeros(1,T);
mline=cell(1,T);

%Load first image
imagename=directory(1).name;
im=imread(imagename);
[imM,imN]=size(im);
labels=zeros(imM,imN,T);
labels2=zeros(imM,imN,T);

%Track cells
for t=1:T
    t
    
    %Load image
    imagename=directory(t).name;

    im=imread(imagename);
    [imM,imN]=size(im);
    
    %De-speckle image
    im=medfilt2(im);
    
    %Normalize images
    ppix=0.5;
    im=norm16bit(im,ppix);
    
    %Enhance contrast
    imc=imcomplement(im);
    if checkhist==1;
        figure,imhist(imc),pause;
    end
    
    if thresh==0;
        [imcounts,bins]=imhist(imc);
        [imcounts,idx]=sort(imcounts);
        bins=bins(idx);
        thresh1=bins(end-1);
    else
        thresh1=thresh;
    end
    imc=imadjust(imc,[thresh1/65535 1],[]);   
     
    %Find edges
    [ed2,thresh2]=edge(imc,'canny',[],sm*sqrt(2));
    
    %Clean image
    cc=bwconncomp(ed2,8);
    stats=regionprops(cc,imc,'Area','MeanIntensity');
    idx=find([stats.Area]>minA&[stats.Area]<1e5&[stats.MeanIntensity]>IntThresh);
    ed2=ismember(labelmatrix(cc),idx);
    
    %Close gaps in edges
    despurred=bwmorph(ed2,'spur');
    spurs=ed2-despurred;
    [spy,spx]=find(spurs);
    for k=1:length(spx)
        ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)=ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)+rot90(ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1),2);
        ed2(spy(k),spx(k))=1;
    end
    ed2=bwmorph(ed2,'bridge'); 
  
    se=strel('disk',dr);
    ed2=imdilate(ed2,se);
    ed2=bwmorph(ed2,'thin',2);
    
    %Identify cells based on size and intensity
    ed3=~ed2;
    ed3(1,:)=ones(1,imN);
    ed3(end,:)=ones(1,imN);
    ed3(:,1)=ones(imM,1);
    ed3(:,end)=ones(imM,1);
    
    cc=bwconncomp(ed3,4);
    stats=regionprops(cc,imc,'Area','MeanIntensity');
    idx=find([stats.Area]>minA&[stats.Area]<1e5&[stats.MeanIntensity]>3e4);
    ed4=ismember(labelmatrix(cc),idx);
    
    %Find cell areas and centroids
    bw=bwmorph(ed4,'thicken');
    [P,bw]=bwboundaries(bw,4,'noholes');
    stats=regionprops(bw,'Area','Centroid','PixelIdxList');
    
    L=bwlabel(bw);    
    labels(:,:,t)=L;
    labels2(:,:,t)=bw;
    
    nc(t)=length(P);
    areas=[stats.Area];
    cents=cat(1,stats.Centroid);
    a(nc(t),t)=0;
    a(1:nc(t),t)=[stats.Area]';
    centroids=cents;
     
    %Calculate smooth cell contours
    for n=1:nc(t)
         
        rP=[P{n}(:,2),P{n}(:,1)];
        px=[rP(1:end-1,1);rP(1:end-1,1);rP(:,1)];
        py=[rP(1:end-1,2);rP(1:end-1,2);rP(:,2)];
        sp=length(rP);
        dS=sqrt(diff(px).^2+diff(py).^2);
        S=[0 cumsum(dS)'];

        px=csaps(S,px,0.05,S);
        py=csaps(S,py,0.05,S);
        
        px=px(sp+1:2*sp);
        py=py(sp+1:2*sp);
        
        px(end)=px(1);
        py(end)=py(1);
        
        dS=sqrt(diff(px).^2+diff(py).^2);
        S=[0 cumsum(dS)];
        ls=length(S);
        DS(n,t)=S(end)/(ls-1);
        Sn=(0:DS(n,t):S(end));
        nx=spline(S,px,Sn);
        ny=spline(S,py,Sn);
        
        boun{n,t}=[nx',ny'];
        pxls{n,t}=stats(n).PixelIdxList;
        
    end
    allcentroids=[allcentroids;centroids];
    tstamp=[tstamp;ones(nc(t),1)*t];
    cellnum=[cellnum;(1:nc(t))'];
    
if vis==1
   figure
   imshow(im)
   hold on
   for k=1:nc(t)
       plot(boun{k,t}(:,1),boun{k,t}(:,2),'-r')
   end
    
  pause
  close all
end
    toc

end

%Calculate cell length, width, etc.
for t=1:T
    t
    for n=1:nc(t)    
         X=boun{n,t}(:,1);
         Y=boun{n,t}(:,2);   
                       
         [sX,~]=size(X);
         
         %Find poles
         [X,Y,pole(n,t)]=polefinder(X,Y);
         
         %Create mesh
         npts=min(pole(n,t),sX-pole(n,t)+1);
         S=(0:DS(n,t):(sX-1)*DS(n,t));
         
         s1=(0:S(pole(n,t))/(npts-1):S(pole(n,t)));
         s2=(S(pole(n,t)):(S(end)-S(pole(n,t)))/(npts-1):S(end));
         xc1=spline(S(1:pole(n,t)),X(1:pole(n,t)),s1);
         xc2=spline(S(pole(n,t):end),X(pole(n,t):end),s2);
         yc1=spline(S(1:pole(n,t)),Y(1:pole(n,t)),s1);
         yc2=spline(S(pole(n,t):end),Y(pole(n,t):end),s2);
         xc2=fliplr(xc2);
         yc2=fliplr(yc2);
         
         %Calculate midline
         mline{n,t}=[(xc1+xc2)'/2,(yc1+yc2)'/2];
         dxy=diff(mline{n,t}).^2;
         dl=sqrt(dxy(:,1)+dxy(:,2));
         l(n,t)=sum(dl);
         
         %Calculate width
         ls=[0 cumsum(dl)'];
         [~,mpos1]=min(abs(ls/l(n,t)-0.25));
         [~,mpos2]=min(abs(ls/l(n,t)-0.75));
         
         widths=sqrt((xc1-xc2).^2+(yc1-yc2).^2);
         w(n,t)=(widths(mpos1)+widths(mpos2))/2;
          
    end
end


%Extract timepoints from metadata if it exists
if exist('metaname')==1
    if exist(metaname)==2
        %Extract timepoints from metadata
        tpoints=metadata(metaname);
        
        %Fix bug where micromanager screws up its timing
        dtime=diff(tpoints(1,:));
        fdt=find(dtime>2*(dtime(1)));
        if isempty(fdt)~=1
            fdt=fdt(1);
            tpoints(:,fdt+1:end)=tpoints(:,fdt+1:end)-tpoints(1,fdt+1)+tpoints(1,fdt)+(tpoints(1,fdt)-tpoints(1,fdt-1));
        end
    else
        tpoints=[0:T-1]*tscale;
    end
else
    tpoints=[0:T-1]*tscale;
end

time=tpoints(1,:);
time2=tpoints(end,:);

%Fix bug where micromanager screws up its timing
dtime=diff(time);
fdt=find(dtime>2*(dtime(1)));
if isempty(fdt)~=1
    fdt=fdt(1);
    time(fdt+1:end)=time(fdt+1:end)-dtime(fdt)+dtime(fdt-1);
    time2(fdt+1:end)=time2(fdt+1:end)-dtime(fdt)+dtime(fdt-1);
end

%Track cells frame to frame
tracks=zeros(size(im));
rcents=round(allcentroids);
linind=sub2ind(size(im),rcents(:,2),rcents(:,1));
tracks(linind)=1;

nhood=[0,1,0;1,1,1;0,1,0];
tracks=imdilate(tracks,strel('disk',cellLink));
overlay1=imoverlay(im,tracks,[.3 1 .3]);

[tracksL,ncells]=bwlabel(tracks,8);

lcell=zeros(ncells,T);
wcell=zeros(ncells,T);
acell=zeros(ncells,T);
pcell=zeros(ncells,T);
B=cell(ncells,T);
pixels=cell(ncells,T);
mlines=cell(ncells,T);
lcents=length(allcentroids);

for i=1:lcents
    cellid=tracksL(linind(i));
    lcell(cellid,tstamp(i))=l(cellnum(i),tstamp(i));
    wcell(cellid,tstamp(i))=w(cellnum(i),tstamp(i));
    acell(cellid,tstamp(i))=a(cellnum(i),tstamp(i));
    B{cellid,tstamp(i)}=boun{cellnum(i),tstamp(i)};
    pixels{cellid,tstamp(i)}=pxls{cellnum(i),tstamp(i)};
    mlines{cellid,tstamp(i)}=mline{cellnum(i),tstamp(i)};
    pcell(cellid,tstamp(i))=pole(cellnum(i),tstamp(i));
end

%Throw away cells with only one or two time points
delind=[];
for i=1:ncells
    if length(nonzeros(lcell(i,:)))<=2
        delind=[delind;i];
    end
end

lcell(delind,:)=[];
wcell(delind,:)=[];
acell(delind,:)=[];
pcell(delind,:)=[];
B(delind,:)=[];
pixels(delind,:)=[];
mlines(delind,:)=[];
[ncells,~]=size(lcell);

lcell(lcell==0)=NaN;
wcell(wcell==0)=NaN;
acell(acell==0)=NaN;
pcell(pcell==0)=NaN;

%Dimsionalize the variables
lcell=lcell*lscale;
wcell=wcell*lscale;
acell=acell*lscale^2;

%Throw away cells that are too short or too fat or too skinny
lcell(lcell<minL|wcell>maxW|wcell<minW)=NaN;
wcell(lcell<minL|wcell>maxW|wcell<minW)=NaN;
acell(lcell<minL|wcell>maxW|wcell<minW)=NaN;

%Calculate circumferential strain
wcell(isnan(wcell))=0;
ew=zeros(size(lcell));
for i=1:ncells
    ew(i,:)=wcell(i,:)/mean(wcell(i,:));
end

%Calculate the growth rate
deltat=time(2:end)-time(1:end-1);
v=(lcell(:,2:end)-lcell(:,1:end-1))./((lcell(:,1:end-1)+lcell(:,2:end))/2);
av=(acell(:,2:end)-acell(:,1:end-1))./((acell(:,1:end-1)+acell(:,2:end))/2);
for i=1:ncells
    v(i,:)=v(i,:)./deltat;
    av(i,:)=av(i,:)./deltat;
end

%Calculate total length of all cells
lcell(isnan(lcell))=0;
acell(isnan(acell))=0;
for t=1:T
    ltotal(t)=sum(nonzeros(l(:,t)));
    atotal(t)=sum(nonzeros(a(:,t)));
end
lcell(lcell==0)=NaN;
acell(acell==0)=NaN;

%Throw out outliers and calculate the average width,strain and strain rate 
%across cells
v(isnan(v))=0;
av(isnan(av))=0;
ew(isnan(ew))=0;

for t=1:T-1
    vav(t)=mean(nonzeros(v(:,t)));
    vstd(t)=std(nonzeros(v(:,t)));
end
vavm=ones(ncells,1)*vav;
vstdm=ones(ncells,1)*vstd;

inddel=abs(v-vavm)>2*vstdm&vstdm~=0;

v(inddel)=0;
av(inddel)=0;
lcell(inddel)=0;
acell(inddel)=0;
wcell(inddel)=0;
ew(inddel)=0;

for t=1:T
    wav(t)=mean(nonzeros(wcell(:,t)));
    wstd(t)=std(nonzeros(wcell(:,t)));
    wste(t)=wstd(t)./length(nonzeros(wcell(:,t)));
    ewav(t)=mean(nonzeros(ew(:,t)));
    ewstd(t)=std(nonzeros(ew(:,t)));
    ewste(t)=ewstd(t)./length(nonzeros(ew(:,t)));
end
for t=1:T-1
    vav(t)=mean(nonzeros(v(:,t)));
    vstd(t)=std(nonzeros(v(:,t)));
    ndp(t)=length(nonzeros(v(:,t)));
    vste(t)=vstd(t)/sqrt(ndp(t));
    avav(t)=mean(nonzeros(av(:,t)));
    avstd(t)=std(nonzeros(av(:,t)));
    avste(t)=avstd(t)/ndp(t);
end

v(v==0)=NaN;
av(av==0)=NaN;
lcell(lcell==0)=NaN;
wcell(wcell==0)=NaN;
acell(acell==0)=NaN;
ew(ew==0)=NaN;

end

%Plot data
figure(1), title('Cell Length vs. Time')
clf
hold on
for i=1:ncells  
    lcell(i,:)=movingaverage2(lcell(i,:),3);
    %indx=isnan(lcell(i,:))~=1;
    %indx=find(indx);
    %plot(time(indx),lcell(i,indx))
    plot(time(1:end),lcell(i,1:end)) 
end
xlabel('Time (s)')
ylabel('Length (\mum)')
fig2pretty

figure(2), title('Cell Width vs. Time')
hold on
for i=1:ncells
    plot(time,wcell(i,:)) 
end
plot(time,wav,'-r','LineWidth',2)
xlabel('Time (s)')
ylabel('Width (/mum)')
fig2pretty

% figure(4), title('Circumferential Strain vs. Time')
% hold on
% for i=1:ncells
%     plot(time,ew(i,:)) 
% end
% plot(time,ewav,'-r','LineWidth',2)
% xlabel('t (s)')
% ylabel('\epsilon_w')
% fig2pretty

tmid=(time(2:end)+time(1:end-1))/2;

figure(5), title('Elongation Rate vs. Time')
hold on
for i=1:ncells
    plot(tmid,v(i,:))
end
plot(tmid,vav,'-r')
xlabel('Time (s)')
ylabel('Elongation Rate (s^{-1})')
fig2pretty

figure(6), title('Elongation Rate vs. Time')
hold on
ciplot((vav-vstd)*3600,(vav+vstd)*3600,tmid,[0.75 0.75 1])
plot(tmid,vav*3600,'-r')
xlabel('Time (s)')
ylabel('Elongation (hr^{-1})')
fig2pretty

cd(dirname);
save([basename '_BT'])
save([basename '_BTlab'],'labels','labels2','-v7.3')
