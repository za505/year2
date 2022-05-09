%LinTrack2.m
%% Part 1: is BacTrack2.m
%Tracks bacterial growth from phase image stacks.  Customized for E. coli.

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
%ndp: vector of lenth T-1 with number of data points averaged over.

%Calls on the following m-files:
%norm16bit.m
%polefinder.m
%cellcurvature.m
%metadata.m
%extrema.m
%EffectiveLength.m
%fig2pretty.m
%movingaverage


clear
close all
cd /Users/dylanfitzmaurice/Documents/MATLAB/data

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='FC2_120920_P1';%Name of the image stack, used to save file.
%dirname=['/Users/dylanfitzmaurice/Documents/MATLAB/Matlab Ready/' basename '/' basename '_1_a'];%Directory that the image stack is saved in
dirname=['/Users/dylanfitzmaurice/Documents/MATLAB/Matlab Ready/' basename '/' basename '_1_a'];%_cell1
%metaname=['/Volumes/Lacie/Matlab Ready/' basename '/metadata.txt'];%Name of metadata file.  Will only work if images were taken with micromanager.
recrunch=0;%Display data from previously crunched data? 0=No, 1=Yes.
checkhist=0;%Display image histogram? 0=No, 1=Yes.
vis=1;%Display cell tracking? 0=No, 1=Yes.
lscale=0.08;%Microns per pixel.
tscale=10;%Frame rate.
minL=0.5;%Minimum cell length
minW=0.6;%Minimum cell width, default 0.6
maxW=2.5;%Maximum cell width, default 2.5
minA=50;%Minimum cell area. default 50
maxA=0.009e5;%Maximum maximum area. default 1e5, 0.009e5 to not count holes
minI=2e4;%Minimum cell intensity. default 2e4, for shock 
thresh=2000;%Threshold used to enhance contrast. Default:35000
thresh2=0.0009;%Threshold used in edge finder. Default:0.01, 0.001 seems to 
% work better for me
wthresh=1;%Threshold used to find cell septa, default=1.
ppix=0.5;% default = 0.5 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==1
    load([basename '_BT'])
else
    
    %Determine number of frames
    curdir=cd;
    cd(dirname);
    directory=dir('*.tif');
    
    T=length(directory);
    %T=26;
    
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
    wvav=zeros(1,T-1);
    avav=zeros(1,T-1);
    avstd=zeros(1,T-1);
    avste=zeros(1,T-1);
    a=zeros(1,T);
    o=zeros(1,T);
    w=zeros(1,T);
    l=zeros(1,T);
    DS=zeros(1,T);
    boun=cell(1,T);
    pole=zeros(1,T);
    mline=cell(1,T);
    total_vol=cell(1,T);
    
    %Load first image
    imagename=directory(1).name;
    im=imread(imagename);
    [imM,imN]=size(im);
    labels=zeros(imM,imN,T);
    labels2=zeros(imM,imN,T);
    
    im=imread(imagename);
    level=graythresh(im);
    bw=im2bw(im,level);
    bw=~bw;
    stats_lab=cell(1,T);
    
    %Find cells in each image
    for t=flip(1:T)
        %Load image
        imagename=directory(t).name;
        im=imread(imagename);
        [imM,imN]=size(im);
        %figure,imhist(im),pause
        
        %De-speckle
        im=medfilt2(im,[3 3]); %change back to 3
        
        %Normalize images
        %ppix=0.5;
        im2=norm16bit(im,ppix);
        %figure, imshow(im2), pause
        
        %Enhance contrast
        imc=imcomplement(im2);
        if checkhist==1;
            figure,imhist(imc),pause;
        end
        imc=imadjust(imc,[thresh/65535 1],[]);
        %figure, imshow(imc), pause
        
        %Find edges
        [ed2,thresh2]=edge(imc,'zerocross',thresh2);
        %figure,imshow(ed2),pause
        
        %Clean image
        cc=bwconncomp(ed2,8);
        stats=regionprops(cc,imc,'Area','MeanIntensity');
        idx=find([stats.Area]>minA&[stats.Area]<maxA&[stats.MeanIntensity]>minI);
        ed2=ismember(labelmatrix(cc),idx);
        %figure,imshow(ed2),pause
        
        %Close gaps
        despurred=bwmorph(ed2,'spur');
        spurs=ed2-despurred;
        [spy,spx]=find(spurs);
        for k=1:length(spx)
            ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)=ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)+rot90(ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1),2);
            ed2(spy(k),spx(k))=1;
        end
        ed2=bwmorph(ed2,'bridge');
        
        %Identify cells based on size and intensity
        ed3=~ed2;
        ed3(1,:)=ones(1,imN);
        ed3(end,:)=ones(1,imN);
        ed3(:,1)=ones(imM,1);
        ed3(:,end)=ones(imM,1);
        
        cc=bwconncomp(ed3,4);
        stats=regionprops(cc,imc,'Area','MeanIntensity');
        %figure,hist([stats.MeanIntensity])
        idx=find([stats.Area]>minA&[stats.Area]<maxA&[stats.MeanIntensity]>minI);
        ed4=ismember(labelmatrix(cc),idx);
        
        ed5=ed4;
        ed5(ed2==1)=1;
        ed6=imdilate(ed5,[0 1 0;1 1 1;0 1 0]);
        ed7=ed6-ed5;
        
        %Find septa
        ed8=bwdist(~ed4);
        ed9=zeros(size(ed8));
        ed9(ed7==1)=100;
        
        H=fspecial('gaussian',4,1);
        ed9=imfilter(ed9,H,'replicate');
        ed9(ed8>wthresh)=-100;
        ed9(ed6==0)=-100;
        %figure,imshow(ed9),pause
        
        %Final watershed
        L=watershed(ed9);
        ed10=L==0;
        ed11=~ed10;
        cc1=bwconncomp(ed11,4);
        stats=regionprops(cc1,imc,'Area','MeanIntensity','Centroid');
        idx=find([stats.Area]>minA&[stats.Area]<maxA&[stats.MeanIntensity]>minI);
        ed12=ismember(labelmatrix(cc1),idx);
        L=bwlabel(ed12);
        stats=regionprops(ed12,'Area','Centroid','Orientation','PixelIdxList');
        
        labels(:,:,t)=L; %important for lineage tracking
        labels2(:,:,t)=bw;
       
        nc(t)=length([stats.Area]);
        areas=[stats.Area];
        cents=cat(1,stats.Centroid);
        %cents_indx{:,t}=cents;
        a(nc(t),t)=0;
        o(nc(t),t)=0;
        a(1:nc(t),t)=[stats.Area]';
        o(1:nc(t),t)=[stats.Orientation]';
        centroids=cents;
        
        %Calculate smooth cell contours
        for n=1:nc(t)
            
            ed13=L==n;
            ed13=imdilate(ed13,[0 1 0;1 1 1;0 1 0]);
            
            P(n)=bwboundaries(ed13,4,'noholes');
            
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
            %figure
            %imshowpair(L_tracking{1,t})
            imshow(im)
            hold on
            for k=1:nc(t)
                plot(boun{k,t}(:,1),boun{k,t}(:,2),'-r')
                hold on
                plot(cents(:,1),cents(:,2),'b*')
                %imshow(L_tracking{1,t})
            end
            
            %pause
            pause(0.1)
        %pause
        clf
            %close all
        end
        
        toc
        
    end
       
    %Calculate cell length, width, etc.
    for t=flip(1:T)
        t;
        for n=1:nc(t)
            X=boun{n,t}(:,1);
            Y=boun{n,t}(:,2);
            
            [sX,~]=size(X);
            
            %Find poles
            [X,Y,pole(n,t)]=polefinder(X,Y);
            
            %Create mesh along cell outline
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
            
%             %calculate volume Dylan Edit
%             %mean_widths = mean([widths(1:end-1);widths(2:end)]);
%             %seg_area = widths * mline;
%             %seg_volume = pi * seg_area .* (mean_widths/2);
            xc1 = xc1';
            yc1 = yc1';
           
            xc2 = xc2';
            yc2 = yc2';
            
            
% r_2 is the radius from midline to cell surface using xc2 and yc2, hence r_2
            xr_2 = (xc2(1:end,1)-(mline{n,t}(1:end,1))).^2;
            yr_2 = (yc2(1:end,1)-(mline{n,t}(1:end,2))).^2;
            
            r_2 = sqrt(xr_2 + yr_2 ); % this is the distance from midline to surface
            av_r_2 = 0.5 * (r_2(1:end-1,1) + r_2(2:end,1));
           
%             av_xr_2 = 0.5 * (xr_2(1:end-1,1) + xr_2(2:end,1));
%             av_yr_2 = 0.5 * (yr_2(1:end-1,1) + yr_2(2:end,1));
%             av_r_2 = [av_xr_2 av_yr_2];
%             r_2 = abs(sqrt((av_r_2(1:end,1).^2)-(av_r_2(1:end,2).^2)));
%             r_2 = r_2 .* lscale;
%             area_2 = r_2 .* dl .* lscale;
            
            area_2 = av_r_2 .* dl;

            vol_2 = area_2 .* pi .* av_r_2;
                    

            % r_1
            xr_1 = (xc1(1:end,1)-(mline{n,t}(1:end,1))).^2;
            yr_1 = (yc1(1:end,1)-(mline{n,t}(1:end,2))).^2;
            
            r_1 = sqrt(xr_1 + yr_1 ); % this is the distance from midline to surface
            av_r_1 = 0.5 * (r_1(1:end-1,1) + r_1(2:end,1));
            
            area_1 = av_r_1 .* dl;

            vol_1 = area_1 .* pi .* av_r_1;

            total_v(n,t) = sum(vol_2)+sum(vol_1);
            
            total_vol = total_v;
%             
% %             xr_1 = abs((sqrt(xc1(1:end,1).^2)-(mline{n,t}(1:end,1).^2)));
% %             yr_1 = abs((sqrt(yc1(1:end,1).^2)-(mline{n,t}(2:end,1).^2)));
% %             r_1 = xr_1 + yr_1;
%             
%             %av_r_2 = mean([r_2(1:end-1);r_2(2:end)]);
%             %seg_area2 = av_r_2.*mline{n,t}(1:end,1);
        end
    end
end

%total_vol = num2cell(total_vol);

%Extract timepoints from metadata
if exist('metaname')==1
    if exist(metaname)==2
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

%[r,c] = find(bwlabel(L_tracking{1,27})>=1)

%Track cells frame to frame
tracks=zeros(size(im));
rcents=round(allcentroids);
linind=sub2ind(size(im),rcents(:,2),rcents(:,1));
%linind=sub2ind(size(im),r,c);%new linind for tracking binary image
tracks(linind)=1;

nhood=[0,1,0;1,1,1;0,1,0];
tracks=imdilate(tracks,strel('disk',2));
%imshow(tracks)

[tracksL,ncells]=bwlabel(tracks);
%imshow(tracksL)

lcell=zeros(ncells,T);
wcell=zeros(ncells,T);
acell=zeros(ncells,T);
pcell=zeros(ncells,T);
B=cell(ncells,T);
pixels=cell(ncells,T);
mlines=cell(ncells,T);
total_vols = zeros(ncells,T);%%% new
%total_vols = cell(ncells,T);%%%
lcents=length(allcentroids);

for i=1:lcents
    cellid=tracksL(linind(i));
    lcell(cellid,tstamp(i))=l(cellnum(i),tstamp(i));
    wcell(cellid,tstamp(i))=w(cellnum(i),tstamp(i));
    acell(cellid,tstamp(i))=a(cellnum(i),tstamp(i));
    B{cellid,tstamp(i)}=boun{cellnum(i),tstamp(i)};
    pixels{cellid,tstamp(i)}=pxls{cellnum(i),tstamp(i)};
    mlines{cellid,tstamp(i)}=mline{cellnum(i),tstamp(i)};
    total_vols(cellid,tstamp(i))=total_vol(cellnum(i),tstamp(i)); %%% new
    pcell(cellid,tstamp(i))=pole(cellnum(i),tstamp(i));
end

%Filter out cells found at only one or two time points
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
total_vols(delind,:)=[];%%%
[ncells,~]=size(lcell);

lcell(lcell==0)=NaN;
wcell(wcell==0)=NaN;
acell(acell==0)=NaN;
pcell(pcell==0)=NaN;
total_vols(total_vols==0)=NaN;%%% New

[lm,ln]=size(l);
tmid=(time(2:end)+time(1:end-1))/2;

%Dimsionalize the variables
lcell=lcell*lscale;
wcell=wcell*lscale;
acell=acell*lscale^2;
total_vols=total_vols*lscale^3;

%Throw away cells that are too short, too fat, or too skinny
lcell(lcell<minL|wcell>maxW|wcell<minW)=NaN;
wcell(lcell<minL|wcell>maxW|wcell<minW)=NaN;
acell(lcell<minL|wcell>maxW|wcell<minW)=NaN;

%Throw away cells that grow too fast
ld=(lcell(:,2:end)-lcell(:,1:end-1));
ld=abs([zeros(ncells,1) ld]);
lcell(ld>1)=NaN;

wd=(wcell(:,2:end)-wcell(:,1:end-1)); %Edit
wd=abs([zeros(ncells,1) wd]); %Edit
wcell(wd>1)=NaN; %Edit

%Throwing out widths that vary to greatly % Dylan Edit
for i = 1:size(wcell,1)
    if nanmean(abs(diff(wcell(i,1:end)))) > 0.03
        wcell(i,1:end) = NaN;
    else wcell(i,1:end) = wcell(i,1:end);
    end
    if nanmean(abs(diff(total_vols(i,1:end)))) > 0.09
        total_vols(i,1:end) = NaN;
    else total_vols(i,1:end) = total_vols(i,1:end);
    end
end

ad=(acell(:,2:end)-acell(:,1:end-1)); %Edit
ad=abs([zeros(ncells,1) ad]); %Edit
acell(ad>1)=NaN; %Edit

%Calculate the elongation rate
deltat=time(2:end)-time(1:end-1);
v=(lcell(:,2:end)-lcell(:,1:end-1))./((lcell(:,1:end-1)+lcell(:,2:end))/2);
av=(acell(:,2:end)-acell(:,1:end-1))./((acell(:,1:end-1)+acell(:,2:end))/2);
wv=(wcell(:,2:end)-wcell(:,1:end-1))./((wcell(:,1:end-1)+wcell(:,2:end))/2); %Edit
for i=1:ncells
    v(i,:)=v(i,:)./deltat;
    av(i,:)=av(i,:)./deltat;
    wv(i,:)=wv(i,:)./deltat; %Edit
end

%Throw away outliers and calculate the average width, and elongation rate across cells
v(isnan(v))=0;
av(isnan(av))=0;
wv(isnan(wv))=0; %EDIT

for t=1:T-1
    vav(t)=mean(nonzeros(v(:,t)));
    vstd(t)=std(nonzeros(v(:,t)));
    %wvav(t)=mean(nonzeros(wv(:,t))); %EDIT
    %wvstd(t)=std(nonzeros(wv(:,t))); %EDIT
end
vavm=ones(ncells,1)*vav;
vstdm=ones(ncells,1)*vstd;
%wvavm=ones(ncells,1)*wvav; %EDIT
%wvstdm=ones(ncells,1)*wvstd; %EDIT

inddel=abs(v-vavm)>2*vstdm&vstdm~=0;
%winddel=abs(wv-wvavm)>2*wvstdm&wvstdm~=0; %EDIT?

v(inddel)=0;
av(inddel)=0;
wv(inddel)=0; %EDIT
lcell(inddel)=0;
acell(inddel)=0;
wcell(inddel)=0;
ew(inddel)=0;


for t=1:T
    wav(t)=nanmean(nonzeros(wcell(:,t)));
    wstd(t)=std(nonzeros(wcell(:,t)));
    wste(t)=wstd(t)./length(nonzeros(wcell(:,t)));
    %ewav(t)=nanmean(nonzeros(ew(:,t)));
    %ewstd(t)=std(nonzeros(ew(:,t)));
    %ewste(t)=ewstd(t)./length(nonzeros(ew(:,t)));
end
for t=1:T-1
    ndp(t)=length(nonzeros(v(:,t)));
    wvav(t)=mean(nonzeros(wv(:,t))); %EDIT
    wvstd(t)=std(nonzeros(wv(:,t))); %EDIT
    wvste(t)=wvstd(t)/sqrt(ndp(t)); %EDIT
    vav(t)=mean(nonzeros(v(:,t)));
    vstd(t)=std(nonzeros(v(:,t)));
    vste(t)=vstd(t)/sqrt(ndp(t));
    avav(t)=mean(nonzeros(av(:,t)));
    avstd(t)=std(nonzeros(av(:,t)));
    avste(t)=avstd(t)/ndp(t);
end

v(v==0)=NaN;
av(av==0)=NaN;
wv(wv==0)=NaN; %EDIT
lcell(lcell==0)=NaN;
wcell(wcell==0)=NaN;
acell(acell==0)=NaN;
ew(ew==0)=NaN;
       

Leff=EffectiveLength(tmid,vav);
Lsmooth=(Leff-movingaverage(Leff,12))./movingaverage(Leff,12);

Weff = EffectiveLength(tmid,wvav); %EDIT
Wsmooth=(Weff-movingaverage(Weff,12))./movingaverage(Weff,12);

Aeff = EffectiveLength(tmid,avav); %EDIT
Asmooth=(Aeff-movingaverage(Aeff,12))./movingaverage(Aeff,12);

% Getting volumes from formula 
cell_radius = wcell./2;
cell_side_length = lcell - wcell;
cell_volume1 = pi.*cell_radius.^2;
cell_volume2 = (4/3).*cell_radius + cell_side_length;
cell_volume = cell_volume1.*cell_volume2;
vcell = cell_volume;

surface_area1 = 2*pi.*cell_radius;
surface_area2 = (2.*cell_radius) + cell_side_length;
surface_area = surface_area1.*surface_area2;
scell = surface_area;

circumference = surface_area1; 
ccell = circumference;

tv=(total_vols(:,2:end)-total_vols(:,1:end-1))./((total_vols(:,1:end-1)+total_vols(:,2:end))/2); %Edit 
vv=(vcell(:,2:end)-vcell(:,1:end-1))./((vcell(:,1:end-1)+vcell(:,2:end))/2); %Edit
sv=(scell(:,2:end)-scell(:,1:end-1))./((scell(:,1:end-1)+scell(:,2:end))/2); %Edit
cv=(ccell(:,2:end)-ccell(:,1:end-1))./((ccell(:,1:end-1)+ccell(:,2:end))/2); %Edit
for i=1:ncells
    tv(i,:)=tv(i,:)./deltat;
    vv(i,:)=vv(i,:)./deltat; %Edit
    sv(i,:)=sv(i,:)./deltat; %Edit
    cv(i,:)=cv(i,:)./deltat; %Edit
end

%Throw away outliers and calculate the average width, and elongation rate across cells
tv(isnan(tv))=0; %EDIT
vv(isnan(vv))=0; %EDIT
sv(isnan(sv))=0; %EDIT
cv(isnan(cv))=0; %EDIT

tv(inddel)=0; %EDIT
vv(inddel)=0; %EDIT
sv(inddel)=0; %EDIT
cv(inddel)=0; %EDIT

for t=1:T-1
    tvav(t)=mean(nonzeros(tv(:,t)));
    
    vvav(t)=mean(nonzeros(vv(:,t))); %EDIT
    vvstd(t)=std(nonzeros(vv(:,t))); %EDIT
    vvste(t)=vvstd(t); %EDIT
    
    svav(t)=mean(nonzeros(sv(:,t))); %EDIT
    svstd(t)=std(nonzeros(sv(:,t))); %EDIT
    svste(t)=svstd(t)/(ndp(t)); %EDIT
    
    cvav(t)=mean(nonzeros(cv(:,t))); %EDIT
    cvstd(t)=std(nonzeros(cv(:,t))); %EDIT
    cvste(t)=cvstd(t)/sqrt(ndp(t)); %EDIT
end

tv(tv==0)=NaN; %EDIT
vv(vv==0)=NaN; %EDIT
sv(vv==0)=NaN; %EDIT
cv(vv==0)=NaN; %EDIT

TVeff = EffectiveLength(tmid,tvav); %EDIT
TVsmooth=(TVeff-movingaverage(TVeff,12))./movingaverage(TVeff,12);

Veff = EffectiveLength(tmid,vvav); %EDIT
Vsmooth=(Veff-movingaverage(Veff,12))./movingaverage(Veff,12);
Seff = EffectiveLength(tmid,svav); %EDIT
Ssmooth=(Seff-movingaverage(Seff,12))./movingaverage(Seff,12);
Ceff = EffectiveLength(tmid,cvav); %EDIT
Csmooth=(Ceff-movingaverage(Ceff,12))./movingaverage(Ceff,12);


%Plot data
figure(1), title('Cell Length vs. Time')
clf
hold on
for i=1:ncells
    %indx=isnan(lcell(i,:))~=1;
    %plot(time(indx)/60,lcell(i,indx))
    %plot(time/60,lcell(i,:))
    plot(1:size(lcell,2),lcell(i,:))
end
xlabel('Time (min)')
ylabel('Length (\mum)')
fig2pretty

figure, title('Effective Length vs. Time')
hold on 
plot(1:length(Leff),Leff)
%legend('Leff','Weff','Aeff','Volume from capsule formula','Veff')
xlabel('Time (min)')
ylabel('Length (\mum)')
fig2pretty

figure(3), title('Cell Width vs. Time')
hold on
for i=1:ncells
    plot(time/60,wcell(i,:))
       
end
plot(time/60,wav,'-r','LineWidth',20)
xlabel('Time (min)')
ylabel('Width (\mum)')
fig2pretty

figure(4), title('Cell Width vs. Time')
hold on
for i=1:ncells
    plot(time/60,wcell(i,:))
       
end
plot(time/60,wav,'-r','LineWidth',20)
xlabel('Time (min)')
ylabel('Width (\mum)')
fig2pretty

figure (5), title('Average Cell Width vs. Time')
% y = movingaverage(y,10); ??? 
plot(time/60, wav,'-r','LineWidth',4)
plot(1:length(wav), wav,'-r','LineWidth',4)
xlabel('Time (min)')
ylabel('Width (\mum)')
fig2pretty

figure(6), title('Elongation Rate vs. Time')
hold on
%plot(tmid/60,vav*60,'-r')
plot(1:length(vav*60),vav*60,'-r')
xlabel('Time (min)')
ylabel('Elongation Rate (min^{-1})')
fig2pretty

% Part 2: Lineage Tracking
% %view binary images
% for i=flip(1:27)
% imshow(labels(:,:,i),[0 2])
% pause
% end

%%
%Get coordinates of cells
for i = flip(1:T)
%for i=27
numofcells(i)=max(max(labels(:,:,i)));% find number of cells
mxnumofcells=max(numofcells);
    for j=1:mxnumofcells
        [r1,c1] = find(bwlabel(labels(:,:,i))==j);   %Find coordinates of cells in frame T

        r1c1{j,i}=[r1 c1]; % concatenate XY coords frame T 
        
    end
end

for i = flip(1:T)
    for j=1:mxnumofcells
        if isempty(r1c1{j,i})==1
            r1c1{j,i} = [NaN NaN];
        end
    end
end

%Align cells from frame to frame, correct r1c1 to track lineage
for i=T
    for j=1:mxnumofcells
        cell_id(j,i)=j; %makes cell id vector
        r2c2{j,i}=r1c1{j,i};%Seed r2c2
    end
end

%Tracks cells backward through time,  
for i = flip(2:T)
    for j=1:mxnumofcells
        if isnan(r2c2{j,i})==1 %if NaN, find previous nonNaN id matrix for tracking
            index=cellfun('length', r2c2(j,:));
            nindex=index>1;
            ncols=min(find(nindex));
            Lia1=ismember(r2c2{j,ncols},r1c1{j,i-1},'rows');
            if 100*sum(Lia1)/length(r1c1{j,i-1})>50==1
                cell_id(j,i-1)=j;
                r2c2{j,i-1}=r1c1{j,i-1};
            else
                cell_id(j,i-1)=NaN;
                r2c2{j,i-1}=NaN;
            end
        else
            nidx=cellfun('length', r2c2(j,:)); %Normal tracking backward 
            nnidx=nidx>1;
            ncols=min(find(nnidx));
            Lia1=ismember(r1c1{j,i-1},r2c2{j,ncols},'rows');
                if 100*sum(Lia1)/length(r2c2{j,ncols})>50==1
                            cell_id(j,i-1)=j;
                            r2c2{j,i-1}=r1c1{j,i-1};
                    else
                        r2c2{j,i-1}=NaN;
                end
        end
     end
end

% Fill NaNs tracking backward
for i = flip(1:T)
for j=1:mxnumofcells
    for k=1:mxnumofcells
        if isnan(r2c2{j,i})==1
            index=cellfun('length', r2c2(j,:));
            nindex=index==1;
            ncols=max(find(nindex));
            Lia1=ismember(r2c2{j,ncols+1},r1c1{k,i},'rows');
            if 100*sum(Lia1)/length(r1c1{k,ncols+1})>50==1
                cell_id(j,i)=k;
                r2c2{j,i}=r1c1{k,i};
            else
                %cell_id(j,i-1)=NaN;
                %r2c2{j,i-1}=NaN;
                end
        end
     end
end
end

%% Organize lcell with lineage tracking
%Track cells frame to frame
tracks=zeros(size(im));
[r2,c2] = find(bwlabel(labels(:,:,end))>=1);
linind=sub2ind(size(im),r2,c2);%new linind for tracking binary image
tracks(linind)=1;

[tracksL(:,:,i),ncells]=bwlabel(tracks);

nlcell=zeros(ncells,T);

for i=1:lcents
    nlcell(cellnum(i),tstamp(i))=l(cellnum(i),tstamp(i));
end

%scale lcell
nlcell=nlcell*lscale;

% Find division time point for dividing lcell
% - Still need to optimize for larger data sets
for k=1:size(cell_id,1)-1
    if  find(cell_id(k+1,:)==k)
        split=find(cell_id(k+1,:)==k);
        lcell(k,1:split)=lcell(k,1:split)/2;
        lcell(k+1,1:split)=lcell(k,1:split);
    end
end

nlcell(nlcell == 0) = NaN;

figure
plot(1:218,nlcell(:,:))
% hold on 
% plot(1:27,lcell(2,:))
    
%% 
% cd /Users/dylanfitzmaurice/Documents/MATLAB/data
% save([basename '_BT'])
% save([basename '_BTlab'],'labels','labels2','-v7.3')

% For saving detergent data in parts
%save([basename '_BT_L1'])
%save([basename '_BTlab_L1'],'labels','labels2','-v7.3')

beep