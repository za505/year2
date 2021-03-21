%PlasmolysisMeasure.m
%Based on BacTrack.m and PlasmolysisCounter.m
%Customized for B. subtilis.

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
basename=["02102021_pbs100_colony1", "02102021_pbs100_colony2","02102021_pbs100_colony3"];%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/pbs_analysis2/'];%Directory that the image stack is saved in.
savedir=[dirname 'finalAnalysis'];%Directory to save the output .mat file to.

lscale=0.08;%%Microns per pixel.
thresh=0;%For default, enter zero.
IntThresh=13000;%Threshold used to enhance contrast. Default:35000
dr=1;%Radius of dilation before watershed %default: 1
sm=2;%Parameter used in edge detection %default: 2
minL=2;%Minimum cell length, default = 2
minW=0.2;%Minimum cell width, default = 0.2
maxW=1.5;%Maximum cell width, default = 1.5
minA=200;%Minimum cell area. default: 50 (also use 200)
crunch=1;%add counts to data and calculate # plasmolysis events/micron 0=No, 1=Yes.
checkhist=0;%Display image histogram? 0=No, 1=Yes.
B=length(basename); %number of main directories to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if crunch==1
    
   %load previously processed data
   x=load([dirname 'plasmolysisCounter/finalAnalysis/BTplasmolysis.mat'])
        
    %load the data you just saved
    y=load([dirname 'plasmolysisCounter/finalAnalysis/BTplasmolysis.mat'])
    
    %load count data
    filename=[ dirname 'README.xlsx'];
    data=xlsread(filename);
    
    %for this analysis, we'll have to average the plasmolysis events for
    %diff colonies of the same experiment
    
    pis=[mean(data(1:3, 9)), mean(data(4:6, 9)), mean(data(7:9, 9)), mean(data(10:12, 9))];
    pil=[mean(data(1:3, 10)), mean(data(4:6, 10)), mean(data(7:9, 10)), mean(data(10:12, 10))];
    pit=[mean(data(1:3, 11)), mean(data(4:6, 11)), mean(data(7:9, 11)), mean(data(10:12, 11))];
    pfs=[mean(data(1:3, 12)), mean(data(4:6, 12)), mean(data(7:9, 12)), mean(data(10:12, 12))];
    pfl=[mean(data(1:3, 13)), mean(data(4:6, 13)), mean(data(7:9, 13)), mean(data(10:12, 13))];
    pft=[mean(data(1:3, 14)), mean(data(4:6, 14)), mean(data(7:9, 14)), mean(data(10:12, 14))];
    
    ltotal=[mean(x.ltotal(1:3)), mean(x.ltotal(4:6)), mean(x.ltotal(7:9)), mean(y.ltotal(1:3))];
    
    %calculate plasmolysis events/micron immediately after shock
    pism=pis./ltotal;
    pilm=pil./ltotal;
    pitm=pit./ltotal;
    pfsm=pfs./ltotal;
    pflm=pfl./ltotal;
    pftm=pft./ltotal;

    %plot
    tiledlayout(1, 3)
    
    %nexttile
    %X=categorical({'immediately after shock', '10 min after shock'});
    %X=reordercats(X, {'immediately after shock', '10 min after shock'});
    X=[0,10];
    createfigure(X, [pitm; pftm])
    ylabel('Plasmolysis Events/Micron Post-Hyperosmotic Shock')
    xlabel({'Time (minutes after hypershock)'})
    legend('pbs 0 min','pbs 1 min', 'pbs 10 min', 'pbs 100 min')
    
    %nexttile
    %X=categorical({'immediately after shock', '10 min after shock'});
    %X=reordercats(X, {'immediately after shock', '10 min after shock'});
    %bar(X, [pism; pfsm])
    %ylabel('Septal Plasmolysis Events/Micron Post-Hyperosmotic Shock')
    %legend('pbs 0 min','pbs 1 min', 'pbs 10 min', 'pbs 100 min')
    
    %nexttile
    %X=categorical({'immediately after shock', '10 min after shock'});
    %X=reordercats(X, {'immediately after shock', '10 min after shock'});
    %bar(X, [pilm; pflm])
    %ylabel('Lateral Plasmolysis Events/Micron Post-Hyperosmotic Shock')
    %legend('pbs 0 min','pbs 1 min', 'pbs 10 min', 'pbs 100 min')

else
    nc=zeros(1,B); %we're interested in the number of cells counted for the controls
    allcentroids=[];
    cellnum=[];

    %Pre-allocate matrices
    a=zeros(1,B);
    w=zeros(1,B);
    l=zeros(1,B);
    DS=zeros(1,B);
    boun=cell(1,B); %we only need the boundaries for the controls
    pole=zeros(1,B);
    mline=cell(1,B);
    
    ltotal=zeros(1, B);
    atotal=zeros(1, B);
    
    for b=1:B
        
        base=char(basename(b))
        
        %Determine number of frames
        curdir=cd;
        cd([dirname base '/' base '_originals']);
        directory=dir('*.tif');
        T=length(directory);

        cd(curdir);
        path(dirname,path)

        %Load first image
        imagename=[directory(1).folder '/' directory(1).name];
        im=imread(imagename);
        [imM,imN]=size(im);
        labels=zeros(imM,imN,T);
        labels2=zeros(imM,imN,T);
    
        for t=1:T
            t
            %Load image
            imagename=[directory(t).folder '/' directory(t).name];

            im=imread(imagename);
            [imM,imN]=size(im);
            
           if t==1
               
               %%First, our inputs
                thresh=0.9;%Parameter used to threshold image, to estimate this, use checkhist.
                dr=6;%Radius to dilate posts by

                tscale=60;%Frame rate in seconds
                lscale=0.098;%Length scale in microns/pixel
                smooth=20;%Number of frames over which to smooth growth rate profiles

                vis=1;%Visualize tracking %vis=1 is tracking, vis=0 no tracking
                checkhist=0;%Look at image histogram
                checkposts=0;%Check microfluidic posts
                recrunch=0;%Load previously cruched data
                
                %%Now for the rest of the FluoColonyTrack code
                %Load first image
                im=imread(imagename);
                im2=imread(imagename);

                %Determine Background
                figure,imshow(im2,[]), hold on, title('Select Background')
                k=waitforbuttonpress;
                set(gcf,'Pointer','cross')
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
                set(gcf,'Pointer','cross')
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

                %Load images
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

%                 %find region props
%                 stats=regionprops(col_edges,'Area','MeanIntensity');
%                 idx=find([stats.Area]>minA&[stats.Area]<1e5&[stats.MeanIntensity]>IntThresh);
%                 
%                 stats = idx;
%                 cc=bwconncomp(ed3,4);
%                 stats=regionprops(cc,imc,'Area','MeanIntensity');
%                 idx=find([stats.Area]>minA&[stats.Area]<1e5&[stats.MeanIntensity]>3e4);
%                 ed4=ismember(labelmatrix(cc),idx);

                %Find cell areas and centroids
                [P,bw]=bwboundaries(bw,4,'noholes');
                stats=regionprops(bw,'Area','Centroid','PixelIdxList');
                idx=find([stats.Area]>minA&[stats.Area]<1e5);
                P=P(idx, :);
                stats=stats(idx, :);
                nc(b)=length(stats);
                
%                 L=bwlabel(bw);    
%                 labels(:,:,1)=L;
%                 labels2(:,:,1)=bw;

                %nc(b)=length(P);
                areas=[stats.Area];
                cents=cat(1,stats.Centroid);
                a(nc(b),b)=0;
                a(1:nc(b),b)=[stats.Area]';
                centroids=cents;

                %Calculate smooth cell contours
                for n=1:nc(b)

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
                    DS(n,b)=S(end)/(ls-1);
                    Sn=(0:DS(n,1):S(end));
                    nx=spline(S,px,Sn);
                    ny=spline(S,py,Sn);

                    boun{n,b}=[nx',ny'];
                    pxls{n,b}=stats(n).PixelIdxList;

                end

                allcentroids=[allcentroids;centroids];
                cellnum=[cellnum;(1:nc(b))'];

                %Calculate cell length, width, etc.
                for n=1:nc(b)
                    if isempty(boun{n,b}) == 0
                     X=boun{n,b}(:,1);
                     Y=boun{n,b}(:,2);   

                     [sX,~]=size(X);

                     %Find poles
                     [X,Y,pole(n,b)]=polefinder(X,Y);

                     %Create mesh
                     npts=min(pole(n,b),sX-pole(n,b)+1);
                     S=(0:DS(n,b):(sX-1)*DS(n,b));

                     s1=(0:S(pole(n,b))/(npts-1):S(pole(n,b)));
                     s2=(S(pole(n,b)):(S(end)-S(pole(n,b)))/(npts-1):S(end));
                     xc1=spline(S(1:pole(n,b)),X(1:pole(n,b)),s1);
                     xc2=spline(S(pole(n,b):end),X(pole(n,b):end),s2);
                     yc1=spline(S(1:pole(n,b)),Y(1:pole(n,b)),s1);
                     yc2=spline(S(pole(n,b):end),Y(pole(n,b):end),s2);
                     xc2=fliplr(xc2);
                     yc2=fliplr(yc2);

                     %Calculate midline
                     mline{n,b}=[(xc1+xc2)'/2,(yc1+yc2)'/2];
                     dxy=diff(mline{n,b}).^2;
                     dl=sqrt(dxy(:,1)+dxy(:,2));
                     l(n,b)=sum(dl);

                     %Calculate width
                     ls=[0 cumsum(dl)'];
                     [~,mpos1]=min(abs(ls/l(n,b)-0.25));
                     [~,mpos2]=min(abs(ls/l(n,b)-0.75));

                     widths=sqrt((xc1-xc2).^2+(yc1-yc2).^2);
                     w(n,b)=(widths(mpos1)+widths(mpos2))/2; 
                    else
                        nc(b) = nc(b) - 1;
                    end %end if statement
                end

                %Dimsionalize the variables
                lcell=l*lscale;
                wcell=w*lscale;
                acell=a*lscale^2;

                %Throw away cells that are too short or too fat or too skinny
                lcell(lcell<minL|wcell>maxW|wcell<minW)=NaN;
                wcell(lcell<minL|wcell>maxW|wcell<minW)=NaN;
                acell(lcell<minL|wcell>maxW|wcell<minW)=NaN;

                %Calculate total length of all cells
                lcell(isnan(lcell))=0;
                acell(isnan(acell))=0;
                ltotal(1, b)=sum(nonzeros(lcell(:,b)));
                atotal(1, b)=sum(nonzeros(acell(:,b)));
                
           else
               %re-declare inputs
                lscale=0.08;%%Microns per pixel.
                thresh=0;%For default, enter zero.
                IntThresh=13000;%Threshold used to enhance contrast. Default:35000
                dr=1;%Radius of dilation before watershed %default: 1
                sm=2;%Parameter used in edge detection %default: 2
                minL=2;%Minimum cell length, default = 2
                minW=0.2;%Minimum cell width, default = 0.2
                maxW=1.5;%Maximum cell width, default = 1.5
                minA=200;%Minimum cell area. default: 50 (also use 200)
                crunch=0;%add counts to data and calculate # plasmolysis events/micron 0=No, 1=Yes.
                checkhist=0;%Display image histogram? 0=No, 1=Yes.
                B=length(basename); %number of main directories to analyze
                
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
            
               cd([dirname base '/' base '_markups']);
               figure
               imshow(im)
               hold on
               for k=1:nc(b)
                   if lcell(k, 1) ~= 0
                    plot(boun{k,b}(:,1),boun{k,b}(:,2),'-r')
                   end
               end
               saveas(gcf, [directory(t).name])

              pause
              close all
           end %end of if t==1 statement
        end %end of t=1:T loop
       
        cd(savedir)
        save([base 'BTtemp.mat'])

    end %end of b=1:B

    %Save data
    cd(savedir)
    save('BTpbs100.mat')
        
    if crunch==-1
        
        %load previously processed data
        y=load([savedir '/BTplasmolysis.mat'])
        
        %load the data you just saved
        x=load([savedir '/BTpbs100.mat'])
    
        % Check to see that both files contain the same variables
        vrs = fieldnames(x);
        if ~isequal(vrs,fieldnames(y))
            error('Different variables in these MAT-files')
        end
        
        % Concatenate data
        for k = 1:length(vrs)
            x.(vrs{k}) = [x.(vrs{k});y.(vrs{k})];
        end
    
        % Save result in a new file
        cd(savedir)
        save('BTplasmolysis_final.mat')
       
    end %end of crunch=-1
    
end %end of if crunch=1