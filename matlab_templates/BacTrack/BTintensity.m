%BTintensity.m  
%Calculates the intensity of a flourescent tracer dye that was imaged 
%simultaneously with a phase stack that was analyzed with BacTrack.m.
%INSTRUCTIONS FOR USE:
%Run BacTrack.m on a phase image stack.  Save the fluorescent image stack 
%(which images a tracer dye) in a separate directory.
%
%INPUT:
%basename: name of the image stack (from BacTrack).
%dirname: full path name of directory in which fluorescent image stack is
%           saved.
%Cmin: concentration of medium without tracer dye.
%Cmax: concentration of medium with tracer dye.
%
%OUTPUT:
%conc: vector of length T with concentration trace.
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
Cmax=500;
Cmin=0;
recrunch=0;
basename='03172021_Exp2_colony1';%Name of the image stack, used to save file.
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03172021_analysis/' basename];
dirname=[[filename '/' basename '_647/' basename '_full']];%Directory that the image stack is saved in.
%basenameData='11092020'
dirnameData=[filename '/' basename '_phase/' basename '_figures/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(dirnameData);
filelist = dir([basename '_BTphase.mat']);
%load data from BacTrack.
load([filelist(1).name],'T','time','tmid','vav','tpoints','lcell','time2','vstd','ncells');
%path(dirname,path);
%time2=tpoints(end,1:T);
curdir=cd;
cd(dirname);
directory=dir('*.tif');
T=length(directory); 
%Tp=T;
int=zeros(1,T);
curdir=cd;
cd(dirname);
directory=dir('*.tif');
cd(curdir);
for t=1:T
    t
    %load image
    imagename=directory(t).name;
    im=imread(imagename);
    %int(t)=mean(nonzeros(im));
    Idouble = im2double(im);
    avg(t) = 50*mean2(Idouble); %why multipy by 50?
end
%calculate Intencity Rate
deltat=time(2:end)-time(1:end-1);
Int_rate=(avg(2:end)-avg(1:end-1))./((avg(1:end-1)+avg(2:end))/2);
Int_rate=Int_rate(:)./deltat;
%Average
%for t=1:T-1
%    intent(t)=mean(nonzeros(Int_rate(t)));
%end
%hold on
figure
plot(tmid,Int_rate)
xlabel('time (S)')
ylabel('Intensity Rate (s^{-1})')
fig2pretty


%plot results
% figure
% ciplot((vav-vstd)*3600,(vav+vstd)*3600,tmid,[0.75 0.75 1]),hold on
% plot(tmid,vav*3600,'-r','LineWidth',1)
% ylabel('Growth Rate (s^{-1})')
% g = get(gca,'Position');
% g(1:2) = g(1:2) + (1-0.9)/2*g(3:4);
% g(3:4) = 0.9*g(3:4);
% set(gca,'Position',g);
% box off
% apos=get(gca,'Position');
% xtic=get(gca,'Xtick');
% %ytic=get(gca,'Xtick');
% lims=get(gca,'Xlim');
% axes('Position',apos,'XTick',xtic,'Color','none','YAxisLocation','right',...
%     'Xlim',lims);
% hold on
% plot(time2,avg,'LineWidth',2)
% xlabel('time (S)')
% ylabel('C (mM)')
% fig2pretty
% %figure,plot(conc,lcell)
% %fig2pretty
% %save
% %load([basename '_BT'])
% figure, title('Cell Length vs. Time')
% clf
% hold on
% for i=1:ncells  
%     lcell(i,:)=movingaverage2(lcell(i,:),3);
%     %indx=isnan(lcell(i,:))~=1;
%     %indx=find(indx);
%     %plot(time(indx),lcell(i,indx))
%     plot(time(1:end),lcell(i,1:end)) 
% end
% hold on
% plot(time,avg,'--k','LineWidth',2)
% xlabel('time (S)')
% ylabel('length (\mum)')
% fig2pretty
% %plot results
% % figure
% % plot(time(48:80),lcell(i,48:80),'-r','LineWidth',2)
% % ylabel('Growth Rate (s^{-1})')
% % g = get(gca,'Position');
% % g(1:2) = g(1:2) + (1-0.9)/2*g(3:4);
% % g(3:4) = 0.9*g(3:4);
% % set(gca,'Position',g);
% % box off
% % apos=get(gca,'Position');
% % xtic=get(gca,'Xtick');
% % lims=get(gca,'Xlim');
% % axes('Position',apos,'XTick',xtic,'Color','none','YAxisLocation','right',...
% %     'Xlim',lims);
% % hold on
% % plot(time,conc)
% % xlabel('time (s)')
% % ylabel('C (mM)')
% % figure
% % plot(time(58:66),lcell(i,58:66),'-r','LineWidth',2)
% % ylabel('Growth Rate (s^{-1})')
% % g = get(gca,'Position');
% % g(1:2) = g(1:2) + (1-0.9)/2*g(3:4);
% % g(3:4) = 0.9*g(3:4);
% % set(gca,'Position',g);
% % box off
% % apos=get(gca,'Position');
% % xtic=get(gca,'Xtick');
% % lims=get(gca,'Xlim');
% % axes('Position',apos,'XTick',xtic,'Color','none','YAxisLocation','right',...
% %     'Xlim',lims);
% % hold on
% % plot(time,conc)
% % xlabel('time (s)')
% % ylabel('C (mM)')
% %fig2pretty
% %save([basename '_BTint'])