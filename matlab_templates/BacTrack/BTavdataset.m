%BTavplotmultiple.m
%Collates and plots multiple growth rate traces tracked with BacTrack.m
%
%INSTRUCTIONS FOR USE:
%Track the data sets interest with BacTrack.m before using.
%
%INPUT:
%basename: root of filenames for files to be compared.  That is, files to
%be collated should be named basename*, where * is a suffix denoting the
%individual file.

clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input 
basename=["03172021_Exp3_colony1","03172021_Exp3_colony2","03172021_Exp3_colony3"];
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03172021_analysis'];
savename=['03172021_Exp3'];
sfac=1;
B=length(basename); %number of main directories to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b=1:B
    base=char(basename(b))
    cd([filename '/' base '/' base '_phase/' base '_figures/'])
    filelist(b) = dir('*_BTphase.mat');
end

%filelist=filelist(end);
lfl=length(filelist);

cmap=jet;
lcmap=floor(64/lfl);

%load([filelist{1} '_BTconn'],'T')
load([filelist(B).name],'T')
vt=[];
lcellt=[];
wcellt=[];
atotalt=zeros(1,T);
ltotalt=zeros(1,T);

for ii=1:lfl
    ii
    base=char(basename(ii))
    cd([filename '/' base '/' base '_phase/' base '_figures/'])
    %load([filelist{ii} '_BTconn'],'v','tmid','T','atotal','ltotal','vav')
    load([filelist(ii).name],'v','tmid','T','vav','B','im','lcell','wcell','time','lscale','tscale')
    tscale
    lcell(isnan(lcell))=0;
    for t=1:T
        ltotal(t)=sum(nonzeros(lcell(:,t)));
    end
    
    vt=[vt;v];
    lcellt=[lcellt;lcell];
    wcellt=[wcellt;wcell];
    ltotalt=ltotalt+ltotal;
    
end
vt=vt;
vt(isnan(vt))=0;
wcellt(isnan(wcellt))=0;

vavt=zeros(1,T-1);
for t=1:T-1
    vavt(t)=mean(nonzeros(vt(:,t)));
    vstd(t)=std(nonzeros(vt(:,t)));
    vnum(t)=length(nonzeros(vt(:,t)));
    vste(t)=vstd(t)/sqrt(vnum(t)*sfac);
    wavt(t)=mean(nonzeros(wcellt(:,t)));
end

%filelist2 = dir([basename '*_BTint.mat']);
%ct=[];
% 
% for ii=1:lfl
%     ii
%     %load([filelist{ii} '_BTconn'],'v','tmid','T','atotal','ltotal','vav')
%     load([filelist2(ii).name],'conc','Cmax','Cmin')
%     
%     ct=[ct;conc];
% end

%cav=mean(ct,1);

Leff=EffectiveLength(tmid,vavt);
Lsmooth=(Leff-movingaverage(Leff,12))./movingaverage(Leff,12);


%plot data
% figure, title('Strain Rate vs. Time')
% hold on
% plot(tmid,vavt*3600,'-r')
% xlabel('t(s)')
% ylabel('Strain Rate (s^{-1})')
% fig2pretty
% [X,Y]=ginput;

cd(filename);

figure, title('Effective Length vs. Time')
plot(tmid,Leff)
xlabel('t (s)')
ylabel('l_{eff} (\mum)')
fig2pretty
xline(0, '--k', 'LB + 647') %frame 1-16
xline(170, '--k', '*PBS + 5% detergent') %frame 17-28
xline(290, '--k', '*PBS + 647') %frame 29-41
xline(420, '--k', '*PBS + 647 + 20 mM NaCl') %frame 42+
ylim([1.5 3.5])
saveas(gcf, [savename,'_effLength.png'])

% [mvt,~]=size(vt);
% figure, title('Strain Rate vs. Time')
% hold on
% for i=1:mvt
%     plot(tmid,vt(i,:)*3600,'-r')
% end
% xlabel('t(s)')
% ylabel('Strain Rate (s^{-1})')
% fig2pretty
% 
% nint=length(X)/2;
% mpos1=zeros(nint,1);
% mpos2=zeros(nint,2);
% 
% for i=1:nint
%     [~,mpos1p(i)]=min(abs(tmid-X(2*i-1)));
%     [~,mpos2p(i)]=min(abs(tmid-X(2*i)));
% end
% 
% vavi2=vav;
% for j=1:nint
%      vavts(mpos1p(j):mpos2p(j))=movingaverage(vavt(mpos1p(j):mpos2p(j)),sfac);
%      vstds(mpos1p(j):mpos2p(j))=movingaverage(vstd(mpos1p(j):mpos2p(j)),sfac);
%      vstes(mpos1p(j):mpos2p(j))=movingaverage(vste(mpos1p(j):mpos2p(j)),sfac);
% end
% % %tmid=tmid/60;
% 
% lcellt(lcellt==0)=NaN;
[lct,~]=size(lcellt);
% 
% % figure,hold on,
% % % for i=1:lct
% % %     plot(tmid,vt(i,:))
% % % end
% % %ciplot((vavts-vstes)*3600,(vavts+vstes)*3600,tmid,[0.75 0.75 1])
% % plot(tmid,vavts*3600,'-r')
% % xlabel('t (s)')
% % ylabel('e (hr^{-1})')
% % ylims=get(gca,'YLim');
% % ylim([0 ylims(2)])
% % g = get(gca,'Position');
% % g(1:2) = g(1:2) + (1-0.9)/2*g(3:4);
% % g(3:4) = 0.9*g(3:4);
% % set(gca,'Position',g);
% % box off
% % apos=get(gca,'Position');
% % xtic=get(gca,'Xtick');
% % %ytic=get(gca,'Xtick');
% % lims=get(gca,'Xlim');
% % axes('Position',apos,'XTick',xtic,'Color','none','YAxisLocation','right',...
% %     'Xlim',lims);
% % hold on
% % plot(time,cav)
% % ylabel('C (mM)')
% % fig2pretty

figure,hold on
for i=1:lct
    %indx=isnan(lcellt(i,:))~=1;
    indx=lcellt(i,:)~=0;
    indx=find(indx);
    plot(time(indx),lcellt(i,indx))
    %plot(time(indx)/60,lcellt(i,indx))
    %plot(time(1:end),lcell(i,1:end)) 
end
xlabel('t (s)')
ylabel('Cell Wall Length (\mum)')
fig2pretty
xline(0, '--k', 'LB + 647') %frame 1-16
xline(170, '--k', '*PBS + 5% detergent') %frame 17-28
xline(290, '--k', '*PBS + 647') %frame 29-41
xline(420, '--k', '*PBS + 647 + 20 mM NaCl') %frame 42+
saveas(gcf, [savename,'_lengthTrace.png'])
save([savename '_BTads'])

% figure
% plot(time,ltotalt)
% xlabel('t (s)')
% ylabel('l (\mum)')
% fig2pretty
% 
% figure,hold on
% ciplot((vavts-vstds)*3600,(vavts+vstds)*3600,tmid,[0.75 0.75 1])
% plot(tmid,vavts*3600,'-r')
% xlabel('t (s)')
% ylabel('e (s^{-1})')
% title('Smoothed Elongation Rate')
% fig2pretty
% 
% figure,hold on
% ciplot((vavts-vstes)*3600,(vavts+vstes)*3600,tmid,[0.75 0.75 1])
% plot(tmid,vavts*3600,'-r')
% xlabel('t (s)')
% ylabel('e (s^{-1})')
% title('Smoothed Elongation Rate')
% fig2pretty
% 
% mean_l=nanmean(nanmean(lcellt(:,1:30)'));
% mean_w=nanmean(nanmean(wcellt(:,1:30)'));
% std_l=nanstd(nanmean(lcellt(:,1:30)'));
% std_w=nanstd(nanmean(wcellt(:,1:30)'));
% 
% save([basename '_BTads'])
