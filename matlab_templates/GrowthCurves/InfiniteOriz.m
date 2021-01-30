%01062019 Paola 
%code to analyze growth curves 

filename=('growth_curves 01042020');
fullpath=('C:\Users\paola\Desktop\Growth Curves\01032020');
xlRange='B50:EL145';
nwells=96
T=141

cond1=[1:12:25];
blank1=[39];
cond2=[2:12:26];
blank2=[38];
cond3=[3:12:27];
blank3=[39];
cond4=[4:12:28];
blank4=[40];
cond5=[5:12:29];
blank5=[41];
cond6=[6:12:30];
blank6=[42];
cond7=[7:12:31];
blank7=[43]
cond8=[8:12:32];
blank8=[44]
cond9=[9:12:33]
blank9=[45]
cond10=[10:12:34]
blank10=[46]
cond11=[11:12:35]
blank11=[47]
cond12=[12:12:36]
blank12=[45]
cond13=[49:12:73];
blank13=[86];
cond114=[50:12:74];
blank14=[86];
cond15=[51:12:75];
blank15=[87];
cond16=[4:12:76];
blank6=[88];
cond17=[5:12:77];
blank17=[89];
cond18=[6:12:78];
blank8=[90];
cond19=[7:12:79];
blank9=[91]
cond20=[8:12:80];
blank20=[92]
cond121=[9:12:81]
blank21=[93]
cond22=[10:12:82]
blank22=[94]
cond23=[11:12:83]
blank23=[95]
cond24=[12:12:84]
blank24=[92]


condlab={'Wt','Wt1.1%','Wt1.2%','Wt1.3%','Wt1.4%','Wt1.5%','Wt1.6%','Wt1.7%','Wt1.8%','Wt1.9%','Wt2%','Wt1%','MurA','MurA1.1%','MurA1.2%','MurA1.3%','MurA1.4%','MurA1.5%','MurA1.6%','MurA1.7%','MurA1.8%','MurA1.9%','MurA2%','MurA1%'}
%upload file 

GCtable=xlsread(filename,xlRange);
OD=GCtable;
tscale=GCtable(1,:)/600

%parse data

WellInd={cond1,blank1;
  cond2,blank2;
  cond3,blank3;
  cond4,blank4;
  cond5,blank5;
  cond6,blank6;
  cond7,blank7;
  cond8,blank8;
  cond9,blank9;
  cond10,blank10;
  cond11,blank11;
  cond12,blank12;
  cond13,blank13;
  cond14,blank14;
  cond15,blank15;
  cond16,blank16;
  cond17,blank17;
  cond18,blank18;
  cond19,blank19;
  cond20,blank20;
  cond21,blank21;
  cond22,blank22;
  cond23,blank23;
  cond24,blank24};
  

[Ncond ~]=size(WellInd); 
wellODi=cell(Ncond,1);
blankODi=cell(Ncond,1);
OD_av=zeros(Ncond,T);
OD_av_bl=zeros(Ncond,T);
OD_norm=zeros(Ncond,T);

 
for i=1:1:Ncond
    wellODi{i,:}=OD(WellInd{i,1},:)
    blankODi{i,:}=OD(WellInd{i,2},:)
end 

for i=1:Ncond
    OD_av(i,:)=mean(wellODi{i},1)
    OD_av_bl(i,:)=mean(blankODi{i},1)
    OD_norm=(OD_av-OD_av_bl)
    %OD_st=std(OD_norm)
end
   
%calculate growthrate  for i=1:Ncond

time=([0:T-1]*600)
dOD=(OD_norm(:,2:2:end))-(OD_norm(:,1:2:end-1))
eOD=dOD./(OD_norm(:,1:2:end-1))/600
eODsmooth=movingaverage(eOD,10)
tmid=(time(2:2:end)+time(1:2:end-1))/2
tim=tmid./3600
time2=([0:T-1]*600)/3600


figure
plot(tim,eODsmooth*3600)
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
legend(condlab)
fig2pretty
saveas(gcf,'/Users/paola/Desktop/growth curves/09052019')
 
figure 
plot(time2,OD_norm)
xlabel('Time (h)')
ylabel('OD(AU)')
legend(condlab)
fig2pretty
saveas(gcf,'/Users/paola/Desktop/growth curves/09052019')
