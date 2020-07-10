
%FCTCompile
%
%Paola Bardetti 
%6/21/19

%Compiles data from several positions tracked with FluoColonyTrack35

clear, close all

%filelist={'060419_P1';'060419_P2'};
filelist={'01162020_Series1_P1','01162020_Series2_P1','01162020_Series3_P1', '01162020_Series4_P1'};
basename='01162020_WTdownshift';%WT B.s. LB->LB/2
samples={'Colony1', 'Colony2', 'Colony3', 'Colony 4', 'Average'};
smooth=20;
folder=['C:\Users\zarin\Downloads\NYU\R2_Rojas Lab\growthRate_01162020_WTdownshift'];
%sz=360; %max lenght of vector

for i=1:length(filelist)
    load([filelist{i},'_FCT'],'Ncells','tmid','time','tscale','ColonyArea','eAsmooth','eA');
    eA_compsmooth{i}=eAsmooth;
    Ncells_comp{i}=Ncells;
    CA_comp{i}=ColonyArea;
    eA_comp{i}=eA;
    tmid_comp{i}=tmid;
end

sz=0;

for i=1:length(filelist)
    if length(eA_compsmooth{1,i}) > sz
        sz = length(eA_compsmooth{1,i});
    end
end

for i=1:length(filelist)
     if length(tmid_comp{1,i}) == sz
        tmid_comp1=tmid_comp(i);
     end
end

tmid_comp1=cell2mat(tmid_comp1);

%growth rate
eA_matrix=cell2mat(eA_compsmooth(1));
eA_matrix2=cell2mat(eA_compsmooth(2));
eA_matrix3=cell2mat(eA_compsmooth(3));
eA_matrix4=cell2mat(eA_compsmooth(4));

%eA_tot=zeros(length(filelist),max+1)
eA_matrix(end:sz)=NaN;
eA_matrix2(end:sz)=NaN;
eA_matrix3(end:sz)=NaN;
eA_matrix4(end:sz)=NaN;

eA_tot=[eA_matrix; eA_matrix2; eA_matrix3; eA_matrix4];

eA_totmean=nanmean(eA_tot,1);
 
colors = [0.85 0.85 0.85;0.75 0.75 0.75; 0.7 0.7 0.7;0.4 0.4 0.4];

figure, hold on
for i=1:length(filelist)
    plot(tmid_comp1/3600,eA_tot(i,:)*3600,'Color',colors(i, :))
end
plot(tmid_comp1/3600,eA_totmean*3600,'--r','LineWidth', 4)
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
legend(samples)
fig2pretty
cd(folder)
savefig([basename 'growthRateAverage.fig'])

for i=1:length(filelist)
     if length(Ncells_comp{1,i}) > sz
        Ncells_comp{1,i}=Ncells_comp{1,i}(1, 1:sz);
    end
end

%convert the colony data to a matrix
Ncells_matrix=cell2mat(Ncells_comp(1));
Ncells_matrix2=cell2mat(Ncells_comp(2));
Ncells_matrix3=cell2mat(Ncells_comp(3));
Ncells_matrix4=cell2mat(Ncells_comp(4));


%assure that matrices are the same length
Ncells_matrix(end:sz)=NaN;
Ncells_matrix2(end:sz)=NaN;
Ncells_matrix3(end:sz)=NaN;
Ncells_matrix4(end:sz)=NaN;

Ncells_tot=[Ncells_matrix; Ncells_matrix2; Ncells_matrix3; Ncells_matrix4];

Ncells_totmean=nanmean(Ncells_tot,1);

figure, hold on
for i=1:length(filelist)
    plot(tmid_comp1/3600,Ncells_tot(i,:)/Ncells_tot(i,1),'Color',colors(i, :))
end
plot(tmid_comp1/3600,Ncells_totmean/Ncells_totmean(1,1),'--r','LineWidth',3)
xlabel('Time (h)')
ylabel('Colony Size (# of cells)')
legend(samples)
fig2pretty
savefig([basename 'colonySizeAverage.fig'])

save([basename '_FCTCompile'])