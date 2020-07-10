%FCTCompile
%
%Rico Rojas 
%6/21/19

%Compiles data from several positions tracked with FluoColonyTrack35

clear, close all

filelist={'060419_P1';'060419_P2'};
filelist={'110819_P1';'110819_P2';'110819_P3';'110819_P4'};,basename='110819';%WT B.s. LB->LB/2
%filelist={'ClpX_downshift';'ClpX_downshift_2';'ClpX_downshift_3'};,basename='ClpX_downshift';%ClpX B.s. LB->LB/2
%filelist={'101719_P1';'101719_P2'};,basename='101719';%ClpP B.s. LB->LB/2
%filelist={'102319_P1_C1';'102319_P1_C2'};,basename='102319';%ClpP B.s. LB->LB/2

smooth=10;

for i=1:length(filelist)
    load([filelist{i},'_FCT'],'Ncells','tmid','time','tscale','ColonyArea','eAsmooth')
    Ncells_comp{i}=Ncells;
    CA_comp{i}=ColonyArea;
    eA_comp{i}=eAsmooth;
    tmid_comp{i}=tmid;
end

Ncells_sum=zeros(size(Ncells));
CA_sum=zeros(size(Ncells));
for i=1:length(filelist)
    CA_sum=CA_sum+CA_comp{i};
end

dA=CA_sum(2:end)-CA_sum(1:end-1);
eA_av=dA./CA_sum(1:end-1)/(tscale)*3600;
eA_av=movingaverage(eA_av,smooth);

figure, hold on
for i=1:length(filelist)
    plot(tmid/3600,eA_comp{i}*3600,'Color',[0.75 0.75 0.75])
end
plot(tmid/3600,eA_av,'-r','LineWidth',3)
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
fig2pretty

save([basename '_FCTCompile'])