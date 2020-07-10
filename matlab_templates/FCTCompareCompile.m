%FCTCompile
%
%Rico Rojas 
%6/21/19

%Compiles data from several positions tracked with FluoColonyTrack35

clear, close all

filelist={'110819';'ClpX_downshift';'101719';'102319'};

for i=1:length(filelist)
    load([filelist{i},'_FCTCompile'],'eA_comp','eA_av','tmid')
    eA_comp_cc{i}=eA_comp;
    eA_av_cc{i}=eA_av;
    t_cc{i}=tmid;
end

lfl=length(filelist);
cmap=colormap;
figure, hold on
for i=1:lfl
    for j=1:length(eA_comp_cc{i})
        plot(t_cc{i}/3600,eA_comp_cc{i}{j}*3600,'Color',[0.75 0.75 0.75],'HandleVisibility','off')
    end
    plot(t_cc{i}/3600,eA_av_cc{i},'LineWidth',3,'Color',cmap(round(i*40/lfl),:))
end
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
title('LB -> LB/2')
fig2pretty
legend({'WT';'ClpX';'ClpP';'ClpP 2'})
