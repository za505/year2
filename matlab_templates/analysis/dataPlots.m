%dataPlots.m Author: Zarina Akbary Date: 30 March 2021 Purpose: to gather
%all the data in one place and plot them to make sound conclusions.

clear, close all
%%%%%%%%%%%%%

%% User Input
%first, let's start with the easier one, 647
maindir=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/data_analysis'];
filename={['03172021_Exp1'] ['03252021_Exp2'] ['03292021_Exp2']};
channel=['_647'];

%we'd never included index as a variable in the previous scripts in our
%workflow, so we have to hard code it. for the Mg2+ gradient experiments, it is this ['first
%frame with dye' 'last frame of first perfusion' 'last frame of second
%perfusion' 'last frame of third perfusion' 'last frame of fourth
%perfusion']
experiment(1).index=[19 30 42 54 61];
experiment(2).index=[20 30 40 54 66];
experiment(3).index=[21 29 42 53 67];

%we only have 4 switches in each experiment
nSwitches=4;

%and this is our range for the Mg2+ gradient
mgRange=[0 6.66 12 12.33 15 20];

%%%%%%%%%%%%%%%%
%% 647

for f=1:length(filename)
    
    %we are going to go through each experiment and find the BTfluoAVG
    %data we'd stored
    file=char(filename(f));
    cd(maindir)
    experiment(f).BTfluoAVG=load([file '_BTfluoAVG' channel '.mat']);
    experiment(f).colonies=length(experiment(f).BTfluoAVG.basename);
    experiment(f).basename=experiment(f).BTfluoAVG.basename;
    experiment(f).mgRange=experiment(f).BTfluoAVG.mgRange;
  
    index=experiment(f).index;
    experiment(f).ratio{1}=mean(experiment(f).BTfluoAVG.ratio(:, index(1):index(2)), 2, 'omitnan');
    experiment(f).ratio{2}=mean(experiment(f).BTfluoAVG.ratio(:, index(2)+1:index(3)), 2, 'omitnan');
    experiment(f).ratio{3}=mean(experiment(f).BTfluoAVG.ratio(:, index(3)+1:index(4)), 2, 'omitnan');
    experiment(f).ratio{4}=mean(experiment(f).BTfluoAVG.ratio(:, index(4)+1:index(5)), 2, 'omitnan');

    temp1_lcell=[];
    temp2_lcell=[];
    temp3_lcell=[];
    temp4_lcell=[];
    
    %now, let's go through each colony and get out BTfluo and BTphase data
    for g=1:experiment(f).colonies
        
        base = char(experiment(f).BTfluoAVG.basename(g));
        %experiment(f).BTfluo{g}=load([base channel '_BTfluo.mat']);
        experiment(f).BTphase{g}=load([base '_BTphase.mat']);
        
        temp1 = experiment(f).BTphase{1, g}.lcell(:, index(1):index(2));
        temp2 = experiment(f).BTphase{1, g}.lcell(:, index(2)+1:index(3));
        temp3 = experiment(f).BTphase{1, g}.lcell(:, index(3)+1:index(4));
        temp4 = experiment(f).BTphase{1, g}.lcell(:, index(4)+1:index(5));
        
        temp1_lcell=[temp1_lcell; temp1];
        temp2_lcell=[temp2_lcell; temp2];
        temp3_lcell=[temp3_lcell; temp3];
        temp4_lcell=[temp4_lcell; temp4];
      
    end
    
    experiment(f).lcell{1}=mean(temp1_lcell, 2, 'omitnan');
    experiment(f).lcell{2}=mean(temp2_lcell, 2, 'omitnan');
    experiment(f).lcell{3}=mean(temp3_lcell, 2, 'omitnan');
    experiment(f).lcell{4}=mean(temp4_lcell, 2, 'omitnan');
    
%     %let's get the range of Mg2+ concentrations used
%     experiment(f).mgRange=cell2mat(struct2cell(load([filelist{f}.name],'mgRange')));
%     %and the number of time points
%     experiment(f).T=cell2mat(struct2cell(load([filelist{f}.name],'T')));
%     %what were the temporal averages of the ratios during each switch?
%     experiment(f).avgMg=cell2mat(struct2cell(load([filelist{f}.name],'avgMg')));
%     %what are the ratios themselves?
%     experiment(f).ratio=cell2mat(struct2cell(load([filelist{f}.name],'ratio')));
%     %what about the BTphase data? We didn't compile them in BTfluoAVG,
%     but %what about here?
%     experiment(f).BTphase=struct2cell(load([filelist{f}.name],'filelist'));
%     experiment(f).BTphase=experiment(f).BTphase{1,1};
%     
%     %alright, let's go through and pull some BTphase data out
%     temp_length=length(experiment(f).BTphase); for g=1:temp_length
%         temp_struct=experiment(f).BTphase{1,g};
%         temp_path=load([temp_struct.folder '/' temp_struct.name],
%         'curdir')
%         
%     end
    
end

avgRatio={};
avgLength={};

for h=1:length(mgRange)
    
    mgRange(h)
    tempRatio=[];
    tempLength=[];
    
    for f=1:length(filename)
        if any(experiment(f).mgRange==mgRange(h))
            idx=find(experiment(f).mgRange==mgRange(h));
            tempRatio=[tempRatio; experiment(f).ratio{idx}];
            tempLength=[tempLength; experiment(f).lcell{idx}];  
        end 
    end
    
    avgRatio{1,h}=tempRatio;
    avgLength{1,h}=tempLength;
end
