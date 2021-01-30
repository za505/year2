%ImageAlign.m
%Registers image stacks using a common reference stack.  
%Calls on dftregistration.m
%
%INSTRUCTIONS FOR USE:
%Save each image stack you want to register in separate directories.  Also
%save the reference image stack (which contains features from which it is 
%possible to track drift) in its own directory.  The program will create
%new directories for each of the image stacks to be registered in the same
%parent directory.  The reference stacks and the stacks to be registered
%should have the same number of images.
%
%INPUT
%basename: choose a name for the directories in which the aligned images
%          will be saved.
%dirname:  cell array containing the paths of the directories in which the
%          image stacks to be registered are.
%regname:  path of reference directory.
%
%CALLS ON:
%dftregistration.m
%imtranslate.m

clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='08052020_hyposhock';
dirname={['/Users/zarina/Downloads/NYU/Lab_2020_Summer/08052020_hyposhock/' basename '_phase'];
      ['/Users/zarina/Downloads/NYU/Lab_2020_Summer/08052020_hyposhock/' basename '_fluor']};
regname=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/08052020_hyposhock/' basename '_crop']; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

workdir=pwd;
ld=length(dirname);

cd(regname)
directory=dir('*.tif');
path(regname,path)
T=length(directory);
tf=T*ld;

imagename=directory(1).name;
RefIm=imread([regname '/' imagename]);

count=0;
for i=1:ld
    
    cd(dirname{i})
    path(dirname{i},path)
    directory2=dir('*.tif');
    cd('../')
    
    mkdir([basename '_' num2str(i) '_a'])
    cd(['./' basename '_' num2str(i) '_a'])
    
    shft=zeros(T,2);
    
    for t=1:T
        count=count+1;
        pd=count/tf*100;
        pds=sprintf(['%4.1f'],pd);
        [pds '%'] %not sure what this part of the code does, or what pds stands for
        
        imagename=directory(t).name;
        TestIm=imread([regname '/' imagename]);
        
        [output NewImFT]=dftregistration(fft2(RefIm),fft2(TestIm),10);
        NewIm=abs(ifft2(NewImFT));
        
        shft(t+1,1)=output(3);
        shft(t+1,2)=output(4);
        
        imagename=directory2(t).name;
        
        I=imread([dirname{i} '/' imagename]);
        [counts,bins]=imhist(I);
        [~,maxpos]=max(counts);
        padcolor=bins(maxpos);
        
        I=imtranslate(I,shft(t+1,:),padcolor);
      
        b=sprintf(['%4.4d'],t);  
        savename=[basename '_a' b '.tif'];
        
        imwrite(I,savename);
    end
end

cd(workdir);