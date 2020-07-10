%poissonDist.m
%
%Zarina Akbary, 6/5/2020
%
%Purpose: To simulate a Poisson distribution. In-class NSC assignment
%get the number of inputs that occur in some time=T
%
%empty m-file
clear, close all

numrep=10^4;
nvals=[];

for T=2:2:20
for i=1:numrep
T=10;
n=0;
k=1;
t=0;
N=1000;
dt=T/N;

while t<T
    x=rand;
    if x<k*dt
        n=n+1;
    end
    t=t+dt;
end
nvals(end+1)=n;
end
nav=mean(nvals);
nstd=std(nvals);
cv=nstd/nav;
cv_vals(end+1)=cv;
end

T=2:2:20;
    
figure
hist(nvals)
plot(T, cv_vals)
plot(T,arrayfun(@(T)1/sqrt(T),T)) %fix