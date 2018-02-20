%executable.m

f1=@McCT;


F=1000; %Generic flowrate
x1f=0.5; %Feed composition
x1d=0.96; %Distillate mole fraction of volatile component
%Note: Increasing x1d above 0.95 will give performance issues; %will be
%rectified in future version


x1b=0.025; %Bottoms mole fraction of volatile component
Rratio=2; %R/Rmin
q=1; %Feed stream quality (q=1 indicates saturated liquid)


[Ntray,Rmin,L,V,Lbar,Vbar,B,D]=McCT(F,x1f,x1d,x1b,Rratio,q);

%Output:
%Ntray = number of trays
%Rmin = minimum reflux ratio
%L and V are the enriching section flows
%Lbar and Vbar are the stripping section flows
%B and D are the bottoms and distillate flows