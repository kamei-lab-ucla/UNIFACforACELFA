%% Tie line estimation: estimated tie lines assuming weight ratio = volume ratio
clear all, clc

%Mechuk parameters for EOPO-SC
a = 71.8490104;
b = -0.461965134;
c = 0.000488689386;

%Initial range of Merchuck function
x_range = linspace(0.001,60,3000);
EOPO_bifunc = @(x)a*exp(b*x.^(.5)-c*x.^3);
EOPO_bi = EOPO_bifunc(x_range);

%Initializing variables
tieline = zeros(1,(length(x_range)));
TLL = zeros(1,1);
fslope = zeros(1,1);
yintercept = zeros(1,1);

xin1 = zeros(1,1);
xin2 = zeros(1,1);
yin1 = zeros(1,1);
yin2 = zeros(1,1);

figure(1) %defines figure 1

%EOPO-SC volume ratio + composition
VR = 3;
p = [10 14.203];

%Initialize range of possible slopes
m_range = linspace(-10,-1/1000,10000);
intercept = p(2)-(p(1)*m_range);
weightratios = zeros(1,10000);
bottomlength = zeros(1,10000);
toplength = zeros(1,10000);
ints = zeros(10000,4);

%Iterate through all possible slopes. For each slope, calculate
%intersection points between tie line and binodal curve. Find ratio of
%lengths between tie line segments. 
for k=1:1:10000
   
EOPO_tiefunc = @(X)m_range(k)*X + intercept(k);
EOPO_tie = EOPO_tiefunc(x_range);

sub = EOPO_bi - EOPO_tie;

neg_slot = find(sub<0);
neg_val = x_range(neg_slot);
xint1 = neg_val(1);
xint2 = neg_val(end);
yint1 = EOPO_tiefunc(xint1);
yint2 = EOPO_tiefunc(xint2);
ints(k,:) = [xint1 xint2 yint1 yint2];

bottomlength(k) = sqrt((yint1-p(2))^2+(p(1)-xint1)^2);
toplength(k) = sqrt((p(2)-yint2)^2+(xint2-p(1))^2);
weightratios(k) = toplength(k)/bottomlength(k);
end

%Find index of tie line slope that corresponds to a ratio of tie line
%segments which is very close to volume ratio
ratiopos = find(abs(weightratios-VR) < 1e-2);
fslope(1) = m_range(ratiopos(1));
yintercept(1) = intercept(ratiopos(1));
xin1(1) = ints(ratiopos(1),1);
xin2(1) = ints(ratiopos(1),2);
yin1(1) = ints(ratiopos(1),3);
yin2(1) = ints(ratiopos(1),4);

%Graph tie line and curve
tieline(1,:) =fslope(1)*(x_range) + yintercept(1);

TLL(1) = bottomlength(ratiopos(1))+toplength(ratiopos(1));

plot(p(1),p(2),'bo')
hold on

plot(x_range,EOPO_bi)
hold on
plot(x_range,tieline)

hold on %plots the phase diagram datapoints

%Print Equation of Tie Line
fprintf("Slope of Tie Line: %d\n", fslope(1));
fprintf("Y-Int of Tie Line: %d\n", yintercept(1));

%Print Intersection Points
fprintf("Polymer Content in Polymer-Rich Phase: %d\n", xin1(1));
fprintf("Salt Content in Polymer-Rich Phase: %d\n", yin1(1));
fprintf("Polymer Content in Salt-Rich Phase: %d\n", xin2(1));
fprintf("Salt Content in Salt-Rich Phase: %d\n", xin2(1));
