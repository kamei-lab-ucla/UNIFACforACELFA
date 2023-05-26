function [gamma] = unifac(x, T, zVal, r, q, R, Q, nComp, nGroups, compNArray, X, Theta, aij, PsiWaterK, PsiKWater)
%File is adapted from Elliott and Lira (https://sourceforge.net/projects/chethermo/)

%This function returns a row vector(1,nComp) which contains the activity
%coefficients of every component within the system

%Note that this function modifies the activity coefficient slightly using
%the delta value
%This is meant to convert the activity coefficient so that it is based on
%an infinite dilution standard state
%Since we are equating mu's in this problem, the standard states will
%cancel, so this does not affect the calculation of the partition
%coefficient at all
%However, the calculations in this file vary slightly from the commonly
%seen form of the UNIFAC equation

%x is a row vector(1,nComp) containing mole fraction of each component
%T is temperature(K) of system
%nComp(int) is the number of components
%nGroups(int) is the total number of groups
%r is a row vector(1,nComp) containing r value for each component
%q is a row vector(1,nComp) containing q value for each component
%R is a row vector(1,nGroups) containing R value for each group
%Q is a row vector(1,nGroups) containing Q value for each group

%compNArray is an array(nGroups,nComp) containing the nu values for each
%group within each component, where each column corresponds to the same
%component within x (rows are groups)

%X is an array(nGroups,nComp) where each column represents the mole
%fractions of each group in a pure solution of just the component that the
%column corresponds to

%Theta is an array(nGroups,nComp) where each column represents the surface
%area fractions of each group within a pure solution of the component that
%column corresponds to

%aij is an array(nGroups,nGroups) which enumerates all the aij values for
%all group interactions with one another

%PsiWaterK is a row vector(1,nGroups) which contains all the interaction
%parameters between water and the other groups, in that order

%PsiKWater is a row vector(1,nGroups) which contains all the interactions
%parameters between the groups in solution and water, in that order

z = zVal;
r_water = 0.92;
q_water = 1.4;

e = 1.60e-19;
Na = 6.6022e23;
kb = 1.381e-23;
e0 = 8.854e-12;
ew = 80;
B = sqrt((2000*Na*e*e)/(e0*ew*kb*T));
A = (e*e*B)/(e0*ew*kb*T);
a0 = 4e-10; %ionic size

%Calculate mole fraction of each group by multiplying nu values of groups
%with mole fractions of components
calc_1 = compNArray*x'; 
Xmix = calc_1/sum(calc_1)'; %normalize

X = [X Xmix]; %add mole fraction of mixed state

%convert from aij to Psi
Psi = exp(-aij/T);
PsiWaterK = exp(-PsiWaterK/T);
PsiKWater = exp(-PsiKWater/T);

%Calculate surface area fraction of each group by multiplying Q of each
%group by mole fraction of group, and dividing by total surface area
calc_2 = Q'.*Xmix;
ThetaMix = calc_2/sum(calc_2)'; %normalize

Theta = [Theta ThetaMix]; %append mixed state

%Calculates sum over groups n of surface area fraction of group n times 
%interaction parameter of n and another group m, where each row 
%corresponds to group m
%Take transpose, since Psi by itself presents parameters as Psi_mn, when we
%want Psi_nm
sigmaThetanPsinm = Psi'*Theta;

%Now we want to take previous sum and use that as the denominator in a new
%formula
%We want to the take the sum over all groups i of the fraction where the
%numerator is surface area fraction of group i times the interaction
%parameter of another group k and i, and the denominator is the previous
%sum
%Here we use PSI by itself, since we want Psi_ki, not Psi_ik
%Comparing to previous sum, index k is the same as index m, 
%and index i is the same as index n
calc_3 = Theta./sigmaThetanPsinm; %intermediate calc
sigmaFraction = Psi*calc_3; %multiplies and takes sum

%Based on the formula, for every group k, we need to multiply Q_k with 1 -
%log of the first sum - the sum of the large fraction
%We use the kron() function to duplicate the column vector Q' so that the
%number of columns matches the number of columns in our other variables
%This is because we are calculating all of these values for solutions of
%pure component for every component and also for the mixed state
lnGammaR = kron(Q',ones(1,nComp+1)).*(1-log(sigmaThetanPsinm) - sigmaFraction);

%Determines same value for just water and k interactions
lnGamma_del = Q.*(1-log(PsiWaterK) - PsiKWater);

% Calculate residual part of gammas, by subtracting the lnGamma of the
% mixed state and the lnGamma of the pure component for each group, and
% then multiplying by the nu value of each group, then summing over all
% groups within that component
% deltaR is the same value for just water and k interactions
lngammaRes= zeros(1,nComp);
deltaR = zeros(1,nComp);
for i = 1:nComp 
    lngammaRes(i) = compNArray(:,i)'*(lnGammaR(:,nComp+1)-lnGammaR(:,i));
    deltaR(i) = compNArray(:,i)'*(lnGamma_del(1,:)'-lnGammaR(:,i));
end

%Calculate the configuration component of gamma based on the formula
theta = x.*q/(x*q');
phi = x.*r/(x*r');
% lngammaC = log(phi./x) + (1 - phi./x) - (z/2) * q .* (log(phi./theta)+ ( 1 - phi./theta));
L = (z/2) * (r - q) - (r - 1);
sum_LX = x*L';
lngammaC = log(phi./x) + (z/2) * q .* log(theta./phi) + L - sum_LX * (phi./x);
L_water = (z/2) * (r_water - q_water) - (r_water - 1);
deltaC = log(r./r_water) + (z/2) * q.* log((q./q_water)./(r./r_water)) + L - L_water*(r./r_water);
delta = -(deltaR + deltaC);

%delta is the conversion factor to an infinite dilution standard state
%Basically removes all interactions with water
gamma = exp(lngammaC + lngammaRes + delta);
%gamma = exp(lngammaC + lngammaRes);

end
