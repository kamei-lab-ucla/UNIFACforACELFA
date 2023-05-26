function [dlngammaDAB] = dlngammadDAB(x, T, zVal, r, q, R, Q, nComp, nGroups, compNArray, X, Theta, aij, PsiWaterK, PsiKWater)
%This function returns the partial derivative of the ln of the activity 
%coefficient of DAB with respect to the mole fraction of DAB

z = zVal;
m=0;%index of first row in compNArray which corresponds to a group in DAB
for i = 1:nGroups
    if compNArray(i,2) > 0
        m = i;
        break
    end
end

Psi = exp(-aij/T);

%Calculate mole fraction of each group by multiplying nu values of groups
%with mole fractions of components
calc_1 = compNArray*x'; 
sum_nu_xj = sum(calc_1);
sum_nu_xj_groupsi = sum(calc_1(m:m+2));
sum_nu_xj_groupsnoti = sum_nu_xj - sum_nu_xj_groupsi;
Xmix = calc_1/sum(calc_1)'; %normalize

%Calculate surface area fraction of each group by multiplying Q of each
%group by mole fraction of group, and dividing by total surface area
calc_2 = Q'.*Xmix;
ThetaMix = calc_2/sum(calc_2)'; %normalize

dX = zeros(1,nGroups);
for i = 1:nGroups
    if(i>=m && i<m+3)
        dX(i) = sum_nu_xj_groupsnoti/x(2);
    else 
        dX(i) = -1*sum_nu_xj_groupsi;
    end 
end
dX = (calc_1' .* dX)/(sum_nu_xj*sum_nu_xj);

sum_Q_X = Q*Xmix;
sum_Q_dX = Q*dX';

dTheta = ((sum_Q_X*Q'.*dX')-(sum_Q_dX*Q'.*Xmix))/(sum_Q_X*sum_Q_X);

%Calculates sum over groups n of surface area fraction of group n times 
%interaction parameter of n and another group m, where each row 
%corresponds to group m
%Take transpose, since Psi by itself presents parameters as Psi_mn, when we
%want Psi_nm
sigmaThetanPsinm = Psi'*ThetaMix;
sigmadThetanPsinm = Psi'*dTheta;

%intermediate calc
calc_3 = ((sigmaThetanPsinm.*dTheta)-(ThetaMix.*sigmadThetanPsinm))./(sigmaThetanPsinm.*sigmaThetanPsinm);
dsigmaFraction = Psi*calc_3;

dlnGammaR = Q'.*(sigmadThetanPsinm./sigmaThetanPsinm - dsigmaFraction);

dlngammaR = compNArray'*dlnGammaR;

%Calculate the configuration component of dgamma based on the formula
theta = x.*q/(x*q');
phi = x.*r/(x*r');
%lngammaC = log(phi./x) + (1 - phi./x) - (z/2) * q .* (log(phi./theta)+ ( 1 - phi./theta));
L = (z/2) * (r - q) - (r - 1);
% dlngammaC = -(phi./x) + (z/2) * q .* ((phi./x)-(theta./x)) + ((phi.*phi)./x).*L - (phi./x).*L;
dlngammaC = -(phi./x) + (phi./x).*(phi./x) - (z/2) * q .* ((theta./x)-2*(phi./x)+(phi.*phi)./(x.*theta));

dlngamma = dlngammaR' + dlngammaC;

dlngammaDAB = dlngamma(2);

end