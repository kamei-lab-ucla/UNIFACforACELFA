
% For EOPO-DAB
%This array decomposes each component into its constituent groups
%Each column represents a different component
%In our case, column 1 is EOPO, column 2 is TMB/DAB, column 3 is water
% disp('EOPO-DAB')
% compArray = {
%           'OH'        2       0       0;
%           'CH2'       205     0       0;
%           'CH2O'      255     0       0;
%           'CH'        52      0       0;
%           'CH3'       52      0       0;
%           'AC'        0       2       0;
%           'ACH'       0       6       0;
%           'ACNH2'     0       4       0;
%           'H2O'       0       0       1;
%           };

% For EOPO-TMB
disp('EOPO-TMB')
compArray = {
          'OH'        2       0       0;
          'CH2'       205     0       0;
          'CH2O'      255     0       0;
          'CH'        52      0       0;
          'CH3'       52      0       0;
          'AC'        0       2       0;
          'ACH'       0       4       0;
          'ACCH3'     0       4       0;
          'ACNH2'     0       2       0;
          'H2O'       0       0       1;
          };

[paramsPure] = unifacSetUp(compArray);%Grab parameters of interest within this system

% This is for 14.203% EOPO, 10% SC Updated
disp('14.203% EOPO, 10% SC')
% Old numbers
% polTop = 18.8847;
% saltTop = 6.6631;
% polBottom = 0.1908;
% saltBottom = 19.9873;
% New numbers, relaxed density assumption
polTop = 19.22937;
saltTop = 6.543072;
rhoTop = 1.068007;%density of top phase
polBottom = 0.2341185;
saltBottom = 19.60721;
rhoBottom = 1.150995;%density of bottom phase
rhoTot = 1.089967;%density of overall ATPS
VR = 3;
MW_pol = 12000;
MW_TriCit = 258.06;
MW_CitAcid = 192.194;
SaltRatio = 2.6;%weight ratio of trisodium citrate to citric acid
[xpolTop, xTop, xpolBottom, xBottom] = compCalc(polTop, saltTop, polBottom, saltBottom, rhoTop, rhoBottom, rhoTot, VR, MW_pol, MW_TriCit, MW_CitAcid, SaltRatio);
xTotal = 2.583778197e-5; %overall mole fraction of TMB/DAB
comp_top = [xpolTop 0 0];
comp_bottom = [xpolBottom 0 0];
mole_frac_top = xTop;%mole fraction of top phase in whole ATPS
mole_frac_bottom = xBottom;%mole fraction of bottom phase in whole ATPS


xT = xTotal;%mole fraction of TMB/DAB in the top phase
xB = xTotal;%mole fraction of TMB/DAB in the bottom phase
des_err = 0.0000000001;
i = 1;%counter variable for iterations
j = 1;%counter variable for iterations
T = 298;%temperature
z = 10;%Number of nearest neighbors

%update the vector with our initial guesses for mole fraction of TMB/DAB
comp_top(2) = xT;
comp_bottom(2) = xB;

%Calculate the mole fraction of water, so that all mole fractions sum to 1
comp_top(3) = 1-(comp_top(1) + comp_top(2));
comp_bottom(3) = 1-(comp_bottom(1) + comp_bottom(2));

%Calculate partial derivative of the natural log of the activity
%coefficient of TMB/DAB with respect to the mole fraction of TMB/DAB 
dlngammaDAB_bottom = dlngammadDAB(comp_bottom, T, z, paramsPure{:});
dlngammaDAB_top = dlngammadDAB(comp_top, T, z, paramsPure{:});

%Calculate activity coefficients
gamma_bottom = unifac(comp_bottom, T, z, paramsPure{:});
gamma_top = unifac(comp_top, T, z, paramsPure{:});

%Convert partial derivative of natural log of activity coefficient to
%partial derivative of activity coefficient
dgammaDAB_bottom = gamma_bottom(2)*dlngammaDAB_bottom;
dgammaDAB_top = gamma_top(2)*dlngammaDAB_top;

%Calculate initial values of f1 and f2
f1 = xT*gamma_top(2) - xB*gamma_bottom(2);
f2 = xT*mole_frac_top + xB*mole_frac_bottom - xTotal;

%If f1 and f2 are not equal to 0, use Jacobian and Gauss-Seidel to find 
%step sizes and update the mole fractions until f1 and f2 are equal to 0
while (sqrt(f1*f1+f2*f2)>des_err)
    i = i+1;
    del_xB = 0;
    del_xT = 0;
    
    %Calculate coefficients in Jacobian
    df1_dxB = -1*(gamma_bottom(2) + comp_bottom(2)*dgammaDAB_bottom);
    df1_dxT = gamma_top(2) + comp_top(2)*dgammaDAB_top;

    df2_dxB = mole_frac_bottom;
    df2_dxT = mole_frac_top;
    
    del_xBnew = (-f1 - del_xT*df1_dxT)/df1_dxB;
    del_xTnew = (-f2 - del_xBnew*df2_dxB)/df2_dxT;
    
    approx_err_1 = abs((del_xTnew-del_xT)/del_xTnew);
    approx_err_2 = abs((del_xBnew-del_xB)/del_xBnew);
    approx_err = max([approx_err_1,approx_err_2]);
    
    del_xB = del_xBnew;
    del_xT = del_xTnew;
    
    while(approx_err>0.1)
        j = j + 1;
        del_xBnew = (-f1 - del_xT*df1_dxT)/df1_dxB;
        del_xTnew = (-f2 - del_xBnew*df2_dxB)/df2_dxT;

        approx_err_1 = abs((del_xTnew-del_xT)/del_xTnew);
        approx_err_2 = abs((del_xBnew-del_xB)/del_xBnew);
        approx_err = max([approx_err_1,approx_err_2]);

        del_xB = del_xBnew;
        del_xT = del_xTnew;
    end
    
    xB = xB + del_xB;
    xT = xT + del_xT;
    
    comp_top(2) = xT;
    comp_bottom(2) = xB;

    comp_top(3) = 1-(comp_top(1) + comp_top(2));
    comp_bottom(3) = 1-(comp_bottom(1) + comp_bottom(2));

    dlngammaDAB_bottom = dlngammadDAB(comp_bottom, T, z, paramsPure{:});
    dlngammaDAB_top = dlngammadDAB(comp_top, T, z, paramsPure{:});

    gamma_bottom = unifac(comp_bottom, T, z, paramsPure{:});
    gamma_top = unifac(comp_top, T, z, paramsPure{:});

    dgammaDAB_bottom = gamma_bottom(2)*dlngammaDAB_bottom;
    dgammaDAB_top = gamma_top(2)*dlngammaDAB_top;

    f1 = xT*gamma_top(2) - xB*gamma_bottom(2);
    f2 = xT*mole_frac_top + xB*mole_frac_bottom - xTotal;
end

K = gamma_bottom(2)/gamma_top(2);
fprintf("Predicted K: %d\n", K);
