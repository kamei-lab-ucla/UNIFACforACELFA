function[xpolTop, xTop, xpolBottom, xBottom] = compCalc(polTop, saltTop, polBottom, saltBottom, rhoTop, rhoBottom, rhoTot, VR, MW_pol, MW_salt1, MW_salt2, SaltRatio)
MW_H2O = 18.02;
%Salt Ratio should describe the weight ratio of salt1:salt2 in the salt
%stock solution

%Calculating masses of top and bottom phases from densities,
%volume ratio and total mass
mTot = 0.5; %experimental mass of ATPSs
VTot = mTot/rhoTot; %total volume
VTop = (VR/(VR+1))*VTot;
VBottom = (1/(VR+1))*VTot;
mTop = rhoTop*VTop;
mBottom = rhoBottom*VBottom;

%Assuming all concentrations are in % values from 0-100
m_polTop = (polTop/100) * mTop;
m_polBottom = (polBottom/100) * mBottom;
m_saltTop = (saltTop/100) * mTop;
m_saltBottom = (saltBottom/100) * mBottom;

n_H2OTop = (mTop - m_polTop - m_saltTop)/MW_H2O;
n_H2OBottom = (mBottom - m_polBottom - m_saltBottom)/MW_H2O;

salt1Frac = SaltRatio/(SaltRatio+1);
salt2Frac = 1/(SaltRatio+1);

m_Salt1Top = m_saltTop*salt1Frac;
m_Salt2Top = m_saltTop*salt2Frac;
m_Salt1Bottom = m_saltBottom*salt1Frac;
m_Salt2Bottom = m_saltBottom*salt2Frac;

n_polTop = m_polTop / MW_pol;
n_polBottom = m_polBottom / MW_pol;
n_Salt1Top = m_Salt1Top / MW_salt1;
n_Salt2Top = m_Salt2Top / MW_salt2;
n_Salt1Bottom = m_Salt1Bottom / MW_salt1;
n_Salt2Bottom = m_Salt2Bottom / MW_salt2;

n_Top = n_polTop + n_Salt1Top + n_Salt2Top + n_H2OTop;
n_Bottom = n_polBottom + n_Salt1Bottom + n_Salt2Bottom + n_H2OBottom;
nTotal = n_Top + n_Bottom;

xTop = n_Top / nTotal;
xBottom = n_Bottom / nTotal;

xpolTop = n_polTop / n_Top;
xpolBottom = n_polBottom / n_Bottom;
end
