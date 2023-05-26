function[xpolTop, xTop, xpolBottom, xBottom] = comp_calc(polTop, saltTop, polBottom, saltBottom, MW_pol, MW_salt1, MW_salt2, SaltRatio)
MW_H2O = 18.02;
%Salt Ratio should describe the weight ratio of salt1:salt2 in the salt
%stock solution

%Assuming overall mass of 0.5 g (doesn't really matter)
%Also assuming all concentrations are in % values from 0-100
m_polTop = polTop * 0.005;
m_polBottom = polBottom * 0.005;
m_saltTop = saltTop * 0.005;
m_saltBottom = saltBottom * 0.005;

n_H2OTop = ((100-polTop-saltTop)* 0.005)/MW_H2O;
n_H2OBottom = ((100-polBottom-saltBottom)* 0.005)/MW_H2O;

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

n_Top = n_polTop + n_Salt1Top+ n_Salt2Top + n_H2OTop;
n_Bottom = n_polBottom + n_Salt1Bottom + n_Salt2Bottom + n_H2OBottom;
nTotal = n_Top + n_Bottom;

xTop = n_Top / nTotal;
xBottom = n_Bottom / nTotal;

xpolTop = n_polTop/n_Top;
xpolBottom = n_polBottom/n_Bottom;
end
