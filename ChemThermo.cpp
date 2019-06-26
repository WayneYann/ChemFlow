#include <cmath>
#include <stdexcept>
#include "ChemThermo.h"

ChemThermo::ChemThermo(const std::string& inputFile, const double& p0)
    : gas_(inputFile,"gas"),
      p0_(p0)
{
    nsp_ = gas_.nSpecies();
}

double ChemThermo::W(const int& k) const
{
    return gas_.molecularWeight(k);
}

double ChemThermo::ha(const double& p, const double& T, const int& k)
{
    gas_.setState_TP(T, p);
    return R*T*gas_.enthalpy_RT_ref()[k];
}

double ChemThermo::cp(const double& p, const double& T, const int& k)
{
    gas_.setState_TP(T, p);
    return R*gas_.cp_R_ref()[k];
}

int ChemThermo::speciesIndex(const std::string& name) const
{
    return gas_.speciesIndex(name);
}

std::string ChemThermo::speciesName(const int& k) const
{
    return gas_.speciesName(k);
}

double ChemThermo::calcHs(const double& T, const double* y)
{
    // This order also serves the purpose of setting gas_ state(T, y)
    gas_.setState_TPY(Tstd, pstd, y);
    const double h0 = gas_.enthalpy_mass();
    gas_.setState_TPY(T, p0_, y);
    const double ha = gas_.enthalpy_mass();

    return ha - h0;
}

void ChemThermo::calcT(Eigen::VectorXd& T,
                       const std::vector<Eigen::VectorXd>& Y,
                       const Eigen::VectorXd& hs)
{
    for (int j=0; j<T.size(); j++) {
        double y[nsp_];
        massFractions(Y, y, j);

        // Newton iterative method
        double Test = T(j);
        double Tnew = T(j);
        double TnewS = T(j);
        double Ttol = T(j)*1e-04;
        double fx, alpha;
        int iter = 0;
        do
        {
            Test = TnewS;
            fx = calcHs(Test, y) - hs(j);
            Tnew = Test - fx / gas_.cp_mass();
            alpha = 1.0;
            TnewS = Tnew;
            while (std::abs(fx) < std::abs(calcHs(TnewS, y) - hs(j))) {
                alpha /= 2.0;
                TnewS = alpha*Tnew + (1.0-alpha)*Test;
            }
            if (iter++ > 100) {
                throw std::runtime_error("Maximum iterations reached");
            }
        } while (std::abs(Tnew - Test) > Ttol);
        T(j) = Tnew;
    }
}


double ChemThermo::sutherland(const double& T) const
{
    const double As = 1.67212e-6;
    const double Ts = 170.672;
    return As*std::sqrt(T) / (1.0+Ts/T);
}

void ChemThermo::updateThermo(const Eigen::VectorXd& T,
                              const std::vector<Eigen::VectorXd>& Y,
                              const double Le, Eigen::VectorXd& rho,
                              Eigen::VectorXd& mu, Eigen::VectorXd& kappa,
                              Eigen::VectorXd& alpha, Eigen::VectorXd& D)
{
    for (int j=0; j<T.size(); j++) {
        double y[nsp_];
        massFractions(Y, y, j);
        gas_.setState_TPY(T(j), p0_, y);
        const double WM = gas_.meanMolecularWeight();  // [kg/kmol]
        const double Cv = gas_.cv_mass();  // [J/kg K]
        const double Cp = gas_.cp_mass();  // [J/kg K]
        const double Rg = R/WM;  // [J/kg K]

        rho(j) = gas_.density();  // [kg/m3]
        mu(j) = sutherland(T(j));  // [kg/m s]
        kappa(j) =  mu(j)*Cv*(1.32 + 1.77*Rg/Cv);  // [W/m K]
        alpha(j) = kappa(j)/(rho(j)*Cp);  // [m2/s2]
        D(j) = alpha(j)/Le;  // [m2/s2]
    }
}

void ChemThermo::massFractions(const std::vector<Eigen::VectorXd>& Y,
                               double* y, const int& j) const
{
    for (int k=0; k<nsp_; k++) {
        y[k] = Y[k](j);
    }
}

