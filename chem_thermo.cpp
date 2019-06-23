#include <cmath>
#include "chem_thermo.h"

ChemThermo::ChemThermo(const std::string& inputFile, const double& p0)
    : gas_(inputFile,"gas"),
      p0_(p0)
{
    nsp_ = gas_.nSpecies();
}

int ChemThermo::speciesIndex(const std::string& name) const
{
    return gas_.speciesIndex(name);
}

std::string ChemThermo::speciesName(const int& k) const
{
    return gas_.speciesName(k);
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

