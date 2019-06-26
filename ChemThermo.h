#ifndef CTF_CHEMTHERMO_H_
#define CTF_CHEMTHERMO_H_

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "cantera/thermo.h"
#include "cantera/IdealGasMix.h"
#include "cantera/thermo/IdealGasPhase.h"

const double R = 8314.47;  // [J/kmol K]
const double Tstd = 298.15;  // [K]
const double pstd = 101325.0;  // [Pa]

// Wrapper class for Cantera IdealGasMix with Sutherland transport
class ChemThermo
{
public:
    ChemThermo(const std::string& inputFile, const double& p0);
    ChemThermo() = default;

    ChemThermo& operator=(const ChemThermo&) = delete;

    Cantera::IdealGasMix gas() const {
        return gas_;
    }

    int nsp() const {
        return nsp_;
    }

    // Molecular weight [kg/kmol]
    double W(const int& k) const;

    // Return thermodynamic enthalpy of specie k [J/kmol]
    double ha(const double& p, const double& T, const int& k);

    // Return heat capacity of specie k [J/kmol K]
    double cp(const double& p, const double& T, const int& k);

    int speciesIndex(const std::string& name) const;

    std::string speciesName(const int& k) const;

    void massFractions(const std::vector<Eigen::VectorXd>& Y, double* y, const int& j) const;

    double calcHs(const double& T, const double* y);

    void calcT(Eigen::VectorXd& T, const std::vector<Eigen::VectorXd>& Y,
               const Eigen::VectorXd& hs);

    double sutherland(const double& T) const;

    void updateThermo(const Eigen::VectorXd& T,
                      const std::vector<Eigen::VectorXd>& Y,
                      const double Le, Eigen::VectorXd& rho,
                      Eigen::VectorXd& mu, Eigen::VectorXd& kappa,
                      Eigen::VectorXd& alpha, Eigen::VectorXd& D);


private:
    Cantera::IdealGasMix gas_;

    int nsp_;

    double p0_;

};
#endif  // CTF_CHEMTHERMO_H_
