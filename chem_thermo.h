#ifndef CTF_ChemThermo_H_
#define CTF_ChemThermo_H_

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "cantera/thermo.h"
#include "cantera/IdealGasMix.h"
#include "cantera/thermo/IdealGasPhase.h"

const double R = 8314.47;  // [J/kmol K]
const double Tstd = 298.15;  // [K]

// Wrapper class of Cantera IdealGasMix with Sutherland transport
class ChemThermo
{
public:
    ChemThermo(const std::string& inputFile, const double& p0);
    ChemThermo() = default;

    ChemThermo(const ChemThermo&) = delete;
    ChemThermo& operator=(const ChemThermo&) = delete;

    int nsp() const {
        return nsp_;
    }

    int speciesIndex(const std::string& name) const;

    std::string speciesName(const int& k) const;

    void massFractions(const std::vector<Eigen::VectorXd>& Y, double* y, const int& j) const;

    double calcHs(const double& T, const double* y);

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
#endif  // CTF_ChemThermo_H_
