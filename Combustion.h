#ifndef CTF_COMBUSTION_H_
#define CTF_COMBUSTION_H_

#include <cmath>
#include <vector>
#include <Eigen/Dense>

#include "ChemThermo.h"
#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/ReactorNet.h"

class Combustion
{
public:
    Combustion() = default;
    Combustion(const ChemThermo& chemThermo);

    Combustion(Combustion&) = delete;
    Combustion& operator=(const Combustion&) = delete;

    // Global solve -- returns integrated reaction rates and heat release rate
    void solve(const double& deltaT, const Eigen::VectorXd& T,
               const std::vector<Eigen::VectorXd>& Y,
               const double& p, std::vector<Eigen::VectorXd>& wdot,
               Eigen::VectorXd& qdot);

    // Solve at each grid point
    void solve(double* y, const double& p,
               const double& xEnd);


private:
    ChemThermo chemThermo_;

    Cantera::IdealGasMix thermo_;
    Cantera::IdealGasConstPressureReactor react_;
    Cantera::ReactorNet ode_;

    int nsp_;

    int maxSteps_;
};


#endif  // CTF_COMBUSTION_H_
