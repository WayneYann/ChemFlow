#include <fstream>
#include "Combustion.h"

Combustion::Combustion(const ChemThermo& chemThermo)
    : thermo_(chemThermo.gas()),
      nsp_(chemThermo.nsp()),
      reactionRate_(chemThermo.nsp()),
      maxSteps_(100),
      react_(),
      ode_()
{
    react_.setThermoMgr(thermo_);
    ode_.addReactor(react_);
}

void Combustion::solve(const double& deltaT, const Eigen::VectorXd& T,
                       const std::vector<Eigen::VectorXd>& Y,
                       const double& p0)
{
    for (int k=0; k<nsp_; k++) {
        reactionRate_[k].resize(T.size());
    }

    for (int j=0; j<T.size(); j++) {
        double y[nsp_+2], y0[nsp_+2];
        y[0] = 1.0;
        y0[0] = 1.0;
        y[1] = T(j);
        y0[1] = T(j);
        for (int k=0; k<nsp_; k++) {
            y[k+2] = Y[k](j);
            y0[k+2] = Y[k](j);
        }

        this->solve(y, p0, deltaT);

        for (int k=0; k<nsp_; k++) {
            reactionRate_[k](j) = (y[k+2] - y0[k+2])/deltaT;
        }
    }
    Eigen::VectorXd::Index loc;
    std::cout << "Max omega of " << thermo_.speciesName(11) << "   "
              << reactionRate_[11].maxCoeff(&loc) << " @ position "
              << loc << std::endl;
}

void Combustion::solve(double* y, const double& p, const double& xEnd)
{
    const double xStart = 0.0;
    thermo_.setState_TP(y[1], p);
    // sub-step
    double dx = (xEnd - xStart)/maxSteps_;
    ode_.reinitialize();
    ode_.setInitialTime(xStart);
    ode_.setMaxTimeStep(dx);
    ode_.updateState(y);

    for (int nStep=1; nStep<maxSteps_+1; nStep++) {
        // Integrate up to x+dx*nStep
        ode_.advance(xStart + dx*nStep);
    }
    ode_.getState(y);
}

