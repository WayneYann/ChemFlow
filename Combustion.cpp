#include "Combustion.h"

Combustion::Combustion(const ChemThermo& chemThermo)
    : chemThermo_(chemThermo),
      thermo_(chemThermo.gas()),
      nsp_(chemThermo.nsp()),
      maxSteps_(100),
      react_(),
      ode_()
{
    react_.insert(thermo_);
    ode_.addReactor(react_);
}

void Combustion::solve(const double& deltaT, const Eigen::VectorXd& T,
                       const std::vector<Eigen::VectorXd>& Y,
                       const double& p0, std::vector<Eigen::VectorXd>& wdot,
                       Eigen::VectorXd& qdot)
{
    qdot.setZero();
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
            wdot[k](j) = (y[k+2] - y0[k+2])/deltaT;
            qdot(j) -= wdot[k](j)*chemThermo_.Hc(k);
        }
    }
}

void Combustion::solve(double* y, const double& p, const double& xEnd)
{
    const double xStart = 0.0;
    thermo_.setState_TP(y[1], p);
    react_.syncState();
    // sub-step
    double dx = (xEnd - xStart)/maxSteps_;
    ode_.setInitialTime(xStart);
    ode_.reinitialize();
    ode_.setMaxTimeStep(dx);
    ode_.updateState(y);

    for (int nStep=1; nStep<maxSteps_+1; nStep++) {
        // Integrate up to x+dx*nStep
        ode_.advance(xStart + dx*nStep);
    }
    ode_.getState(y);
}

