#include <stdexcept>
#include "Combustion.h"

Combustion::Combustion(const ChemThermo& thermo)
    : thermo_(thermo),
      nsp_(thermo.nsp()),
      c_(thermo.nsp(), 0.0),
      c0_(thermo.nsp(), 0.0),
      reactionRate_(thermo.nsp()),
      maxSteps_(10000)
{
}

double Combustion::solve(const double& deltaT, const Eigen::VectorXd& T,
                         const std::vector<Eigen::VectorXd>& Y,
                         const Eigen::VectorXd& rho, const double& p0)
{
    double deltaTMin = 1e3;
    double p = p0;
    for (int k=0; k<nsp_; k++) {
        reactionRate_[k].resize(T.size());
    }
    deltaTChem_.resize(T.size());

    for (int j=0; j<T.size(); j++) {
        double Tj = T(j);
        const double rhoj = rho(j);
        // Calculate molar concentrations [kmol/m3] and store initial values
        for (int k=0; j<nsp_; k++) {
            c_[k] = rhoj*Y[k](j)/thermo_.W(k);
            c0_[k] = c_[k];
        }

        double timeLeft = deltaT;
        while (timeLeft > 1e-10) {
            double dt = timeLeft;
            this->solve(c_, Tj, p, dt, deltaTChem_(j));
            timeLeft -= dt;
        }
        deltaTMin = std::min(deltaTMin, deltaTChem_(j));

        for (int k=0; k<nsp_; k++) {
            reactionRate_[k](j) = (c_[k] - c0_[k])*thermo_.W(k)/deltaT;
        }
    }
    return deltaTMin;
}

void Combustion::solve(std::vector<double>& c, double& T,
                       double& p, const double& xEnd, double& dxTry) const
{
    std::vector<double> y(nsp_+2);
    for (int k=0; k<nsp_; k++) {
        y[k] = c[k];
    }
    y[nsp_] = T;
    y[nsp_+1] = p;

    // Start from x = 0.0
    const double xStart = 0.0;

    stepState step(dxTry);
    double x = xStart;

    int nStep=0;
    for (; nStep<maxSteps_; nStep++) {
        double dxTry0 = step.dxTry;
        step.reject = false;

        // Check if this is a truncated step
        if ((x + step.dxTry - xEnd)*(x + step.dxTry - xStart) > 0) {
            step.last = true;
            step.dxTry = xEnd - x;
        }

        // Integrate as far as possible up to step.dxTry
        solve(x, y, step);

        // Check if reached xEnd
        if ((x - xEnd)*(xEnd - xStart) >= 0) {
            if (nStep>0 && step.last) {
                step.dxTry = dxTry0;
            }
            dxTry = step.dxTry;
            break;
        }
        step.first = false;
        if (step.reject) {
            step.prevReject = true;
        }
    }
    if (nStep >= maxSteps_) throw std::runtime_error("Maximum integration steps reached");

    for (int k=0; k<nsp_; k++) {
        c[k] = std::max(0.0, y[k]);
    }
    T = y[nsp_];
    p = y[nsp_+1];
}

void Combustion::solve(double& x, std::vector<double>& y, stepState& step) const
{
}


