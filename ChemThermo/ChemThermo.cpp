#include "ChemThermo.h"

using namespace Foam;

ChemThermo::ChemThermo(fvMesh& mesh, Time& runTime, const double& p0)
    : pThermo_(rhoReactionThermo::New(mesh)),
      thermo_(pThermo_()),
      pChemistry_(BasicChemistryModel<rhoReactionThermo>::New(thermo_)),
      rho_
      (
          IOobject
          (
              "rho",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
          ),
          thermo_.rho()
      ),
      mu_
      (
          IOobject
          (
              "mu",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
          ),
          thermo_.mu()
      ),
      kappa_
      (
          IOobject
          (
              "kappa",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
          ),
          thermo_.kappa()
      ),
      alphahe_
      (
          IOobject
          (
              "alphahe",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
          ),
          thermo_.alphahe()
      ),
      U_
      (
          IOobject
          (
              "U",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh,
          dimensionedVector("zero", dimVelocity, Zero)
      ),
      chemistry_(pChemistry_()),
      composition_(thermo_.composition()),
      Y_(composition_.Y()),
      p0_(p0)
{
    nsp_ = Y_.size();
}

int ChemThermo::speciesIndex(const std::string& name) const
{
    return composition_.species()[name];
}

std::string ChemThermo::speciesName(const int& k) const
{
    return composition_.species()[k];
}

double ChemThermo::calcHs(const double& T, const double* y) const
{
    double hs = 0.0;
    for (int k=0; k<nsp_; k++) {
        hs += y[k]*composition_.Hs(k, p0_, T);
    }
    return hs;
}

void ChemThermo::calcT(Eigen::VectorXd& T,
                       const std::vector<Eigen::VectorXd>& Y,
                       const Eigen::VectorXd& hs)
{
    for (int j=0; j<T.size(); j++) {
        double y[nsp_];
        massFractions(Y, y, j);
        setY(y);
        thermo_.p() = dimensionedScalar("p", dimPressure, p0_);
        thermo_.he() = dimensionedScalar("h", dimEnergy/dimMass, hs(j));
        thermo_.correct();
        T(j) = thermo_.T()[0];
    }

}

void ChemThermo::massFractions(const std::vector<Eigen::VectorXd>& Y,
                               double* y, const int& j) const
{
    for (int k=0; k<nsp_; k++) {
        y[k] = Y[k](j);
    }
}

void ChemThermo::updateThermo(const Eigen::VectorXd& hs,
                              const std::vector<Eigen::VectorXd>& Y,
                              const double Le, Eigen::VectorXd& rho,
                              Eigen::VectorXd& mu, Eigen::VectorXd& kappa,
                              Eigen::VectorXd& alpha, Eigen::VectorXd& D)
{
    for (int j=0; j<hs.size(); j++) {
        double y[nsp_];
        massFractions(Y, y, j);
        setY(y);
        thermo_.p() = dimensionedScalar("p", dimPressure, p0_);
        thermo_.he() = dimensionedScalar("h", dimEnergy/dimMass, hs(j));
        thermo_.correct();
        syncState();

        rho(j) = rho_[0];
        mu(j) = mu_[0];
        kappa(j) = kappa_[0];
        alpha(j) = alphahe_[0] / rho(j);
        D(j) = alpha(j) / Le;
    }
}

double ChemThermo::solve(const double& deltaT, const Eigen::VectorXd& hs,
                         const std::vector<Eigen::VectorXd>& Y,
                         std::vector<Eigen::VectorXd>& wdot, Eigen::VectorXd& qdot)
{
    double tc = 1.0;
    for (int j=0; j<hs.size(); j++) {
        double y[nsp_];
        massFractions(Y, y, j);
        setY(y);
        thermo_.p() = dimensionedScalar("p", dimPressure, p0_);
        thermo_.he() = dimensionedScalar("h", dimEnergy/dimMass, hs(j));
        thermo_.correct();

        tc = min(tc, chemistry_.solve(deltaT));
        // qdot(j) = chemistry_.Qdot()()[0];
        for (int k=0; k<nsp_; k++) {
            wdot[k](j) = chemistry_.RR(k)[0];
        }
    }
    this->filter(wdot);
    // Compute qdot
    for (int j=0; j<hs.size(); j++) {
        qdot(j) = 0.0;
        for (int k=0; k<nsp_; k++) {
            qdot(j) -= composition_.Hc(k)*wdot[k](j);
        }
    }
    return max(tc, small);
}

void ChemThermo::setY(const double* y)
{
    for (int k=0; k<nsp_; k++) {
        Y_[k] = y[k];
    }
}

void ChemThermo::syncState()
{
    rho_ = thermo_.rho();
    mu_ = thermo_.mu();
    kappa_ = thermo_.kappa();
    alphahe_ = thermo_.alphahe();
}

void ChemThermo::filter(std::vector<Eigen::VectorXd>& wdot) const
{
    std::vector<Eigen::VectorXd> wdotOrig(nsp_);
    for (int k=0; k<nsp_; k++) {
        wdotOrig[k].resize(wdot[k].size());
        for (int j=0; j<wdotOrig[k].size(); j++) {
            wdotOrig[k](j) = wdot[k](j);
        }
    }
    // 5-point averaging
    for (int k=0; k<nsp_; k++) {
        for (int j=2; j<wdot[k].size()-2; j++) {
            wdot[k][j] = 0.08*wdotOrig[k](j-2) + 0.17*wdotOrig[k](j-1)
                        + 0.5*wdotOrig[k](j)
                        + 0.17*wdotOrig[k](j+1) + 0.08*wdotOrig[k](j+2);
        }
    }
}
