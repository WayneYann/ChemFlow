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

        const scalarField rhoC = thermo_.rho();
        const scalarField muC = thermo_.mu();
        const scalarField kappaC = thermo_.kappa();
        const scalarField alphaheC = thermo_.alphahe();
        rho(j) = rhoC[0];
        mu(j) = muC[0];
        kappa(j) = kappaC[0];
        alpha(j) = alphaheC[0] / rho(j);
        D(j) = alpha(j) / Le;
    }
}

void ChemThermo::setY(const double* y)
{
    for (int k=0; k<nsp_; k++) {
        Y_[k] = y[k];
    }
}
