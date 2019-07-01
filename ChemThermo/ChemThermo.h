#ifndef CTF_CHEMTHERMO_H_
#define CTF_CHEMTHERMO_H_
#include <vector>
#include <Eigen/Dense>
#include "fvCFD.H"
#include "rhoReactionThermo.H"
#include "BasicChemistryModel.H"
#include "reactingMixture.H"
#include "chemistrySolver.H"
#include "OFstream.H"
#include "thermoPhysicsTypes.H"
#include "basicSpecieMixture.H"
#include "cellModeller.H"


// Wrapper class for OpenFOAM BasicChemistryModel<rhoReactionThermo>
class ChemThermo
{
public:
    ChemThermo(Foam::fvMesh& mesh, Foam::Time& runTime, const double& p0);

    ChemThermo& operator=(const ChemThermo&) = delete;

    int nsp() const {
        return nsp_;
    }

    // Call OpenFOAM species() that returns the hash table of species
    int speciesIndex(const std::string& name) const;
    std::string speciesName(const int& k) const;

    // Return sensible enthalpy of the mixture [J/kg]
    // with given temperature and species mass fractions
    double calcHs(const double& T, const double* y) const;

    // Call OpenFOAM correct() to compute temperature
    // with given sensible enthalpy and species mass fractions
    void calcT(Eigen::VectorXd& T, const std::vector<Eigen::VectorXd>& Y,
               const Eigen::VectorXd& hs);

    // Transfer mass fractions field into a single C array
    // at point j for all species
    void massFractions(const std::vector<Eigen::VectorXd>& Y,
                       double* y, const int& j) const;

    // Update thermophysical properties
    void updateThermo(const Eigen::VectorXd& hs,
                      const std::vector<Eigen::VectorXd>& Y,
                      const double Le, Eigen::VectorXd& rho,
                      Eigen::VectorXd& mu, Eigen::VectorXd& kappa,
                      Eigen::VectorXd& alpha, Eigen::VectorXd& D);

    // Solve the stiff chemistry and return the chemical time scale
    double solve(const double& deltaT, const Eigen::VectorXd& hs,
                 const std::vector<Eigen::VectorXd>& Y,
                 std::vector<Eigen::VectorXd>& wdot, Eigen::VectorXd& qdot);


private:
    void setY(const double* y);

    void syncState();

    void filter(std::vector<Eigen::VectorXd>& wdot) const;

    Foam::autoPtr<Foam::rhoReactionThermo> pThermo_;
    Foam::rhoReactionThermo& thermo_;
    Foam::autoPtr<Foam::BasicChemistryModel<Foam::rhoReactionThermo>> pChemistry_;
    // Mixture properties
    Foam::volScalarField rho_;  // Density [kg/m3]
    Foam::volScalarField mu_;  // Dynamic viscosity [kg/m s]
    Foam::volScalarField kappa_;  // Thermal conductivity [J/m s K]
    Foam::volScalarField alphahe_;  // Thermal diffusivity for energy [kg/m/s]
    Foam::volVectorField U_;
    Foam::BasicChemistryModel<Foam::rhoReactionThermo>& chemistry_;
    Foam::basicSpecieMixture& composition_;
    Foam::PtrList<Foam::volScalarField>& Y_;

    int nsp_;
    double p0_;
};
#endif  // CTF_CHEMTHERMO_H_
