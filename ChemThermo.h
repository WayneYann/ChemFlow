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

    int speciesIndex(const std::string& name) const;

    std::string speciesName(const int& k) const;

    double calcHs(const double& T, const double* y) const;

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


private:
    void setY(const double* y);

    Foam::autoPtr<Foam::rhoReactionThermo> pThermo_;
    Foam::rhoReactionThermo& thermo_;
    Foam::autoPtr<Foam::BasicChemistryModel<Foam::rhoReactionThermo>> pChemistry_;
    Foam::volScalarField rho_;
    Foam::volVectorField U_;
    Foam::BasicChemistryModel<Foam::rhoReactionThermo>& chemistry_;
    Foam::basicSpecieMixture& composition_;
    Foam::PtrList<Foam::volScalarField>& Y_;

    int nsp_;
    double p0_;
};
#endif  // CTF_CHEMTHERMO_H_
