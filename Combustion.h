#ifndef CTF_COMBUSTION_H_
#define CTF_COMBUSTION_H_

// Direct integration of chemistry following OpenFOAM

#include <cmath>
#include <vector>
#include <Eigen/Dense>

#include "ChemThermo.h"

class Combustion
{
public:
    Combustion() = default;
    Combustion(const ChemThermo& thermo);

    struct stepState
    {
        const bool forward;
        double dxTry;
        double dxDid;
        bool first;
        bool last;
        bool reject;
        bool prevReject;
    
        stepState(const double dx)
        :
            forward(dx > 0 ? true : false),
            dxTry(dx),
            dxDid(0),
            first(true),
            last(false),
            reject(false),
            prevReject(false)
        {}
    };

    // Global solve
    double solve(const double& deltaT, const Eigen::VectorXd& T,
                 const std::vector<Eigen::VectorXd>& Y,
                 const Eigen::VectorXd& rho, const double& p);

    // Solve at each grid point
    void solve(std::vector<double>& c, double& T, double& p,
               const double& xEnd, double& dxTry) const;

    // Solve within each sub step
    void solve(double& x, std::vector<double>& y, stepState& step) const;

private:
    ChemThermo thermo_;

    int nsp_;

    std::vector<double> c_;
    std::vector<double> c0_;

    std::vector<Eigen::VectorXd> reactionRate_;

    Eigen::VectorXd deltaTChem_;

    double maxSteps_;
};


#endif  // CTF_COMBUSTION_H_
