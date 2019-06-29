// ChemFlow
// -- A Segregated Solution to the Quasi-1D Counterflow Flame
// -- xu-zhang@sjtu.edu.cn
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <string>
#include <Eigen/Dense>

#include "tdma.h"
#include "ChemThermo/ChemThermo.h"

void write(const double& time, const ChemThermo& gas, const Eigen::VectorXd& x,
           const Eigen::VectorXd& u, const Eigen::VectorXd& V,
           const Eigen::VectorXd& rho, const Eigen::VectorXd& D,
           const Eigen::VectorXd& T, const std::vector<Eigen::VectorXd>& Y,
           const Eigen::VectorXd& qdot, const std::vector<Eigen::VectorXd>& wdot)
{
    std::stringstream ss;
    ss << time;
    std::ofstream fout("output-"+ss.str()+".csv");
    std::ofstream rout("reaction-"+ss.str()+".csv");
    // Output
    fout << "x (m),u (m/s),V (1/s),rho (kg/m3),D (m2/s2),T (K)";
    for (int k=0; k<gas.nsp(); k++) {
        fout << "," << gas.speciesName(k);
    }
    fout << std::endl;
    for (int j=0; j<x.size(); j++) {
        fout << std::setprecision(6) << x(j) << ","
             << std::setprecision(6) << u(j) << ","
             << std::setprecision(6) << V(j) << ","
             << std::setprecision(6) << rho(j) << ","
             << std::setprecision(6) << D(j) << ","
             << std::setprecision(6) << T(j);
        for (int k=0; k<gas.nsp(); k++) {
            fout << "," << std::setprecision(6) << Y[k](j);
        }
        fout << std::endl;
    }
    // Output reactions related quantities (source terms)
    rout << "x (m),Qdot (J/m3 s)";
    for (int k=0; k<gas.nsp(); k++) {
        rout << "," << gas.speciesName(k);
    }
    rout << std::endl;
    for (int j=0; j<x.size(); j++) {
        rout << std::setprecision(6) << x(j) << ","
             << std::setprecision(6) << qdot(j);
        for (int k=0; k<gas.nsp(); k++) {
            rout << "," << std::setprecision(6) << wdot[k](j);
        }
        rout << std::endl;
    }
}

int main(int argc, char *argv[])
{
    // Discretize space and time
    const int nx = 201;
    const double XBEG = 0.0;
    const double XEND = 0.02;
    const double dx = (XEND - XBEG) / (nx - 1);
    Eigen::VectorXd x(nx);
    const double TBEG = 0.0;
    const double TEND = 0.2;
    const double dtMax = 0.5e-4;
    const double tprecision = 1e-10;
    double dtChem = 0.5e-4;  // initial chemical time scale
    double dt = dtChem;
    double time = TBEG;

    // BC
    const double a = 50.0;  // prescribed strain rate
    const double VL = a*0;
    const double VR = a*0;
    const double TI = 300;
    const double TL = 300;
    const double TR = 300;
    const double YO2Air = 0.23197;
    const double YN2Air = 0.75425;
    const double YARAir = 0.01378;
    const double YFUEL = 1.0;
    std::vector<double> YL;
    std::vector<double> YR;


    // Output
    const size_t WIDTH = 18;

    // Solution and initial conditions
    // Both std::vector and Eigen::VectorXd are used for the purpose of
    // distinguishing between containers for species and spatial grid points
    const double rhoInf = 1.173;
    const double lambda = rhoInf*a*a;  // -1/y*dp/dy [Pa/m2]
    const double p0 = 101325.0;
    #include "createFields.H"
    ChemThermo gas(mesh, runTime, p0);
    const int nsp = gas.nsp();  // number of species
    Eigen::VectorXd u(nx);  // x-direction velocity [m/s]
    Eigen::VectorXd V(nx);  // v/y = dv/dy [1/s]
    Eigen::VectorXd T(nx);  // temperature [K]
    Eigen::VectorXd hs(nx);  // sensible enthalpy [J/kg]
    std::vector<Eigen::VectorXd> Y(nsp);  // species mass fractions [-]
    std::vector<Eigen::VectorXd> wdot(nsp);  // reaction rates [kg/m3 s]
    Eigen::VectorXd qdot(nx);  // heat source [J/m3 s]
    YL.resize(nsp, 0.0);
    YR.resize(nsp, 0.0);
    YL[gas.speciesIndex("CH4")] = YFUEL;
    YR[gas.speciesIndex("O2")] = YO2Air;
    YR[gas.speciesIndex("N2")] = YN2Air;
    YR[gas.speciesIndex("AR")] = YARAir;
    const double hsL = gas.calcHs(TL, YL.data());
    const double hsR = gas.calcHs(TR, YR.data());

    for (int j=0; j<nx; j++) {
        x(j) = XBEG + dx*j;
        u(j) = -a*(x(j) - 0.5*(XEND-XBEG));
        V(j) = a;
        T(j) = TI;
        for (int k=0; k<nsp; k++) {
            Y[k].resize(nx);
            wdot[k].resize(nx);
            Y[k](j) = 0.0;
            wdot[k](j) = 0.0;
        }
        Y[gas.speciesIndex("O2")](j) = YO2Air;
        Y[gas.speciesIndex("N2")](j) = YN2Air;
        Y[gas.speciesIndex("AR")](j) = YARAir;
        double y[nsp];
        gas.massFractions(Y, y, j);
        hs(j) = gas.calcHs(T(j), y);
        qdot(j) = 0.0;
    }

    // Properties
    const double Le = 1.0;
    Eigen::VectorXd rho(nx);
    Eigen::VectorXd rhoPrev(nx);
    Eigen::VectorXd mu(nx);
    Eigen::VectorXd kappa(nx);
    Eigen::VectorXd alpha(nx);
    Eigen::VectorXd D(nx);
    gas.updateThermo(hs, Y, Le, rho, mu, kappa, alpha, D);
    rhoPrev = rho;


    // Time marching
    clock_t startTime, endTime;
    startTime = std::clock();
    Eigen::MatrixXd A(nx,nx);
    Eigen::MatrixXd b(nx,1);
    Eigen::VectorXd m(nx);  // conservative form for continuity equation
    Eigen::VectorXd::Index loc;
    int iter = 0;
    while (true) {
        if (iter++%1000 == 0) write(time, gas, x, u, V, rho, D, T, Y, qdot, wdot);
        // V equation
        A.setZero();
        b.setZero();
        for (int j=1; j<nx-1; j++) {
            const double mul = 0.5*(mu(j)+mu(j-1));
            const double mur = 0.5*(mu(j)+mu(j+1));
            A(j,j-1) = -mul*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? -dt*u(j)/dx : 0.0);
            A(j,j) = 1.0 + dt*V(j) + mul*dt/(rho(j)*dx*dx) + mur*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? dt*u(j)/dx : -dt*u(j)/dx);
            A(j,j+1) = -mur*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? 0.0 : dt*u(j)/dx);
            b(j) = lambda*dt/rho(j) + V(j);
        }
        A(0,0) = 1.0;
        A(nx-1,nx-1) = 1.0;
        b(0) = VL;
        b(nx-1) = VR;
        V = tdma(A,b);
        std::cout << std::setw(WIDTH) << "V.max "
                  << std::setw(WIDTH) << V.maxCoeff(&loc) << " @ position "
                  << loc << std::endl;

        // Continuity equation
        // Propagate from left to right
        // TODO: limit drhodt
        m.setZero();
        m(0) = rho(0) * u(0);
        for (int j=1; j<nx; j++) {
            double drhodt0 = (rho(j-1) - rhoPrev(j-1))/dt;
            double drhodt1 = (rho(j) - rhoPrev(j))/dt;
            m(j) = m(j-1) + dx*(-0.5*(drhodt0+drhodt1) - 0.5*(rho(j-1)*V(j-1)+rho(j)*V(j)));
        }
        const double rhouOffset = (-rho(0)*m(nx-1) - rho(nx-1)*m(0)) / (rho(0) + rho(nx-1));
        m = m.array() + rhouOffset;
        u = m.cwiseQuotient(rho);
        std::cout << std::setw(WIDTH) << "u.max "
                  << std::setw(WIDTH) << u.maxCoeff(&loc) << " @ position "
                  << loc << std::endl;

        dtChem = gas.solve(dt, hs, Y, wdot, qdot);
        // Y equations
        for (int k=0; k<nsp; k++) {
            A.setZero();
            b.setZero();
            for (int j=1; j<nx-1; j++) {
                const double rhoDl = 0.5*(rho(j)*D(j)+rho(j-1)*D(j-1));
                const double rhoDr = 0.5*(rho(j)*D(j)+rho(j+1)*D(j+1));
                A(j,j-1) = -rhoDl*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? -dt*u(j)/dx : 0.0);
                A(j,j) = 1.0 + rhoDl*dt/(rho(j)*dx*dx) + rhoDr*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? dt*u(j)/dx : -dt*u(j)/dx);
                A(j,j+1) = -rhoDr*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? 0.0 : dt*u(j)/dx);
                b(j) = Y[k](j) + dt*wdot[k](j)/rho(j);
            }
            A(0,0) = 1.0;
            A(nx-1,nx-1) = 1.0;
            b(0) = YL[k];
            b(nx-1) = YR[k];
            Y[k] = tdma(A,b);
            std::cout << std::setw(WIDTH) << "Y-" + gas.speciesName(k) + ".max "
                      << std::setw(WIDTH) << Y[k].maxCoeff(&loc) << " @ position "
                      << loc << std::endl;
        }

        // Energy eqaution
        A.setZero();
        b.setZero();
        for (int j=1; j<nx-1; j++) {
            const double rhoAlphal = 0.5*(rho(j)*alpha(j)+rho(j-1)*alpha(j-1));
            const double rhoAlphar = 0.5*(rho(j)*alpha(j)+rho(j+1)*alpha(j+1));
            A(j,j-1) = -rhoAlphal*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? -dt*u(j)/dx : 0.0);
            A(j,j) = 1.0 + rhoAlphal*dt/(rho(j)*dx*dx) + rhoAlphar*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? dt*u(j)/dx : -dt*u(j)/dx);
            A(j,j+1) = -rhoAlphar*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? 0.0 : dt*u(j)/dx);
            b(j) = hs(j) + dt*qdot(j)/rho(j);
        }
        A(0,0) = 1.0;
        A(nx-1,nx-1) = 1.0;
        b(0) = hsL;
        b(nx-1) = hsR;
        // Ignition
        if (time > 0.05 && time < 0.07) {
            A(nx/2,nx/2-1) = 0.0;
            A(nx/2,nx/2) = 1.0;
            A(nx/2,nx/2+1) = 0.0;
            b(nx/2) = 2e6;
        }
        hs = tdma(A,b);
        gas.calcT(T, Y, hs);
        std::cout << std::setw(WIDTH) << "T.max "
                  << std::setw(WIDTH) << T.maxCoeff(&loc) << " @ position "
                  << loc << std::endl;

        rhoPrev = rho;
        gas.updateThermo(hs, Y, Le, rho, mu, kappa, alpha, D);

        time += dt;
        // Adjustable time step according to chemical time scale
        std::cout << "Time =  " << time << std::setprecision(6) << std::endl;
        dt = std::min(dtChem, dtMax);
        if (time+tprecision > TEND) break;
        if (time+dt > TEND) dt = TEND - time;
        std::cout << std::endl;
    }
    std::cout << "End" << std::endl;
    endTime = std::clock();
    std::cout << "Run time   " << double(endTime - startTime) / CLOCKS_PER_SEC
              << std::setprecision(6) << " s" << std::endl;

    return 0;
}
