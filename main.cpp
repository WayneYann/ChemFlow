#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <string>
#include <Eigen/Dense>

#include "tdma.h"
#include "ChemThermo.h"
#include "Combustion.h"

using namespace std;
using namespace Eigen;

int main()
{
    // Discretize space and time
    const int nx = 101;
    const double XBEG = 0.0;
    const double XEND = 0.02;
    const double dx = (XEND - XBEG) / (nx - 1);
    VectorXd x(nx);
    const int nt = 401;
    const double TBEG = 0.0;
    const double TEND = 0.2;
    const double dt = (TEND - TBEG) / (nt - 1);

    // BC
    const double a = 50.0;  // initial strain rate
    const double VL = a*0;
    const double VR = a*0;
    const double TI = 1200;
    const double TL = 1200;
    const double TR = 1200;
    const double YO2Air = 0.23197;
    const double YN2Air = 0.75425;
    const double YARAir = 0.01378;
    const double YFUEL = 1.0;
    vector<double> YL;
    vector<double> YR;


    // Output
    ofstream fout("output.csv");
    ofstream omegaOut("omegaOut.csv");
    const size_t WIDTH = 18;

    // Solution and initial conditions
    const double rhoInf = 1.173;
    const double p0 = 101325.0;
    ChemThermo gas("Ethanol_31.cti", p0);
    Combustion combustion(gas);
    const int nsp = gas.nsp();  // number of species
    VectorXd u(nx);  // x-direction velocity [m/s]
    VectorXd V(nx);  // v/y = dv/dy [1/s]
    const double lambda = rhoInf*a*a;  // -1/y*dp/dy [Pa/m2]
    VectorXd T(nx);  // temperature [K]
    VectorXd hs(nx);  // sensible enthalpy [J/kg]
    vector<VectorXd> Y(nsp);  // species mass fractions [-]
    YL.resize(nsp, 0.0);
    YR.resize(nsp, 0.0);
    YL[gas.speciesIndex("C2H5OH")] = YFUEL;
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
            Y[k](j) = 0.0;
        }
        Y[gas.speciesIndex("O2")](j) = YO2Air;
        Y[gas.speciesIndex("N2")](j) = YN2Air;
        Y[gas.speciesIndex("AR")](j) = YARAir;
        double y[nsp];
        gas.massFractions(Y, y, j);
        hs(j) = gas.calcHs(T(j), y);
    }

    // Properties
    const double Le = 1.0;
    VectorXd rho(nx);
    VectorXd rhoPrev(nx);
    VectorXd mu(nx);
    VectorXd kappa(nx);
    VectorXd alpha(nx);
    VectorXd D(nx);
    gas.updateThermo(T, Y, Le, rho, mu, kappa, alpha, D);
    rhoPrev = rho;


    // Time marching
    clock_t startTime, endTime;
    startTime = clock();
    MatrixXd A(nx,nx);
    MatrixXd b(nx,1);
    VectorXd m(nx);  // conservative form for continuity equation
    VectorXd::Index loc;
    for (int i=0; i<nt; i++) {
        cout << "Time =  " << TBEG+i*dt << setprecision(4) << endl;
 
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
        cout << setw(WIDTH) << "V.max "
             << setw(WIDTH/2) << V.maxCoeff(&loc) << " @ position "
             << loc << endl;

        // Continuity equation
        // Propagate from left to right
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
        cout << setw(WIDTH) << "u.max "
             << setw(WIDTH/2) << u.maxCoeff(&loc) << " @ position "
             << loc << endl;

        combustion.solve(dt, T, Y, p0);
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
                b(j) = Y[k](j);
            }
            A(0,0) = 1.0;
            A(nx-1,nx-1) = 1.0;
            b(0) = YL[k];
            b(nx-1) = YR[k];
            Y[k] = tdma(A,b);
            cout << setw(WIDTH) << "Y-" + gas.speciesName(k) + ".max "
                 << setw(WIDTH/2) << Y[k].maxCoeff(&loc) << " @ position "
                 << loc << endl;
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
            b(j) = hs(j);
        }
        A(0,0) = 1.0;
        A(nx-1,nx-1) = 1.0;
        b(0) = hsL;
        b(nx-1) = hsR;
        hs = tdma(A,b);
        gas.calcT(T, Y, hs);
        cout << setw(WIDTH) << "T.max "
             << setw(WIDTH/2) << T.maxCoeff(&loc) << " @ position "
             << loc << endl;

        rhoPrev = rho;
        gas.updateThermo(T, Y, Le, rho, mu, kappa, alpha, D);

        cout << endl;
    }
    cout << "End" << endl;
    endTime = clock();
    cout << "Run time   " << double(endTime - startTime) / CLOCKS_PER_SEC
         << setprecision(6) << " s" << endl;


    // Output
    fout << "x (m),u (m/s),V (1/s),rho (kg/m3),D (m2/s2),T (K)";
    for (int k=0; k<nsp; k++) {
        fout << "," << gas.speciesName(k);
    }
    fout << endl;
    for (int j=0; j<nx; j++) {
        fout << setprecision(6) << x(j) << ","
             << setprecision(6) << u(j) << ","
             << setprecision(6) << V(j) << ","
             << setprecision(6) << rho(j) << ","
             << setprecision(6) << D(j) << ","
             << setprecision(6) << T(j);
        for (int k=0; k<nsp; k++) {
            fout << "," << setprecision(6) << Y[k](j);
        }
        fout << endl;
    }

    return 0;
}
