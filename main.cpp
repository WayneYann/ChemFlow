#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <string>
#include <Eigen/Dense>

#include "tdma.h"
#include "chem_thermo.h"

using namespace std;
using namespace Eigen;

int main()
{
    // Discretize space and time
    const int nx = 201;
    const double XBEG = 0.0;
    const double XEND = 0.02;
    const double dx = (XEND - XBEG) / (nx - 1);
    VectorXd x(nx);
    const int nt = 201;
    const double TBEG = 0.0;
    const double TEND = 0.1;
    const double dt = (TEND - TBEG) / (nt - 1);

    // BC
    const double a = 50.0;  // strain rate
    const double uL = 0.5*a*(XEND - XBEG);
    const double uR = -0.5*a*(XEND - XBEG);
    const double VL = a;
    const double VR = a;
    const double TI = 300;
    const double TL = 300;
    const double TR = 300;
    const double YO2Air = 0.23197;
    const double YN2Air = 0.75425;
    const double YARAir = 0.01378;
    const double YFUEL = 1.0;
    const double YL = 1.0;
    const double YR = 0.0;

    // Output
    ofstream fout("output.csv");

    // Solution and initial conditions
    const double p0 = 101325.0;
    ChemThermo gas("Ethanol_31.cti", p0);
    const int nsp = gas.nsp();  // number of species
    VectorXd u(nx);  // x-direction velocity [m/s]
    VectorXd V(nx);  // v/y [1/s]
    VectorXd T(nx);  // temperature [K]
    vector<VectorXd> Y(nsp);  // species mass fractions [-]
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
    }

    // Properties
    const double rhoInf = 1.173;
    const double Le = 1.0;
    VectorXd rho(nx);
    VectorXd mu(nx);
    VectorXd kappa(nx);
    VectorXd alpha(nx);
    VectorXd D(nx);
    gas.updateThermo(T, Y, Le, rho, mu, kappa, alpha, D);


    // Time marching
    clock_t startTime, endTime;
    startTime = clock();
    MatrixXd A(nx,nx);
    MatrixXd b(nx,1);
    for (int i=0; i<nt; i++) {
        cout << "Time =  " << TBEG+i*dt << setprecision(4) << endl;
 
        // V equation
        A.setZero();
        b.setZero();
        for (int j=1; j<nx-1; j++) {
            double a1 = -mu(j)*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? -dt*u(j)/dx : 0.0);
            double a2 = 1.0 + dt*V(j) + 2.0*mu(j)*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? dt*u(j)/dx : -dt*u(j)/dx);
            double a3 = -mu(j)*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? 0.0 : dt*u(j)/dx);
            A(j,j-1) = a1;
            A(j,j) = a2;
            A(j,j+1) = a3;
            b(j) = rhoInf*a*a*dt/rho(j) + V(j);
        }
        A(0,0) = 1.0;
        A(nx-1,nx-1) = 1.0;
        b(0) = VL;
        b(nx-1) = VR;
        V = tdma(A,b);
        VectorXd::Index loc;
        cout << "V.max   " << V.maxCoeff(&loc) << " @ position "
             << loc << endl;

        // Continuity equation
        // TODO: try different V-BC (not zero)
        // Using staggered grid
        MatrixXd Ac(nx-1,nx-1);
        MatrixXd bc(nx-1,1);
        Ac.setZero();
        bc.setZero();
        for (int j=1; j<nx-1; j++) {
            Ac(j,j-1) = -1.0;
            Ac(j,j) = 1.0;
            bc(j) = -dx*V(j);
        }
        Ac(0,0) = 1.0;
        bc(0) = uL - 0.5*dx*VL;
        VectorXd uStag = tdma(Ac,bc); // staggered grid of size nx-1
        for (int j=1; j<nx-1; j++) {
            u(j) = 0.5*(uStag(j-1)+uStag(j));
        }
        u(0) = uL; u(nx-1) = uR;

        // Upwind differencing
        // for (int j=1; j<nx-1; j++) {
        //     A(j,j-1) = (u(j) > 0.0 ? -1.0 : 0.0);
        //     A(j,j) = (u(j) > 0.0 ? 1.0 : -1.0);
        //     A(j,j+1) = (u(j) > 0.0 ? 0.0 : 1.0);
        //     b(j) = -dx*V(j);
        // }
        // A(0,0) = 1;
        // A(nx-1,nx-1) = 1;
        // b(0) = uL;
        // b(nx-1) = uR;
        // u = tdma(A,b);
        //
        // Central differencing
        // FAIL
        // A.setZero();
        // b.setZero();
        // for (int j=1; j<nx-1; j++) {
        //     A(j,j-1) = -1.0;
        //     A(j,j+1) = 1.0;
        //     b(j) = -2*dx*V(j);
        // }
        // A(0,0) = 1.0;
        // A(nx-1,nx-1) = 1.0;
        // b(0) = uL;
        // b(nx-1) = uR;
        // u = tdma(A,b);

        cout << "u.max   " << u.maxCoeff(&loc) << " @ position "
             << loc << endl;

        // Y equation
        A.setZero();
        b.setZero();
        for (int j=1; j<nx-1; j++) {
            double a1 = -D(j)*dt/(dx*dx) + (u(j) > 0.0 ? -dt*u(j)/dx : 0.0);
            double a2 = 1.0 + 2.0*D(j)*dt/(dx*dx) + (u(j) > 0.0 ? dt*u(j)/dx : -dt*u(j)/dx);
            double a3 = -D(j)*dt/(dx*dx) + (u(j) > 0.0 ? 0.0 : dt*u(j)/dx);
            A(j,j-1) = a1;
            A(j,j) = a2;
            A(j,j+1) = a3;
            b(j) = Y[0](j);
        }
        A(0,0) = 1.0;
        A(nx-1,nx-1) = 1.0;
        b(0) = YL;
        b(nx-1) = YR;
        Y[0] = tdma(A,b);
        cout << "Y.max   " << Y[0].maxCoeff(&loc) << " @ position "
             << loc << endl;

        // T eqaution
        A.setZero();
        b.setZero();
        for (int j=1; j<nx-1; j++) {
            double a1 = -alpha(j)*dt/(dx*dx) + (u(j) > 0.0 ? -dt*u(j)/dx : 0.0);
            double a2 = 1.0 + 2.0*alpha(j)*dt/(dx*dx) + (u(j) > 0.0 ? dt*u(j)/dx : -dt*u(j)/dx);
            double a3 = -alpha(j)*dt/(dx*dx) + (u(j) > 0.0 ? 0.0 : dt*u(j)/dx);
            A(j,j-1) = a1;
            A(j,j) = a2;
            A(j,j+1) = a3;
            b(j) = T(j);
        }
        A(0,0) = 1.0;
        A(nx-1,nx-1) = 1.0;
        b(0) = TL;
        b(nx-1) = TR;
        T = tdma(A,b);
        cout << "T.max   " << T.maxCoeff(&loc) << " @ position "
             << loc << endl;

        gas.updateThermo(T, Y, Le, rho, mu, kappa, alpha, D);


        cout << endl;
    }
    cout << "End" << endl;
    endTime = clock();
    cout << "Run time   " << double(endTime - startTime) / CLOCKS_PER_SEC
         << setprecision(6) << " s" << endl;


    // Output
    fout << "x (m),u (m/s), V (1/s), T (K), Y" << endl;
    for (int j=0; j<nx; j++) {
        fout << x(j) << setprecision(6) << ","
             << u(j) << setprecision(6) << ","
             << V(j) << setprecision(6) << ","
             << T(j) << setprecision(6) << ","
             << Y[0][j] << setprecision (6) << endl;
    }

    return 0;
}
