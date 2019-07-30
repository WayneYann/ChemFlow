// ChemFlow
// -- StrainedChemFlow
// -- A Segregated Solution to the Quasi-1D Counterflow Flame
// -- xu-zhang@sjtu.edu.cn
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <stdexcept>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>

#include "lininterp.h"
#include "tdma.h"
#include "ChemThermo/ChemThermo.h"

const double numlimSmall = 1e-8;
const double numlimSmallSmall = 1e-14;
const double numlimGreat = 1e12;

void restore(const ChemThermo& gas, const Eigen::VectorXd& x,
             Eigen::VectorXd& u, Eigen::VectorXd& V,
             Eigen::VectorXd& T, Eigen::VectorXd& hs,
             std::vector<Eigen::VectorXd>& Y)
{
    std::vector<double> x0;
    std::vector<double> u0;
    std::vector<double> V0;
    std::vector<double> T0;
    std::vector<std::vector<double> > Y0(gas.nsp());
    std::ifstream fin("initial_solution.csv");
    std::string line, str;
    std::getline(fin, line);
    if ((std::count(line.begin(),line.end(),',')+1) != gas.nsp()+6)
        throw std::runtime_error("Input x,u,V,rho,D,T,Y0,...,Ynsp-1");
    while (std::getline(fin, line)) {
        int n = 0;
        std::istringstream buffer(line);
        while(std::getline(buffer, str, ',')) {
            switch (n) {
                case 0:
                    x0.push_back(std::stold(str));
                    break;
                case 1:
                    u0.push_back(std::stold(str));
                    break;
                case 2:
                    V0.push_back(std::stold(str));
                    break;
                case 3:
                    break;
                case 4:
                    break;
                case 5:
                    T0.push_back(std::stold(str));
                    break;
                default:
                    Y0[n-6].push_back(std::stold(str));
                    break;
            }
            ++n;
        }
    }

    for (int j=0; j<x.size(); j++) {
        u(j) = lininterp(x(j), x0, u0);
        V(j) = lininterp(x(j), x0, V0);
        T(j) = lininterp(x(j), x0, T0);
        for (int k=0; k<gas.nsp(); k++) {
            Y[k](j) = lininterp(x(j), x0, Y0[k]);
        }
        double y[gas.nsp()];
        gas.massFractions(Y, y, j);
        hs(j) = gas.calcHs(T(j), y);
    }
}

void write(const int& iter, const ChemThermo& gas, const Eigen::VectorXd& x,
           const Eigen::VectorXd& u, const Eigen::VectorXd& V,
           const Eigen::VectorXd& rho, const Eigen::VectorXd& D,
           const Eigen::VectorXd& T, const std::vector<Eigen::VectorXd>& Y,
           const Eigen::VectorXd& qdot, const std::vector<Eigen::VectorXd>& wdot)
{
    std::stringstream ss;
    ss << iter;
    std::ofstream fout("data/output-"+ss.str()+".csv");
    std::ofstream rout("data/reaction-"+ss.str()+".csv");
    // Output
    fout << "x (m),u (m/s),V (1/s),rho (kg/m3),D (m2/s2),T (K)";
    for (int k=0; k<gas.nsp(); k++) {
        fout << "," << gas.speciesName(k);
    }
    fout << std::endl;
    for (int j=0; j<x.size(); j++) {
        fout << std::setprecision(8) << x(j) << ","
             << std::setprecision(8) << u(j) << ","
             << std::setprecision(8) << V(j) << ","
             << std::setprecision(8) << rho(j) << ","
             << std::setprecision(8) << D(j) << ","
             << std::setprecision(8) << T(j);
        for (int k=0; k<gas.nsp(); k++) {
            fout << "," << std::setprecision(8) << Y[k](j);
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

void write(const double& time, const ChemThermo& gas, const Eigen::VectorXd& x,
           const Eigen::VectorXd& u, const Eigen::VectorXd& V,
           const Eigen::VectorXd& rho, const Eigen::VectorXd& D,
           const Eigen::VectorXd& T, const std::vector<Eigen::VectorXd>& Y,
           const Eigen::VectorXd& qdot, const std::vector<Eigen::VectorXd>& wdot)
{
    std::stringstream ss;
    ss << time;
    std::ofstream fout("data/output-"+ss.str()+".csv");
    std::ofstream rout("data/reaction-"+ss.str()+".csv");
    // Output
    fout << "x (m),u (m/s),V (1/s),rho (kg/m3),D (m2/s2),T (K)";
    for (int k=0; k<gas.nsp(); k++) {
        fout << "," << gas.speciesName(k);
    }
    fout << std::endl;
    for (int j=0; j<x.size(); j++) {
        fout << std::setprecision(8) << x(j) << ","
             << std::setprecision(8) << u(j) << ","
             << std::setprecision(8) << V(j) << ","
             << std::setprecision(8) << rho(j) << ","
             << std::setprecision(8) << D(j) << ","
             << std::setprecision(8) << T(j);
        for (int k=0; k<gas.nsp(); k++) {
            fout << "," << std::setprecision(8) << Y[k](j);
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
    // input
    std::ifstream finp("input.txt");
    if (!finp) throw std::runtime_error("input.txt NOT FOUND!");
    std::map<std::string, std::string> dict;
    std::string name;
    std::string value;
    while (finp >> name) {
        finp >> value;
        dict[name] = value;
    }
    // Discretize space and time
    const int nx = std::stoi(dict["nPoints"]);
    const double XBEG = std::stod(dict["XBEG"]);
    const double XEND = std::stod(dict["XEND"]);
    const double TBEG = std::stod(dict["TBEG"]);
    const double TEND = std::stod(dict["TEND"]);
    const double dtMax = std::stod(dict["dtMax"]);
    // BC
    double a = std::stod(dict["strainRate"]);  // prescribed strain rate
    const double VL = a*0;
    const double VR = a*0;
    const double TI = std::stod(dict["TI"]);
    const double TL = std::stod(dict["TL"]);
    const double TR = std::stod(dict["TR"]);
    const double YO2Air = 0.23197;
    const double YN2Air = 0.75425;
    const double YARAir = 0.01378;
    const double YFUEL = 1.0;
    std::vector<double> YL;
    std::vector<double> YR;
    const std::string FUELNAME = dict["fuelName"];
    // IC
    const bool rstr = (dict["restore"] == "true" ? true : false);
    const bool ign = (dict["ignition"] == "true" ? true : false);
    const bool strain = (dict["strain"] == "true" ? true : false);
    const double ignBEGt = std::stod(dict["ignBEGt"]);
    const double ignENDt = std::stod(dict["ignENDt"]);
    const double ignHs = std::stod(dict["ignHs"]);
    const double aBEG = std::stod(dict["aBEG"]);
    const double aEND = std::stod(dict["aEND"]);
    const int writeFreq = std::stoi(dict["writeFreq"]);
    const double rhoInf = std::stod(dict["rhoInf"]);
    const double p0 = std::stod(dict["pressure"]);

    const double dx = (XEND - XBEG) / (nx - 1);
    const double tprecision = numlimSmall;
    Eigen::VectorXd x(nx);
    double dtChem = 1e-6;  // initial chemical time scale
    double dt = dtChem;
    double time = TBEG;

    // Output
    std::ofstream fm("data/monitor.csv");
    fm << "time (s),temperature (K)" << std::endl;
    const size_t WIDTH = 18;

    // Solution and initial conditions
    // Both std::vector and Eigen::VectorXd are used for the purpose of
    // distinguishing between containers for species and spatial grid points
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
    YL[gas.speciesIndex(FUELNAME)] = YFUEL;
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

    if (rstr) restore(gas, x, u, V, T, hs, Y);
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
        if (strain) {
            a = aBEG + (aEND-aBEG)/(TEND-TBEG)*(time-TBEG);  // strain the flame
            if (iter++%writeFreq == 0) write(iter, gas, x, u, V, rho, D, T, Y, qdot, wdot);
            //sleep(1);
        } else if (iter++%writeFreq == 0) {
            write(time, gas, x, u, V, rho, D, T, Y, qdot, wdot);
        }
        fm << std::setprecision(10) << time << ","
           << std::setprecision(10) << T.maxCoeff(&loc) << std::endl;
        // V equation
        A.setZero();
        b.setZero();
        for (int j=1; j<nx-1; j++) {
            const double mul = 0.5*(mu(j)+mu(j-1));
            const double mur = 0.5*(mu(j)+mu(j+1));
            A(j,j-1) = -mul*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? -dt*u(j)/dx : 0.0);
            A(j,j) = 1.0 + dt*V(j) + mul*dt/(rho(j)*dx*dx) + mur*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? dt*u(j)/dx : -dt*u(j)/dx);
            A(j,j+1) = -mur*dt/(rho(j)*dx*dx) + (u(j) > 0.0 ? 0.0 : dt*u(j)/dx);
            b(j) = rhoInf*a*a*dt/rho(j) + V(j);
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
        m.setZero();
        m(0) = rho(0) * u(0);
        for (int j=1; j<nx; j++) {
            double drhodt0 = (numlimGreat*(rho(j-1) - rhoPrev(j-1)))/(numlimGreat*dt);
            double drhodt1 = (numlimGreat*(rho(j) - rhoPrev(j)))/(numlimGreat*dt);
            drhodt0 = (dt > tprecision ? drhodt0 : 0.0);
            drhodt1 = (dt > tprecision ? drhodt1 : 0.0);
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
        // Correct
        for (int j=0; j<nx; j++) {
            double sumY = 0.0;
            for (int k=0; k<nsp; k++) {
                Y[k](j) = (Y[k](j) > 0.0 ? Y[k](j) : 0.0);
                Y[k](j) = (Y[k](j) < 1.0 ? Y[k](j) : 1.0);
                sumY += Y[k](j);
            }
            for (int k=0; k<nsp; k++) {
                Y[k](j) /= sumY;
            }
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
        if (ign && time > ignBEGt && time < ignENDt) {
            A(nx/2,nx/2-1) = 0.0;
            A(nx/2,nx/2) = 1.0;
            A(nx/2,nx/2+1) = 0.0;
            b(nx/2) = ignHs;
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
        std::cout << "Time =  " << time << std::setprecision(10) << std::endl;
        dt = std::min(dtChem, dtMax);
        if (time+tprecision > TEND) break;
        if (time+dt > TEND) dt = TEND - time;
        std::cout << std::endl;
    }
    std::cout << "End" << std::endl;
    endTime = std::clock();
    std::cout << "Run time   " << double(endTime - startTime) / CLOCKS_PER_SEC
              << std::setprecision(6) << " s" << std::endl;
    write(time, gas, x, u, V, rho, D, T, Y, qdot, wdot);

    return 0;
}
