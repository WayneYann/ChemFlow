#include <iostream>
#include <cmath>

#include "tdma.h"

// TODO: implement sparse matrix for A
Eigen::VectorXd tdma(const Eigen::MatrixXd& A, const Eigen::MatrixXd& d)
{
    int nx = d.size();
    Eigen::VectorXd x(nx);
    Eigen::VectorXd a(nx);
    Eigen::VectorXd b(nx);
    Eigen::VectorXd c(nx);
    Eigen::VectorXd u(nx);
    Eigen::VectorXd q(nx);

    // TODO: check if A is tridiagonal matrix
    // Group coefficients
    a(0) = 0.0;
    b(0) = A(0,0);
    c(0) = A(0,1);
    a(nx-1) = A(nx-1,nx-2);
    b(nx-1) = A(nx-1,nx-1);
    c(nx-1) = 0.0;
    for (int i=1; i<nx-1; i++) {
        a(i) = A(i,i-1);
        b(i) = A(i,i);
        c(i) = A(i,i+1);
    }

    bool DIAGDOM = true;
    for (int i=0; i<nx; i++) {
        if (std::abs(b(i)) < (std::abs(a(i))+std::abs(c(i)))) DIAGDOM = false;
    }
    if (!DIAGDOM) {
        std::cout << "Using QR decomposition with column pivoting..."
                  << std::endl;
        return A.colPivHouseholderQr().solve(d);
    }

    // Step1
    u(0) = c(0)/b(0);
    q(0) = d(0)/b(0);
    for (int i=1; i<nx; i++) {
        double m = 1.0/(b(i)-u(i-1)*a(i));
        u(i) = c(i)*m;
        q(i) = (d(i)-q(i-1)*a(i))*m;
    }

    // Step2
    x(nx-1) = q(nx-1);
    for (int i=nx-2; i>=0; i--) {
        x(i) = q(i) - u(i)*x(i+1);
    }

    return x;
}

