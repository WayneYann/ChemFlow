#ifndef CTF_TDMA_H_
#define CTF_TDMA_H_

#include <Eigen/Dense>

// The Tridiagonal Matrix Algorithm (TDMA)
// A * x = d
// [b0 c0 0  0  0  0  0 ] [x0] = [d0]
// [a1 b1 c1 0  0  0  0 ] [x1] = [d1]
// [0  a2 b2 c2 0  0  0 ] [x2] = [d2]
// [0  0  a3 b3 c3 0  0 ] [x3] = [d3]
// [0  0  0  a4 b4 c4 0 ] [x4] = [d4]
// [0  0  0  0  a5 b5 c5] [x5] = [d5]
// [0  0  0  0  0  a6 b6] [x6] = [d6]
// Diagonally dominant : |bi| > |ai| + |ci|
Eigen::VectorXd tdma(const Eigen::MatrixXd& A, const Eigen::MatrixXd& d);

#endif  // CTF_TDMA_H_
