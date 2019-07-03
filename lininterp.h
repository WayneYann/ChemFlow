#ifndef CTF_LININTERP_H_
#define CTF_LININTERP_H_

#include <vector>
#include <Eigen/Dense>

template<class Type>
Type lininterp(const double x, const std::vector<double>& xOld,
               const std::vector<Type>& yOld);

double lininterp(const double x, const Eigen::VectorXd& xOld,
                 const Eigen::VectorXd& yOld);

#ifdef NoRepository
    #include "lininterp.cpp"
#endif

#endif  // CTF_LININTERP_H_
