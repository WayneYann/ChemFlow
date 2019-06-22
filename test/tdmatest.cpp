#include <iostream>
#include <Eigen/Dense>
#include "../tdma.h"
using namespace std;
using namespace Eigen;

int main()
{
    MatrixXd A(5,5);
    VectorXd b(5);

    A << 1,0,0,0,0,
         1,2,3,0,0,
         0,1,3,5,0,
         0,0,5,3,-1,
         0,0,0,0,1;
    b << 1,0,-1,0,10;

    VectorXd x(5);
    x = tdma(A,b);
    cout << x << endl;
    return 0;
}
