#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main()
{
    VectorXd x(5);
    x << 1, 2, 3, 4, 5;
    MatrixXd y(3,3);
    y << 1,3,5,
         4,9,7,
         5,2,8;
    const double* a = y.data();
    for (int i=0; i<y.size(); i++) {
        cout << a[i] << " ";
    }
    cout << endl;
    return 0;
}
