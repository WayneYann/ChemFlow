#include <iostream>
#include <Eigen/Dense>

using namespace std;

class abc
{
public:
    typedef int label;
    typedef int label2;
    typedef Eigen::VectorXd scalarField;
    double dot(double x) const;

};

double abc::dot(double x) const
{
    for (label i=0; i<2; i++) {
        x*=x;
    }
    for (label2 i=0; i<2; i++) {
        x/=2;
    }
    return x;
}

int main()
{
    abc func;
    double x =2.0;
    cout << x << " " << func.dot(x) << endl;

    return 0;
}
