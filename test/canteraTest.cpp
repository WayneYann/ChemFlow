#include <iostream>
#include "../ChemThermo.h"

using namespace std;
using namespace Cantera;

int main()
{
    ChemThermo gas("Ethanol_31.cti", 101325.0);

    double T = 300;
    cout << gas.ha(T, 0) << " "
         << gas.ha(T, 1) << " "
         << gas.ha(700, 9) << " "
         << endl;

    return 0;
}
