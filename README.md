# ChemFlow

**ChemFlow** solves the transient planar counterflow flame with complex chemistry. The set of reduced 1D equations are discretized using finite difference, which are then efficiently solved with the *Tridiagonal Matrix Algorithm (TDMA)*. For the stiff chemistry system, the chemistry solver of [OpenFOAM](https://openfoam.org/) is used, providing chemical time scales needed for adjustable time stepping.

## Theory

> Kee, R. J., Coltrin, M. E., Glarborg, P., and Zhu, H., *Chemically Reacting Flow: Theory and Practice*, 2nd ed., John Wiley and Sons, 2017, Chap. 6.

## Installation

This solver depends on [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [OpenFOAM-6](https://openfoam.org/), which can be easily installed using *apt-get* on Ubuntu.

```bash
git clone https://github.com/ZX114/ChemFlow.git YOUR_WORKING_DIR/ChemFlow
./Allgxx
```

## Sample

A standard OpenFOAM case directory with *chemkin*, *constant* and *system* is still needed for holding BasicChemistryModel&lt;rhoReactionThermo&gt;. The reaction mechanism should be supplied in *chemkin*, either in CHEMKIN or OpenFOAM format.

To run the program

```bash
YOUR_WORKING_DIR/ChemFlow/chemflow
```

A simple script *plotOmega.py* is provided to visualize the reaction source terms

```bash
python plotOmega.py
```
