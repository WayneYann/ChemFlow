#include "lininterp.h"

template<class Type>
Type lininterp(const double x, const std::vector<double>& xOld,
               const std::vector<Type>& yOld)
{
    int n = xOld.size();

    int lo = 0;
    for (lo=0; lo<n && xOld[lo]>x; ++lo)
    {}

    int low = lo;
    if (low < n)
    {
        for (int i=low; i<n; ++i)
        {
            if (xOld[i] > xOld[lo] && xOld[i] <= x)
            {
                lo = i;
            }
        }
    }

    int hi = 0;
    for (hi=0; hi<n && xOld[hi]<x; ++hi)
    {}

    int high = hi;
    if (high < n)
    {
        for (int i=high; i<n; ++i)
        {
            if (xOld[i] < xOld[hi] && xOld[i] >= x)
            {
                hi = i;
            }
        }
    }


    if (lo<n && hi<n && lo != hi)
    {
        return yOld[lo]
            + ((x - xOld[lo])/(xOld[hi] - xOld[lo]))*(yOld[hi] - yOld[lo]);
    }
    else if (lo == hi)
    {
        return yOld[lo];
    }
    else if (lo == n)
    {
        return yOld[hi];
    }
    else
    {
        return yOld[lo];
    }
}

double lininterp(const double x, const Eigen::VectorXd& xOld,
                 const Eigen::VectorXd& yOld)
{
    int n = xOld.size();

    int lo = 0;
    for (lo=0; lo<n && xOld[lo]>x; ++lo)
    {}

    int low = lo;
    if (low < n)
    {
        for (int i=low; i<n; ++i)
        {
            if (xOld[i] > xOld[lo] && xOld[i] <= x)
            {
                lo = i;
            }
        }
    }

    int hi = 0;
    for (hi=0; hi<n && xOld[hi]<x; ++hi)
    {}

    int high = hi;
    if (high < n)
    {
        for (int i=high; i<n; ++i)
        {
            if (xOld[i] < xOld[hi] && xOld[i] >= x)
            {
                hi = i;
            }
        }
    }


    if (lo<n && hi<n && lo != hi)
    {
        return yOld[lo]
            + ((x - xOld[lo])/(xOld[hi] - xOld[lo]))*(yOld[hi] - yOld[lo]);
    }
    else if (lo == hi)
    {
        return yOld[lo];
    }
    else if (lo == n)
    {
        return yOld[hi];
    }
    else
    {
        return yOld[lo];
    }
}

