//
//  Gaussian.h
//  MyLibrary
//
//  Created by Alexandre HUMEAU on 27/06/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <vector>

#ifndef GAUSSIAN_H_INCLUDED
#define GAUSSIAN_H_INCLUDED

namespace RandomNumbers
{
    //  The following class is used to generate gaussian random numbers
    class Gaussian1D
    {
    public:
        Gaussian1D();
        Gaussian1D(double dMean, double dStdDev, size_t iNRealisations, int iAntitheticVariables);
        ~Gaussian1D();

        virtual void GenerateGaussian();
        virtual std::vector<double> GetRealisations() const
        {
            return dRealisations_;
        };
        virtual size_t GetNbRealisations() const
        {
            return iNRealisations_ * (iAntitheticVariables_ ? 2 : 1);
        };

    protected:
        double dMean_;
        double dStdDev_;
        size_t iNRealisations_;
        int iAntitheticVariables_;
        std::vector<double> dRealisations_;

    };
}

#endif // GAUSSIAN_H_INCLUDED
