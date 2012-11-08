//
//  Gaussian.cpp
//  MyLibrary
//
//  Created by Alexandre HUMEAU on 27/06/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <vector>
#include "Gaussian.h"
//  For Emile
//#include <random>
//  For Alexandre
#include <tr1/random> 

namespace RandomNumbers {

    //  Default constructor
    Gaussian1D::Gaussian1D() : dMean_(0.0), dStdDev_(1.0), iNRealisations_(1), iAntitheticVariables_(0)
    {
        dRealisations_.resize(iNRealisations_);
    }

    //  Constructor
    Gaussian1D::Gaussian1D(double dMean, double dStdDev, size_t iNRealisations, int iAntitheticVariables):
    dMean_(dMean),
    dStdDev_(dStdDev),
    iNRealisations_(iNRealisations),
    iAntitheticVariables_(iAntitheticVariables)
    {}

    //  Destructor
    Gaussian1D::~Gaussian1D()
    {
        dRealisations_.clear();
    }

    void Gaussian1D::GenerateGaussian()
    {
        std::tr1::ranlux64_base_01 eng; // core engine class
        eng.seed(time(NULL));
        
        //  Generation of normal variables via this distribution
        std::tr1::normal_distribution<double> dist(dMean_, dStdDev_ /** dStdDev_*/);
        dist.reset(); // discard any cached values
        
        std::vector<double> dRealisations(iNRealisations_, 0.0), dRealisationsAntithetic(iNRealisations_, 0.0);
        
        for (std::size_t i = 0 ; i < dRealisations.size() ; ++i)
        {
            dRealisations[i] = dist(eng);
            if (iAntitheticVariables_)
            {
                //  Antithetic variables 
                dRealisationsAntithetic[i] = 2 * dMean_ - dRealisations[i];
            }
        }
        
        dRealisations_.insert(dRealisations_.begin(), dRealisations.begin(), dRealisations.end());
        dRealisations_.insert(dRealisations_.end(), dRealisationsAntithetic.begin(), dRealisationsAntithetic.end());
    }
}
