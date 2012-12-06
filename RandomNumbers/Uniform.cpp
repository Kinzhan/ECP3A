//
//  Uniform.cpp
//  MyLibrary
//
//  Created by Alexandre HUMEAU on 27/06/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Uniform.h"
#include <tr1/random>
#include <time.h>

namespace RandomNumbers {

    //  Default constructor
    Uniform::Uniform() : dLeft_(0.0), dRight_(1.0), iNRealisations_(1), iAntitheticVariables_(0)
    {
        dRealisations_.resize(iNRealisations_);
    }

    //  Constructor
    Uniform::Uniform(double dLeft, double dRight, size_t iNRealisations, int iAntitheticVariables)
    {
        dLeft_          = dLeft;
        dRight_         = dRight;
        iAntitheticVariables_ = iAntitheticVariables;
        iNRealisations_ = iNRealisations * (iAntitheticVariables_ ? 2 : 1);
        dRealisations_.resize(iNRealisations_);
    }

    //  Destructor
    Uniform::~Uniform()
    {
        dRealisations_.clear();
    }

    void Uniform::GenerateUniform()
    {
        std::tr1::ranlux64_base_01 eng; // core engine class
        eng.seed(time(NULL));
        
        // Generation of uniforms via this distribution
        std::tr1::uniform_real<double> dist(dLeft_, dRight_);
        dist.reset();
        
        for (std::size_t i = 0 ; i < iNRealisations_ ; ++i)
        {
            dRealisations_[i] = dist(eng);
            if (iAntitheticVariables_)
            {
                //  Antithetic variables
                i++;
                dRealisations_[i] = dLeft_ + (dRight_ - dLeft_) * (1 - dRealisations_[i - 1]);
            }
        }
    }

}
