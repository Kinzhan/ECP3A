//
//  StochasticBasisSpread.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 09/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "StochasticBasisSpread.h"
#include <cmath>

namespace Processes {
    StochasticBasisSpread::StochasticBasisSpread() : dSigma_(0), dLambda_(0)
    {}
    
    StochasticBasisSpread::~StochasticBasisSpread()
    {}
    
    double StochasticBasisSpread::Sigma_B(double t, double T) const
    {
        return dSigma_ * exp(-dLambda_ * (T - t));
    }
}