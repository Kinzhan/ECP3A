//
//  HullWhiteTS.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 16/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "HullWhiteTS.h"
#include "MathFunctions.h" // for BETAOUTHRESHOLD

namespace Maths {
    HullWhiteTS::HullWhiteTS(const Finance::TermStructure<double,double> & sTermStructure, double dLambda) : TermStructureIntegral(sTermStructure), dLambda_(dLambda)
    {}
    
    HullWhiteTS::~HullWhiteTS()
    {}
    
    double HullWhiteTS::SubIntegral(double dA, double dB) const
    {
        if (dLambda_ < BETAOUTHRESHOLD)
        {
            return dB - dA;
        }
        else
        {
            return (exp(-dLambda_ * dA) - exp(-dLambda_ * dB)) / dLambda_;
        }
    }
}