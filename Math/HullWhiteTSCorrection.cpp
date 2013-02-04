/*
 *  HullWhiteTSCorrection.cpp
 *  Seminaire
 *
 *  Created by Emile on 2/4/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include "HullWhiteTSCorrection.h"
#include "MathFunctions.h" // for BETAOUTHRESHOLD

namespace Maths {
    HullWhiteTSCorrection::HullWhiteTSCorrection(const Finance::TermStructure<double,double> & sTermStructure, const double dLambdaDiscount, const double dLambdaForward/*, const double dRho*/, const double dT1, const double dT2) : TermStructureIntegral(sTermStructure), dLambdaDiscount_(dLambdaDiscount), dLambdaForward_(dLambdaForward)/*, dRho_(dRho)*/, dT1_(dT1), dT2_(dT2)
    {}
    
    HullWhiteTSCorrection::~HullWhiteTSCorrection()
    {}
    
    double HullWhiteTSCorrection::SubIntegral(const double dA, const double dB) const
    {
        if (dLambdaDiscount_ + dLambdaForward_ < BETAOUTHRESHOLD)
        {
            return /*dRho_ * */(dT1_ - dT2_) * 0.5 * ((dT2_ - dA)*(dT2_ - dA) - (dT2_ - dB)*(dT2_ - dB));
        }
        else 
        {
            return /*dRho_ * */1.0 / dLambdaForward_ / dLambdaDiscount_ * (MathFunctions::SumExp(dLambdaForward_+dLambdaDiscount_, dT2_ - dB, dT2_ - dB)
															 + MathFunctions::SumExp(dLambdaForward_, dT1_ - dB, dT1_ - dA)
															 - MathFunctions::SumExp(dLambdaForward_, dT1_ - dB, dT1_ - dA)
															 - MathFunctions::SumExp(dLambdaForward_+dLambdaDiscount_, dA, dB)*exp(-dLambdaForward_*dT1_-dLambdaDiscount_*dT2_));
        }
    }
}