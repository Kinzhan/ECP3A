//
//  YieldCurve.cpp
//  MyLibrary
//
//  Created by Alexandre HUMEAU on 30/08/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "YieldCurve.h"
#include "VectorUtilities.h"

namespace Finance {
    
    YieldCurve::YieldCurve()
    {
        eInterpolationType_ = Utilities::Interp::SPLINE_CUBIC;
    }
    
    YieldCurve::YieldCurve(const std::string & cCCY, const std::string & cName, const std::vector<std::pair<double, double> > & YC, Utilities::Interp::InterExtrapolationType eInterExtrapolationType) : 
    
    cCCY_(cCCY),
    cName_(cName)
    
    {
        std::pair<std::vector<double>, std::vector<double> > YC0 = Utilities::GetPairOfVectorFromVectorOfPair(YC);
        dVariables_ = YC0.first;
        dValues_ = YC0.second;
        eInterpolationType_ = eInterExtrapolationType;
    }
    
    YieldCurve::~YieldCurve()
    {}
    
    double YieldCurve::YC(double t) const
    {
        Utilities::require(t > 0, "t is not positive in YieldCurve::YC");
        return Interp1D(t);
    }
    
    std::string YieldCurve::GetCurrency() const
    {
        return cCCY_;
    }
    
    std::string YieldCurve::GetName() const
    {
        return cName_;
    }
}