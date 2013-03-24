//
//  HJM.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 24/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "HJM.h"

namespace Processes {
    HeathJarrowMorton::HeathJarrowMorton()
    {}
    
    HeathJarrowMorton::HeathJarrowMorton(const Finance::YieldCurve & sDiscountCurve) : 
    sDiscountCurve_(sDiscountCurve)
    {}
    
    HeathJarrowMorton::HeathJarrowMorton(const Finance::YieldCurve & sDiscountCurve, const Finance::YieldCurve & sForwardCurve) : 
    sDiscountCurve_(sDiscountCurve)
    {
        sForwardCurve_ = sForwardCurve;
    }
    
    HeathJarrowMorton::~HeathJarrowMorton()
    {}
    
    double HeathJarrowMorton::Libor(const double dt, const double dStart, const double dEnd, const double dX, const CurveName & eCurveName) const
    {
        //  Must change coverage to take into account real basis
        return 1.0 / (dEnd - dStart) * (BondPrice(dt, dStart, dX, eCurveName) / BondPrice(dt, dEnd, dX, eCurveName) - 1.0);
    }
}