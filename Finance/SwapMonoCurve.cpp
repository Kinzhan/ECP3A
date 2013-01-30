//
//  SwapMonoCurve.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 30/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "SwapMonoCurve.h"
#include "Date.h"
#include "DiscountFactor.h"

namespace Finance{
    SwapMonoCurve::SwapMonoCurve(const Utilities::Date::MyDate & sStartSwap, const Utilities::Date::MyDate & sEndSwap, MyFrequency eFixedLegFrequency, MyBasis eBasis, const YieldCurve & sYieldCurve)
    : 
    Annuity(sStartSwap, sEndSwap, eBasis, eFixedLegFrequency, sYieldCurve), 
    sYieldCurve_(sYieldCurve)
    
    {}
    
    SwapMonoCurve::~SwapMonoCurve()
    {}
    
    double SwapMonoCurve::ComputeSwap() const
    {
        DF sDF(sYieldCurve_);
        
        return (sDF.DiscountFactor(sStart_) - sDF.DiscountFactor(sEnd_)) / ComputeAnnuity();
    }
}