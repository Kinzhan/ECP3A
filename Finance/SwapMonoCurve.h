//
//  SwapMonoCurve.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 30/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_SwapMonoCurve_h
#define Seminaire_SwapMonoCurve_h

#include "Annuity.h"
#include "YieldCurve.h"

namespace Finance {
    class SwapMonoCurve : public Annuity
    {
    protected:
        YieldCurve sYieldCurve_;
    public:
        SwapMonoCurve(const Utilities::Date::MyDate & sStartSwap, const Utilities::Date::MyDate & sEndSwap, MyFrequency eFixedLegFrequency, MyBasis eBasis, const YieldCurve & sYieldCurve);
        virtual ~SwapMonoCurve();
        
        virtual double ComputeSwap() const;
    
    };
}

#endif
