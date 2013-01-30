//
//  Annuity.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 23/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_Annuity_h
#define Seminaire_Annuity_h

#include "Date.h"
#include "Basis.h"
#include "YieldCurve.h"
#include "Frequency.h"

namespace Finance {
    
    class Annuity
    {
    public:
        Annuity(const Utilities::Date::MyDate &   sStart,
                const Utilities::Date::MyDate &   sEnd,
                const MyBasis                     eBasis,
                const MyFrequency                 eFrequency,
                const YieldCurve &                sYieldCurve);
        virtual ~Annuity();
        
        virtual double ComputeAnnuity() const;
        
    protected:
        Utilities::Date::MyDate sStart_, sEnd_;
        MyBasis eBasis_;
        MyFrequency eFrequency_;
        YieldCurve sYieldCurve_;
        
    };
}

#endif
