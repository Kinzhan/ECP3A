//
//  YieldCurve.h
//  FinanceTools
//
//  Created by Alexandre HUMEAU on 04/06/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef FinanceTools_YieldCurve_h
#define FinanceTools_YieldCurve_h

#include "InterExtrapolation.h"
#include "TermStructure.h"

namespace Finance {
    
    class YieldCurve : public Utilities::Interp::InterExtrapolation1D
    {
    protected:
        std::string cCCY_;
        std::string cName_;
        
    public:
        YieldCurve();
        YieldCurve(const std::string & cCCY, const std::string & cName, const std::vector<std::pair<double, double> > & YC, Utilities::Interp::InterExtrapolationType eInterExtrapolationType = Utilities::Interp::SPLINE_CUBIC);
        virtual ~YieldCurve();
        
        virtual std::string GetCurrency() const;
        virtual std::string GetName() const;
        
        virtual double YC(double t) const;
        
        virtual YieldCurve operator + (const YieldCurve & sYieldCurve);
        virtual YieldCurve operator = (double dValue);
        
    };
    
}

#endif
