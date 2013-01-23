/*
 *  2DHullWhiteTS.h
 *  Seminaire
 *
 *  Created by Emile on 1/23/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef Seminaire_2DHullWhiteTS_h
#define Seminaire_2DHullWhiteTS_h

#include "Integral.h"
#include "Require.h"

namespace Maths {
    class TwoDimHullWhiteTS : public TermStructureIntegral
    {
    protected:
        //double dLambda1_;
		//double dLambda2_;
    public:
        TwoDimHullWhiteTS(const Finance::TermStructure<double,double> & sTermStructure1, const Finance::TermStructure<double,double> & sTermStructure2/*, const double dLambda1, const double dLambda2*/);
        virtual ~TwoDimHullWhiteTS();
        
        //  primitive of function f() define in the above integral
		virtual double TwoDimSubIntegral(const double dA, const double dB, const double dS1, const double dS2, const double dLambda1, const double dLambda2) const ;
        
    };
}

#endif
