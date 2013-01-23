/*
 *  2DHullWhiteTS.cpp
 *  Seminaire
 *
 *  Created by Emile on 1/23/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include "2DHullWhiteTS.h"
#include "MathFunctions.h" // for BETAOUTHRESHOLD

namespace Maths {
    TwoDimHullWhiteTS::TwoDimHullWhiteTS(const Finance::TermStructure<double,double> & sTermStructure1, const Finance::TermStructure<double,double> & sTermStructure2/*, const double dLambda1, const double dLambda2*/)
    {
		Finance::TermStructure <double, double> sTermStructureProduct, sTermStructure2_ ;
        sTermStructureProduct = sTermStructure1;
        //  Little trick to make it build
        sTermStructure2_ = sTermStructure2;
		sTermStructureProduct *= sTermStructure2_ ;
		
        //	Initialize data member of termstructure integral
		
		UValues_ = sTermStructureProduct.GetValues() ;
		TVariables_ = sTermStructureProduct.GetVariables();
	}
    
    TwoDimHullWhiteTS::~TwoDimHullWhiteTS()
    {}
    
    double TwoDimHullWhiteTS::TwoDimSubIntegral(const double dA, const double dB, const double dS1, const double dS2, const double dLambda1, const double dLambda2) const
    {
		
		if (dLambda1 + dLambda2 < BETAOUTHRESHOLD)
        {
            return dB - dA;
        }
        else
        {
            return 1 / dLambda1 / dLambda2 *(dB - dA - (exp(dLambda1*(dS1-dA))-exp(dLambda1*(dS1-dB)))/dLambda1
											 - (exp(dLambda2*(dS1-dA))-exp(dLambda2*(dS1-dB)))/dLambda2
											 + (exp(dLambda1*dS1+dLambda2*dS2-(dLambda1+dLambda2)*dA)-exp(dLambda1*dS1+dLambda2*dS2-(dLambda1+dLambda2)*dB))/(dLambda1+dLambda2)) ;
        }
    }
}