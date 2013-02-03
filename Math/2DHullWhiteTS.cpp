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
            return (dB - dA) * dS1 * dS2 - (dS1 + dS2) * 0.5 * (dB * dB - dA * dA) + 1 / 3. * (dB * dB * dB - dA * dA * dA);
        }
        else
        {
            return 1 / dLambda1 / dLambda2 *(dB - dA - (exp(-dLambda1*(dS1-dA))-exp(-dLambda1*(dS1-dB)))/dLambda1
											 - (exp(dLambda2*(-dS1-dA))-exp(dLambda2*(-dS1-dB)))/dLambda2
											 + (exp(-dLambda1*dS1-dLambda2*dS2+(dLambda1+dLambda2)*dA)-exp(-dLambda1*dS1-dLambda2*dS2+(dLambda1+dLambda2)*dB))/(dLambda1+dLambda2)) ;
        }
    }
	
	double TwoDimHullWhiteTS::Integral(const double dT1, const double dT2, const double dS1, const double dS2, const double dLambda1, const double dLambda2) const
	{
		std::vector<double> dTSVariables = GetVariables(), dTSValues = GetValues();
		Utilities::require(dT1 < dT2, "First boundary must be smaller than second boundary.");
		std::size_t iSize = dTSValues.size();
		
		if (iSize == 1) {
			return dTSValues[0] * TwoDimSubIntegral(dT1, dT2, dS1, dS2, dLambda1, dLambda2);
		}
		else {
			double dInf = 0.0, dSup = 0.0, dIntegral = 0.0 ;
			
			// beginning
			if (dT1 < dTSVariables[0]) {
				dIntegral += dTSValues[0] * TwoDimSubIntegral(dT1, std::min(dTSVariables[0], dT2), dS1, dS2, dLambda1, dLambda2);
			}
			
			// middle
			for (std::size_t iTS = 1; iTS < iSize; ++iTS) {
				dInf = std::max(dTSVariables[iTS-1], dT1);
				dSup = std::min(dTSVariables[iTS], dT2);
				if (dInf < dSup) {
					dIntegral += dTSValues[iTS-1] * TwoDimSubIntegral(std::max(dTSVariables[iTS-1], dT1), std::min(dTSVariables[iTS], dT2), dS1, dS2, dLambda1, dLambda2);
				}
			}
			
			// end
			if (dT2 > dTSVariables[iSize - 1]) {
				dIntegral += dTSValues[iSize - 1] * TwoDimSubIntegral(std::max(dTSVariables[iSize - 1], dT1), dT2, dS1, dS2, dLambda1, dLambda2);
			}
			return dIntegral;
		}
	}
	
	double TwoDimHullWhiteTS::SubIntegral(const double dA, const double dB) const
	{
		return 0.0 ;
	}
}