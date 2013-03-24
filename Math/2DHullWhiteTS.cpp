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
		//  this function now computes : 
        // \int_{A}^{B} (1 - e^{-\lambda_1 (S_1 - u)}) (1 - e^{-\lambda_2 (S_2 - u)}) du / (\lambda_1 * \lambda_2)
		if (dLambda1 + dLambda2 < BETAOUTHRESHOLD)
        {
            std::cout << "Lambda1+ LAmbda2 is negative" << std::endl;
            return (dB - dA) * dS1 * dS2 - (dS1 + dS2) * 0.5 * (dB * dB - dA * dA) + 1 / 3. * (dB * dB * dB - dA * dA * dA);
        }
        else
        {
            double dFirstTerm = dB - dA,
            dSecondTerm = (exp(-dLambda1 * (dS1 - dB)) - exp(-dLambda1 * (dS1 - dA))) / dLambda1,
            dThirdTerm = (exp(-dLambda2 * (dS2 - dB)) - exp(-dLambda2 * (dS2 - dA))) / dLambda2,
            dFourthTerm = exp(-dLambda1 * dS1 - dLambda2 * dS2) * (exp((dLambda1 + dLambda2) * dB) - exp((dLambda1 + dLambda2) * dA)) / (dLambda1 + dLambda2);
            return 1 / dLambda1 / dLambda2 * (dFirstTerm - dSecondTerm - dThirdTerm + dFourthTerm);
        }
    }
	
	double TwoDimHullWhiteTS::Integral(const double dT1, const double dT2, const double dS1, const double dS2, const double dLambda1, const double dLambda2) const
	{
        //  This function now computes \int_{T1}^{T2} \Gamma_1(u,S_1) \Gamma_2(u, S_2) du
        
		std::vector<double> dTSVariables = GetVariables(), dTSValues = GetValues();
        //  if T1 and T2 are too close the integral should be 0
        //  A.H. 20.02.2012
        if (std::abs(dT1 - dT2) < 1e-07)
        {
            return 0.0;
        }
		Utilities::require(dT1 < dT2, "First boundary must be smaller than second boundary.");
		std::size_t iSize = dTSValues.size();
		
		if (iSize == 1) 
        {
			return dTSValues[0] * TwoDimSubIntegral(dT1, dT2, dS1, dS2, dLambda1, dLambda2);
		}
		else 
        {
			double dInf = 0.0, dSup = 0.0, dIntegral = 0.0 ;
			
			// beginning
			if (dT1 < dTSVariables[0]) 
            {
				dIntegral += dTSValues[0] * TwoDimSubIntegral(dT1, std::min(dTSVariables[0], dT2), dS1, dS2, dLambda1, dLambda2);
			}
			
			// middle
			for (std::size_t iTS = 1; iTS < iSize; ++iTS) 
            {
				dInf = std::max(dTSVariables[iTS-1], dT1);
				dSup = std::min(dTSVariables[iTS], dT2);
				if (dInf < dSup) 
                {
					dIntegral += dTSValues[iTS-1] * TwoDimSubIntegral(dInf, dSup, dS1, dS2, dLambda1, dLambda2);
				}
                else 
                {
                    break;
                }
			}
			
			// end
			if (dT2 > dTSVariables[iSize - 1]) 
            {
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