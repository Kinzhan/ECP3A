//
//  HullWhiteTS.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 16/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Integral.h"
#include "Require.h"

namespace Maths {
	TermStructureIntegral::TermStructureIntegral() {}
	
	TermStructureIntegral::TermStructureIntegral(const Finance::TermStructure<double,double> & sTermStructure) {
		SetValues(sTermStructure.GetValues());
		SetVariables(sTermStructure.GetVariables());
	}
	
	TermStructureIntegral::~TermStructureIntegral() {}
	
	// computes the integral of TermStructure * f on [dT1, dT2]
	double TermStructureIntegral::Integral(const double dT1, const double dT2) const {
		std::vector<double> dTSVariables = GetVariables(), dTSValues = GetValues();
		Utilities::require(dT1 < dT2, "First boundary must be smaller than second boundary.");
		std::size_t iSize = dTSValues.size();
		
		if (iSize == 1) {
			return dTSValues[0] * SubIntegral(dT1, dT2);
		}
		else {
			double dInf = 0.0, dSup = 0.0, dIntegral = 0.0 ;
			
			// beginning
			if (dT1 < dTSVariables[0]) {
				dIntegral += dTSValues[0] * SubIntegral(dT1, std::min(dTSVariables[0], dT2));
			}
			
			// middle
			for (std::size_t iTS = 1; iTS < iSize; ++iTS) {
				dInf = std::max(dTSVariables[iTS-1], dT1);
				dSup = std::min(dTSVariables[iTS], dT2);
				if (dInf < dSup) {
					dIntegral += dTSValues[iTS-1] * SubIntegral(std::max(dTSVariables[iTS-1], dT1), std::min(dTSVariables[iTS], dT2));
				}
			}
			
			// end
			if (dT2 > dTSVariables[iSize - 1]) {
				dIntegral += dTSValues[iSize - 1] * SubIntegral(std::max(dTSVariables[iSize - 1], dT1), dT2);
                
			}
			
			return dIntegral;
		}
	}
}