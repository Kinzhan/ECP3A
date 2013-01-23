/*
 *  Weights.cpp
 *  Seminaire
 *
 *  Created by Emile on 1/23/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include "Weights.h"
#include "Require.h"

namespace Finance {
    Weights::Weights()
    {}
    
    Weights::Weights(const YieldCurve & sInitialYieldCurve, const std::vector<double> & dS) : DF(sInitialYieldCurve)
    {
		std::size_t iSizeS = dS_.size() ;
		//std::size_t iSizeCoverage = dCoverage_.size() ;
		double dAnnuity = 0. ;
		double dDF = 0. ;
		
		Utilities::require(iSizeS > 1, "Need more fixing dates.");
		//Utilities::require(iSizeS == iSizeCoverage + 1, "Wrong vector sizes.");
		
		dS_ = dS ;
		//dCoverage_ = dCoverage ;
		dWeights_.resize(iSizeS) ;
		
		for (std::size_t iFixing = 1; iFixing < iSizeS; ++iFixing) {
			// We will first assume that the dCoverage is computed with respect to the ACT / 365 convention
			dDF = (dS_[iFixing] - dS_[iFixing-1]) * DiscountFactor(dS_[iFixing]) ;
			dAnnuity += dDF ;
			dWeights_[iFixing] = dDF ;
		}
		
		for (std::size_t iFixg = 1; iFixg < iSizeS; ++iFixg) {
			dWeights_[iFixg] /= dAnnuity ;
		}
		
	}
    
    Weights::~Weights()
    {
		dS_.clear() ;
		//dCoverage_.clear() ;
		dWeights_.clear() ;
	}
    
    double Weights::GetWeight(const std::size_t iFixing) const
    {
        return dWeights_[iFixing] ;
    }
	
	std::vector <double> Weights::GetWeights() const
    {
        return dWeights_ ;
    }
}