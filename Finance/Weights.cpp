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
#include "Schedule.h"

namespace Finance {
    Weights::Weights()
    {}
    
    Weights::Weights(const YieldCurve & sInitialYieldCurve, const std::vector<double> & dS) : DF(sInitialYieldCurve), dS_(dS)
    {
		std::size_t iSizeS = dS_.size() ;
		double dAnnuity = 0. ;
		
		Utilities::require(iSizeS > 1, "Need more fixing dates.");
		
		dWeights_.resize(iSizeS - 1) ;
		
		for (std::size_t iFixing = 1; iFixing < iSizeS; ++iFixing) {
			// We will first assume that the dCoverage is computed with respect to the ACT / 365 convention
            double dDF = (dS_[iFixing] - dS_[iFixing - 1]) * DiscountFactor(dS_[iFixing]) ;
			dAnnuity += dDF ;
			dWeights_[iFixing - 1] = dDF ;
		}
		
		for (std::size_t iFixing = 0; iFixing < dWeights_.size() ; ++iFixing) {
			dWeights_[iFixing] /= dAnnuity ;
		}
		
	}
	
	Weights::Weights(const YieldCurve & sInitialYieldCurve, const std::vector<double> & dT, const std::vector<double> & dS) : DF(sInitialYieldCurve), dS_(dT)
    {
		std::size_t iSizeT = dS_.size() ;
		std::size_t iSizeS = dS.size() ;
		double dAnnuity = 0. ;
		
		Utilities::require(iSizeT > 1 && iSizeS > 1, "Need more fixing dates.");
		
		dWeights_.resize(iSizeT) ;
		
		for (std::size_t iFixing = 1; iFixing < iSizeT; ++iFixing) {
			// We will first assume that the dCoverage is computed with respect to the ACT / 365 convention
			double dDF = (dT[iFixing] - dT[iFixing-1]) * DiscountFactor(dT[iFixing]) ;
			dWeights_[iFixing] = dDF ;
		}
		
		for (std::size_t iFixing = 1; iFixing < iSizeS; ++iFixing) {
			dAnnuity += (dS_[iFixing] - dS_[iFixing-1]) * DiscountFactor(dS_[iFixing]) ;
		}
		
		for (std::size_t iFixing = 1; iFixing < iSizeT; ++iFixing) {
			dWeights_[iFixing] /= dAnnuity ;
		}
		
	}
    
    Weights::Weights(const YieldCurve & sInitialYieldCurve, double dStart, double dEnd, MyFrequency eFrequency, MyBasis eBasis) : DF(sInitialYieldCurve)
    {
        Utilities::require(dEnd > dStart, "Start is after end");
        Utilities::Date::MyDate sStart(dStart), sEnd(dEnd);
        
        Schedule sSchedule(sStart, sEnd, sInitialYieldCurve, eBasis, eFrequency);
        std::vector<EventOfSchedule> sVectOfEventOfSchedule = sSchedule.GetSchedule();
        dWeights_.resize(sVectOfEventOfSchedule.size());
        double dAnnuity = 0;
        
        std::size_t iSizeOfSchedule = sVectOfEventOfSchedule.size();
        
        for (std::size_t i = 0 ; i < iSizeOfSchedule ; ++i)
        {
            double dDF = sVectOfEventOfSchedule[i].GetCoverage() * sVectOfEventOfSchedule[i].GetPayingDateDF();
            dAnnuity += dDF;
            dWeights_[i] = dDF;
        }
        
        for (std::size_t i = 0 ; i < iSizeOfSchedule ; ++i)
        {
            dWeights_[i] /= dAnnuity;
        }
    }
    
    Weights::~Weights()
    {
		dS_.clear() ;
		dWeights_.clear() ;
	}
    
    double Weights::GetWeight(std::size_t iFixing) const
    {
        return dWeights_[iFixing] ;
    }
	
	std::vector <double> Weights::GetWeights() const
    {
        return dWeights_ ;
    }
}