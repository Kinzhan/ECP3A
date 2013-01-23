/*
 *  Weights.h
 *  Seminaire
 *
 *  Created by Emile on 1/23/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef Seminaire_Weights_h
#define Seminaire_Weights_h

#include "DiscountFactor.h"
#include <vector>
#include "Basis.h"
#include "Frequency.h"

namespace Finance {
	class Weights: public Finance::DF
	{
	public:
		Weights();
		// We will first assume that the dCoverage is computed with respect to the ACT / 365 convention
		Weights(const YieldCurve & sInitialYieldCurve, const std::vector<double> & dS);
        Weights(const YieldCurve & sInitialYieldCurve, double dStart, double dEnd, MyFrequency eFrequency, MyBasis eBasis);
		virtual ~Weights();
		virtual double GetWeight(const std::size_t iFixing) const;
		virtual std::vector <double> GetWeights() const;
	private:
		std::vector <double> dS_ ;
		//std::vector <double> dCoverage_ ;
		std::vector <double> dWeights_ ;
	};
}

#endif