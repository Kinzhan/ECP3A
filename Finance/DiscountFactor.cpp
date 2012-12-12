/*
 *  DiscountFactor.cpp
 *  Seminaire
 *
 *  Created by Emile on 12/12/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "DiscountFactor.h"
#include "MathFunctions.h"

namespace Finance {
	DF::DF()
	{}
	
	DF::DF(const YieldCurve & sInitialYieldCurve) : YieldCurve(sInitialYieldCurve)
	{}
	
	DF::~DF()
	{}
	
	double DF::DiscountFactor(const double dT) const
	{
		return exp(-dT * YC(dT));
	}
}
