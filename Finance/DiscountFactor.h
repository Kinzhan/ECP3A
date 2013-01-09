/*
 *  DiscountFactor.h
 *  Seminaire
 *
 *  Created by Emile on 12/12/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef Seminaire_DiscountFactor_h
#define Seminaire_DiscountFactor_h

#include "YieldCurve.h"

namespace Finance {
	class DF: public YieldCurve
	{
	public:
		DF();
		DF(const YieldCurve & sInitialYieldCurve);
		virtual ~DF();
		virtual double DiscountFactor(const double dDate) const;
	private:	
	};
}

#endif