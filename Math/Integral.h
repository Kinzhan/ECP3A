/*
 *  Integral.h
 *  Seminaire
 *
 *  Created by Emile on 1/16/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef Seminaire_Integral_h
#define Seminaire_Integral_h

#include "TermStructure.h"

namespace Maths {
	class TermStructureIntegral: public Finance::TermStructure<double,double> {
	protected:
		//double dT1_;
		//double dT2_;
		//double fSubIntegral(const double dA, const double dB) const;
	public:
		TermStructureIntegral();
		TermStructureIntegral(const Finance::TermStructure<double,double> & sTermStructure);
		virtual ~TermStructureIntegral();
        
        // compute \int_{T1}{T2} \sigma(u) f(u) du
		virtual double Integral(const double dT1, const double dT2) const;
        
        //  primitive of function f() define in the above integral
		virtual double SubIntegral(const double dA, const double dB) const = 0;
        
	};
}
#endif
