//
//  StochasticBasisSpread.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 09/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_StochasticBasisSpread_h
#define Seminaire_StochasticBasisSpread_h

namespace Processes {
    //  This class model a stochastic basis spread
    class StochasticBasisSpread
    {
    protected:
        double dSigma_, dLambda_;
    public:
        StochasticBasisSpread();
        virtual ~StochasticBasisSpread();
        
        virtual double Sigma_B(double t, double T) const;
    };
}

#endif
