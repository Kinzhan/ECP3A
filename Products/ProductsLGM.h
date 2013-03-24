//
//  ProductsLGM.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 31/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_ProductsLGM_h
#define Seminaire_ProductsLGM_h

#include <iostream>
#include "HullWhite.h"

namespace Products {
    class ProductsLGM : public Processes::LinearGaussianMarkov
    {
    public:
        
        ProductsLGM(const Processes::LinearGaussianMarkov & sLGMProcess, double dEpsilonMaturity = 0.001);
        virtual ~ProductsLGM();
        
        virtual std::vector<double> RiskNeutralDiscountFactor(std::size_t iPath, const Finance::SimulationData & sSimulationData) const;
        
        virtual std::vector<double> Caplet(double dStart, double dEnd, double dPay, double dStrike, const Finance::SimulationData & sSimulationData, const Processes::CurveName & eCurveName, double dQA = 1) const;
    private:
        //  about 1 day
        double dEpsilonMaturity_;
    };
}

#endif
