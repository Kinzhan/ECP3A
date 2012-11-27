//
//  ProductsLGM.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 31/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "ProductsLGM.h"

namespace Products {
    
    ProductsLGM::ProductsLGM(const Processes::LinearGaussianMarkov & sLGMProcess, const double dEpsilonMaturity) : dEpsilonMaturity_(dEpsilonMaturity)
    {
        // Added by Emile, 21Nov12
		sInitialYieldCurve_ = sLGMProcess.GetYieldCurve();
		//sInitialYieldCurve_ = sLGMProcess.sInitialYieldCurve_;
        dLambda_ = sLGMProcess.GetLambda();
        dSigma_ = sLGMProcess.GetSigma();
    }
    
    std::vector<double> ProductsLGM::Caplet(const double dStart, const double dEnd, const double dStrike, const Finance::SimulationData & sSimulationData) const
    {
        std::size_t iWhere = Utilities::FindInVector(sSimulationData.GetDateList(), static_cast<long>(dStart * 365)) + 1;
        
        if (iWhere != sSimulationData.GetDateList().size())
        {
            std::vector<std::vector<double> > dMatrixEndFactor = sSimulationData.GetData().second[iWhere];
        
            std::size_t iNPaths = dMatrixEndFactor.size();
            std::vector<double> dResults (iNPaths);
        
            for (std::size_t iPath = 0 ; iPath < iNPaths ; ++iPath)
            {
                //  Only one factor which is simulated for now
                //  Fixing of the libor at start date of the period
                dResults[iPath] = std::max(Libor(dStart, dStart, dEnd, dMatrixEndFactor[iPath][0], Processes::RISK_NEUTRAL) - dStrike, 0.0);
            }
        
            return dResults;
        }
        else
        {
            std::cout<< "Start date not found in simulation" << std::endl;
            return std::vector<double>(0,0.0);
        }
    }
    
    std::vector<double> ProductsLGM::RiskNeutralDiscountFactor(const std::size_t iPath, const Finance::SimulationData &sSimulationData) const
    {
        /*std::vector<double> dResults, dFactors, dShortRate, dSimulationDates;
        std::vector<long> lDates = sSimulationData.GetDateList();
        Finance::SimulationData::Data sData = sSimulationData.GetData();
        
        std::size_t iNDates = lDates.size();
        dFactors.resize(iNDates);
        dShortRate.resize(iNDates);
        dSimulationDates.resize(iNDates);
        dResults.resize(iNDates);
       
        for (std::size_t iDate = 0 ; iDate < iNDates ; ++iDate)
        {
            dFactors[iDate] = sData.second[iDate][iPath][0];
            dSimulationDates[iDate] = lDates[iDate] / 365.0;
        
            dShortRate[iDate] = - (log(BondPrice(dSimulationDates[iDate], dSimulationDates[iDate] + dEpsilonMaturity_, dFactors[iDate])) - log(BondPrice(dSimulationDates[iDate], dSimulationDates[iDate] - dEpsilonMaturity_, dFactors[iDate]))) / (2.0 * dEpsilonMaturity_);
        }
        
        dResults[0] = exp(-dSimulationDates[0] * dShortRate[0]);
        for (std::size_t iDate = 1 ; iDate < iNDates ; ++iDate)
        {
            dResults[iDate] = dResults[iDate - 1] * exp((dShortRate[iDate] + dShortRate[iDate - 1]) * 0.5 * (dSimulationDates[iDate - 1] - dSimulationDates[iDate]));
        }
        
        return dResults;*/
        return std::vector<double>(0,0.0);
    }
    
}