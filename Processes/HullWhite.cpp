//
//  HullWhite.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 19/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "HullWhite.h"
#include "Require.h"
#include "Gaussian.h"
#include "MathFunctions.h"

namespace Processes {
    
    LinearGaussianMarkov::LinearGaussianMarkov()
    {}
    
    LinearGaussianMarkov::LinearGaussianMarkov(const Finance::YieldCurve & sInitialYieldCurve, const double dLambda, const Finance::TermStructure<double, double> & dSigma) : dLambda_(dLambda)
    {
        sInitialYieldCurve_ = sInitialYieldCurve;
        dSigma_ = dSigma;
    }
    
    LinearGaussianMarkov::~LinearGaussianMarkov()
    {}
    
    void LinearGaussianMarkov::Simulate(const std::size_t iNRealisations,
                                        const std::vector<double> & dSimulationTenors,
                                        Finance::SimulationData & sSimulationData,
                                        bool bIsStepByStepMC) const
    {
        Utilities::require(!dSimulationTenors.empty(), "Simulation Times is empty");
        Utilities::require(iNRealisations > 0, "Number of paths has to be positive");
        
        
        sSimulationData.SetDates(dSimulationTenors);
        
        std::vector<double> dSimulationTenorsCopy;
        dSimulationTenorsCopy.insert(dSimulationTenorsCopy.begin(), 0.0);
        dSimulationTenorsCopy.insert(dSimulationTenorsCopy.end(), dSimulationTenors.begin(), dSimulationTenors.end());
        
        //  In this function, we will simulate the factor \int_{s}^{t} a(u) dW^Q_u for s,t in the simulation tenors vector
        std::size_t iNTenors = dSimulationTenorsCopy.size();
        RandomNumbers::Gaussian1D sGaussian(0.0, 1.0, iNRealisations * iNTenors, 0);
        sGaussian.GenerateGaussian();
        
        if (bIsStepByStepMC)
        {
            //  Step by Step Monte Carlo
            for (std::size_t iSimulationTenor = 0 ; iSimulationTenor < iNTenors - 1 ; ++iSimulationTenor)
            {
                double dVariance = 0.0;
                if (!dSigma_.IsTermStructure())
                {
                    dVariance= dSigma_.GetValues()[0] * dSigma_.GetValues()[0] * MathFunctions::Beta_OU(-2.0 * dLambda_, dSimulationTenorsCopy[iSimulationTenor + 1]);
                }
                else
                {
                    //  Computation of the variance with term structure
                    /*for (std::size_t iTS = 0 ; dSigma_.GetVariables()[iTS] <= dSimulationTenorsCopy[iSimulationTenor] ; ++iTS)
                    {
                        dVariance += dSigma_.GetValues()[iTS] * dSigma_.GetValues()[iTS] * (MathFunctions::Beta_OU(-2.0 * dLambda_, std::min(dSimulationTenorsCopy[iSimulationTenor + 1],dSigma_.GetVariables()[iTS + 1])) - MathFunctions::Beta_OU(-2.0 * dLambda_, std::max(dSimulationTenorsCopy[iSimulationTenor],dSigma_.GetVariables()[iTS])));
                    }*/
                    std::cout<<"Not yet implemented"<<std::endl;
                }
                
                for (std::size_t iPath = 0 ; iPath < iNRealisations ; ++iPath)
                {
                    double dCurrentValue = sGaussian.GetRealisations()[iPath + iNRealisations * iSimulationTenor] * sqrt(dVariance);
                    //sSimulationData.Put(iSimulationTenor, iPath , std::vector<double>(1, sGaussian.GetRealisations()[iPath + iNRealisations * iSimulationTenor] * sqrt(dVariance) + (iSimulationTenor == 0 ? 0.0 : sSimulationData.Get(iSimulationTenor - 1, iPath)[0])));
                    sSimulationData.Put(iSimulationTenor, iPath, std::vector<double>(1, dCurrentValue));
                    if (1) // Antithetic variables
                    {
                        sSimulationData.Put(iSimulationTenor, iNRealisations + iPath, std::vector<double>(1, -dCurrentValue));
                    }
                }
            }
            //std::cout<<"Not yet implemented. To be checked"<<std::endl;
        }
        else
        {
            //  Path by Path Monte Carlo
            
            //  Compute the standard deviation
            std::vector<double> dStdDev(iNTenors - 1, 0.0);
            if (!dSigma_.IsTermStructure())
            {
                /*for (std::size_t iSimulationTenor = 0 ; iSimulationTenor < iNTenors - 1 ; ++iSimulationTenor)
                {
                    dStdDev[iSimulationTenor] = dSigma_.GetValues()[0] * sqrt(MathFunctions::Beta_OU(-2.0 * dLambda_, dSimulationTenorsCopy[iSimulationTenor + 1]) - MathFunctions::Beta_OU(2 * dLambda_, dSimulationTenorsCopy[iSimulationTenor]));
                }*/
                std::cout << "Not yet implemented"<<std::endl;
            }
            else
            {
                /*for (std::size_t iSimulationTenor = 0 ; iSimulationTenor < iNTenors - 1; ++ iSimulationTenor)
                {
                    for (std::size_t iTS = 0 ; dSigma_.GetVariables()[iTS] <= dSimulationTenorsCopy[iSimulationTenor] ; ++iTS)
                    {
                        dStdDev[iSimulationTenor] += dSigma_.GetValues()[iTS] * dSigma_.GetValues()[iTS] * (MathFunctions::Beta_OU(-2.0 * dLambda_, std::min(dSimulationTenorsCopy[iSimulationTenor + 1],dSigma_.GetVariables()[iTS + 1])) - MathFunctions::Beta_OU(-2.0 * dLambda_, std::max(dSimulationTenorsCopy[iSimulationTenor],dSigma_.GetVariables()[iTS])));
                    }
                    if (dStdDev[iSimulationTenor] < 0.0)
                    {
                        std::cout<< "Try to compute square root of negative number for simulation tenor " << iSimulationTenor << std::endl;
                    }
                    else
                    {
                        dStdDev[iSimulationTenor] = sqrt(dStdDev[iSimulationTenor]);
                    }
                }*/
                std::cout<<"Not yet implemented"<<std::endl;
            }
            //  Begin the simulation
            for (std::size_t iPath = 0 ; iPath < iNRealisations ; ++iPath)
            {
                //  We simulate the wanted factor
                double dCurrentValue = 0.0;
                for (std::size_t iSimulationTenor = 0 ; iSimulationTenor < iNTenors - 1; ++iSimulationTenor)
                {
                    dCurrentValue = sGaussian.GetRealisations()[iPath + iSimulationTenor * iNRealisations] * dStdDev[iSimulationTenor];
                    sSimulationData.Put(iSimulationTenor, iPath, std::vector<double>(1,dCurrentValue));
                    if (1) // Antithetic variables
                    {
                        sSimulationData.Put(iSimulationTenor, iPath + iNRealisations, std::vector<double>(1,- dCurrentValue));
                    }
                }
                
            }
        }
    }
    
    void LinearGaussianMarkov::ChangeOfProbability(const double dT, 
                                                   const Finance::SimulationData &sSimulationDataRiskNeutral, 
                                                   Finance::SimulationData &sSimulationDataTForward) const
    {
        //  Create a new simulated data object w.r.t. the T-forward neutral probability
        //  The results factors sSimulationDataTForward are martingales under the T-forward neutral probability and the input factors sSimulationDataRiskNeutral are martingales under the risk neutral probability
        
        sSimulationDataTForward.SetDates(sSimulationDataRiskNeutral.GetDateList());
        std::vector<long> lDates = sSimulationDataRiskNeutral.GetDateList();
        for (std::size_t iDate = 0 ; iDate < lDates.size() ; ++iDate)
        {
            double dDate = lDates[iDate] / 365.0;
            double dDrift = DriftChangeOfProbability(dDate, dT);
            
            for (std::size_t iPath =  0 ; iPath < sSimulationDataRiskNeutral.GetData().second[iDate].size() ; ++iPath)
            {
                std::vector<double> dTForwardValues;
                for (std::size_t iVar = 0 ; iVar < sSimulationDataRiskNeutral.GetData().second[iDate][iPath].size() ; ++iVar)
                {
                    dTForwardValues.push_back(sSimulationDataRiskNeutral.GetData().second[iDate][iPath][iVar] + dDrift);
                }
                sSimulationDataTForward.Put(iDate, iPath, dTForwardValues);
            }
        }
    }
    
    double LinearGaussianMarkov::DriftChangeOfProbability(const double dt, const double dT) const
    {
        //  Compute the following integral \int_{0}^{t} ((\beta(T) - \beta(s)) * a(s)) ds where \beta(t) = \int_{0}^{t} b(s) ds
        if (dSigma_.IsTermStructure())
        {
            std::size_t iNSigmaTS = dSigma_.GetValues().size();
            double dResult = 0.0;
            for (std::size_t iTS = 0 ; iTS < iNSigmaTS ; ++iTS)
            {
                double dSigma = dSigma_.GetValues()[iTS];
                double dt2 = dSigma_.GetVariables()[iTS], dt1 = iTS == 0 ? 0 : dSigma_.GetVariables()[iTS - 1];
                dResult += dSigma * dSigma * MathFunctions::Beta_OU(-dLambda_, dT) * MathFunctions::SumExp(dLambda_, dt1, dt2) - 1.0 / dLambda_ * (MathFunctions::SumExp(2.0 * dLambda_, dt1, dt2) - MathFunctions::SumExp(dLambda_, dt1, dt2));
            }
            return dResult;
        }
        else
        {
            //  No term-structure in the volatiliy
            double dSigma = dSigma_.GetValues()[0];
            return dSigma * dSigma * (MathFunctions::Beta_OU(dLambda_, dT) * MathFunctions::Beta_OU(-2.0 * dLambda_, dt) - 1.0 / dLambda_ * (MathFunctions::Beta_OU(-2.0 * dLambda_, dt) - MathFunctions::Beta_OU(-dLambda_, dt)));
        }
    }
    
    double LinearGaussianMarkov::A(double t) const
    {
        return dSigma_.Interpolate(t) * exp(dLambda_ * t);
    }
    
    double LinearGaussianMarkov::B(double t) const
    {
        return exp(-dLambda_ * t);
    }
    
    double LinearGaussianMarkov::BondPrice(const double dt, const double dT, const double dX, const SimulationProbability eProbability) const
    {
        // What is dX here ?
        // B(0,T) / B(0,t) * exp(DeterministPart)	
        double aux1 = exp(-sInitialYieldCurve_.YC(dT) * dT) / exp(-sInitialYieldCurve_.YC(dt) * dt) * exp(MathFunctions::DeterministPartZCBHW1F(dLambda_, dt, dT, dSigma_)) ;
        // sum(b(u)du, u=t..T)
        double aux2 = MathFunctions::SumExp(dLambda_, dt, dT) ;
        
        return aux1 * exp(- dX * aux2 + (eProbability == T_FORWARD_NEUTRAL ? DriftChangeOfProbability(dt, dT) : 0.0));
    }
}