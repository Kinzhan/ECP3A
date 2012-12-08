//
//  main.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 14/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <time.h>
#include "Uniform.h"
#include "Gaussian.h"
#include "Test.h"
#include "MathFunctions.h"

//  Testing of Hull White model
#include "HullWhite.h"
#include "InterExtrapolation.h"
#include "VectorUtilities.h"
#include "YieldCurve.h"
#include "SimulationData.h"
#include "ProductsLGM.h"

#include "Statistics.h"
#include "PrintInFile.h"

void CapletPricingInterface(const double dMaturity, const double dTenor, const double dStrike, std::size_t iNPaths);

void CapletPricingInterface(const double dMaturity, const double dTenor, const double dStrike, std::size_t iNPaths)
{
    /*
     dMaturity : maturity of caplet (should be positive and in years)
     dTenor : Tenor of caplet (should be positive and in years)
     dStrike : Strike of caplet (not-necessarly positive
     iNPaths : number of paths for Monte-Carlo pricing
     */
    
    //  Initialization of LGM parameters 
    std::pair<std::vector<double>, std::vector<double> > dSigma;
    std::vector<std::pair<double, double> > dInitialYC;
    for (std::size_t i = 0 ; i < 4 * (dMaturity + dTenor) ; ++i)
    {
        dInitialYC.push_back(std::make_pair(0.25 * (i + 1), 0.03)); // YieldCurve = 3%
    }
    dSigma.first.push_back(0);
    dSigma.second.push_back(0.01); // volatility = 1%
    double dLambda = 0.05; // Mean reversion = 5%
    
    //  Initialization of classes
    Finance::TermStructure<double, double> sSigmaTS(dSigma.first, dSigma.second);
    Finance::YieldCurve sInitialYC("", "", dInitialYC, Utilities::Interp::LIN); // to change to SPLINE_CUBIC
    Processes::LinearGaussianMarkov sLGM(sInitialYC, dLambda, sSigmaTS);
    Finance::SimulationData sSimulationData;
    std::vector<double> dSimulationTenors;
    dSimulationTenors.push_back(dMaturity);
    
    clock_t start = clock();
    sLGM.Simulate(iNPaths, dSimulationTenors, sSimulationData, /*Step by Step MC ?*/ true);
    std::cout << "Simulation Time : "<< (double)(clock()-start)/CLOCKS_PER_SEC <<" sec" <<std::endl;
    start = clock();
    
    // PRINT FACTORS AFTER CHANGE OF PROBA
    Finance::SimulationData sSimulationDataTForward;
    sLGM.ChangeOfProbability(dMaturity + dTenor, sSimulationData, sSimulationDataTForward); 
    std::cout << "Change of probability time : " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << std::endl;
    
    //COMPUTATION ON THE PRICE
    Products::ProductsLGM sProductLGM(sLGM);
    std::vector<double> dPayoff = sProductLGM.Caplet(dMaturity, dMaturity + dTenor, dStrike, sSimulationDataTForward);
    
    double dMCPrice = 0;
    std::size_t iPath = 100;
    std::vector<double> dMeanPrice;
    iNPaths = dPayoff.size();
    for (std::size_t iLoop = 0 ; iLoop < iNPaths ; ++iLoop)
    {
        dMCPrice += dPayoff[iLoop];
        if (iLoop % iPath == 0 && iLoop != 0)
        {
            dMeanPrice.push_back(dMCPrice * exp(-sInitialYC.YC(dMaturity + dTenor) * (dMaturity + dTenor)) / iLoop);
        }
    }
    Utilities::PrintInFile sPrint("/Users/alexhum49/Desktop/TextCaplet.txt", false, 6);
    sPrint.PrintDataInFile(dMeanPrice);
    std::cout << "Print in file : done ! "<< std::endl;
    
    std::cout << "Final PV : " << dMeanPrice.back() << std::endl;
    
    //  Black-Scholes Price 
    if (!sSigmaTS.IsTermStructure())
    {
        //  Sigma of HW1F model
        double dSigmaVol = dSigma.second[0];
        
        //  Integrated variance of fwd Zero-coupon bond
        double dVolSquareModel = (MathFunctions::Beta_OU(dLambda, dMaturity + dTenor) - MathFunctions::Beta_OU(dLambda, dMaturity)) * (MathFunctions::Beta_OU(dLambda, dMaturity + dTenor) - MathFunctions::Beta_OU(dLambda, dMaturity)) * dSigmaVol * dSigmaVol * (exp(2.0 * dLambda * (dMaturity)) - 1.0) / (2.0 * dLambda);
        
        //  B(t,T,T+\delta) = B(t,T+\delta) / B(t,T) --> forward discount factor
        double dForward = exp(-sInitialYC.YC(dMaturity + dTenor) * (dMaturity + dTenor)) / exp(-sInitialYC.YC(dMaturity) * (dMaturity));
        
        //  Output Black-Scholes result
        double dDFPaymentDate = exp(-sInitialYC.YC(dMaturity + dTenor) * (dMaturity + dTenor));
        std::cout << "Black-Scholes Price : " << dDFPaymentDate * (1.0 + dTenor * dStrike) * MathFunctions::BlackScholes(dForward, 1.0 / (1.0 + dTenor * dStrike), sqrt(dVolSquareModel), Finance::PUT) << std::endl;
    }
    
    std::cout<<"Total Time elapsed : " << (double)(clock()-start)/CLOCKS_PER_SEC <<" sec"<< std::endl;
}

int main()
{
    std::vector<double> dRealisations;
    
    std::cout << "Hello" << std::endl;
    std::size_t iNRealisations = 10, iAntitheticVariables = 1, iChoice;
    
    std::cout<< "Which variables do you want to simulate " << std::endl;
    std::cout << "1-  Uniform" << std::endl;
    std::cout << "2-  Gaussian" << std::endl;
    std::cout << "3-  Interpolation" << std::endl;
    std::cout << "4-  YieldCurve" << std::endl;
    std::cout << "5-  SimulationData" << std::endl;
    std::cout << "6-  AccCumNorm" << std::endl;
    std::cout << "75- HullWhite" << std::endl;
    std::cout << "76- Test" << std::endl;
    std::cout << "77- Martingality of Bond Price" << std::endl;
    std::cin >> iChoice;
    if (iChoice == 1 || iChoice == 2)
    {
        std::cout<< "Enter the number of realisations : ";
        std::cin >> iNRealisations;
        std::cout<< "Do you want to use antithetic variables ?";
        std::cin >> iAntitheticVariables;

        if (iChoice == 1)
        {
            //  Test for uniform variables
            RandomNumbers::Uniform sUnif(0.0,1.0,iNRealisations,static_cast<int>(iAntitheticVariables));
            sUnif.GenerateUniform();
            dRealisations = sUnif.GetRealisations();
        }
        else if (iChoice == 2)
        {
            std::cout<<"Enter standard deviation :";
            double dStdDev = 1.0;
            std::cin>>dStdDev;
            // Test for gaussian variables
            RandomNumbers::Gaussian1D sGaussian(0.0, dStdDev, iNRealisations, static_cast<int>(iAntitheticVariables));
            sGaussian.GenerateGaussian();
            dRealisations = sGaussian.GetRealisations();
        }
    }
    else if (iChoice == 3)
    {
        //  Inter-Extrapolation test
        std::vector<std::pair<double, double> > dVariablesAndValues;
        dVariablesAndValues.push_back(std::make_pair(0.0, 0.0));
        dVariablesAndValues.push_back(std::make_pair(1.0, 1.0));
        dVariablesAndValues.push_back(std::make_pair(4.0, 2.0));
        dVariablesAndValues.push_back(std::make_pair(9.0, 3.0));
        dVariablesAndValues.push_back(std::make_pair(16.0, 4.0));
        
        std::pair<std::vector<double>, std::vector<double> > dVariablesAndValues0 = Utilities::GetPairOfVectorFromVectorOfPair(dVariablesAndValues);
        Utilities::Interp::InterExtrapolation1D sInterp1d(dVariablesAndValues0.first, dVariablesAndValues0.second, Utilities::Interp::SPLINE_CUBIC);
        
        for (std::size_t i = 0 ; i < 150 ; ++i)
        {
            std::cout << i / 5.0 << ";" << sInterp1d.Interp1D(i/5.0) << std::endl;
        }
    }
    else if (iChoice == 4)
    {
        //  Test for YieldCurve
        std::vector<std::pair<double, double> > dVariablesAndValues;
        dVariablesAndValues.push_back(std::make_pair(0.0, 1.0));
        dVariablesAndValues.push_back(std::make_pair(1.0, 6.0));
        dVariablesAndValues.push_back(std::make_pair(3.5, 3.0));
        dVariablesAndValues.push_back(std::make_pair(5.0, 7.0));
        dVariablesAndValues.push_back(std::make_pair(9.0, 8.0));
        Finance::YieldCurve sYieldCurve("EUR", "EUROIS", dVariablesAndValues, Utilities::Interp::SPLINE_CUBIC);
        
        double dt = 0;
        std::cout<<"Enter value : ";
        std::cin>>dt;
        std::cout << "Yield Curve value is : " << sYieldCurve.YC(dt);

    }
    else if (iChoice == 5)
    {
        //  Test for Simulation Data
        std::string cFile = "/Users/Kinz/Desktop/MyFile.txt";
        Finance::SimulationData sSimulatedData;
        sSimulatedData.ReadFromFile(cFile);
        
        sSimulatedData.PrintInFile("/Users/Kinz/Desktop/MyFile.txt", false);
        
        std::cout<< "Well done it is working fine !";
    }
    else if (iChoice == 6)
    {
        std::vector<double> dAccCumNorm;
        for (std::size_t i = 0 ; i < 500 ; ++i)
        {
            dAccCumNorm.push_back(MathFunctions::AccCumNorm(-2.5 + i * 0.01));
        }
        Utilities::PrintInFile sPrint("/Users/alexhum49/Desktop/AccCumNorm.txt", false, 6);
        sPrint.PrintDataInFile(dAccCumNorm);
        
        std::cout << "Print in file : succeeded" << std::endl;
    }
    else if (iChoice == 75)
    {
        //  Test for Hull-White model
        
        std::cout << "Caplet Pricing by simulation" << std::endl;
        std::size_t iNPaths = 1000000;
        double dMaturity = 2.0, dTenor = 0.5, dStrike = 0.01;
        //std::size_t iStepbyStepMC = false;
        
        //  Input some variables
        std::cout << "Enter the number of paths : ";
        std::cin >> iNPaths;
        /*std::cout << "Enter the maturity of the caplet (in years): ";
        std::cin >> dMaturity;
        Utilities::require(dMaturity > 0 , "Maturity is negative");*/
		std::cout << "Enter Strike of caplet : ";
        std::cin >> dStrike;
        /*Utilities::require(dStrike > 0, "Strike is negative");
        std::cout << "Enter Tenor (in years) : ";
        std::cin >> dTenor;
        Utilities::require(dTenor > 0, "Tenor is negative");
        Utilities::require(dTenor < dMaturity, "Tenor is higher than caplet maturity");*/
        /*std::cout << "Step by Step MC : ";
        std::cin >> iStepbyStepMC;*/

        /*for (std::size_t iStrike = 0 ; iStrike < 11 ; ++iStrike)
        {
            std::cout << "Strike : " << iStrike * 0.01 << std::endl;
            CapletPricingInterface(dMaturity, dTenor, iStrike * 0.01, iNPaths);
        }*/
        CapletPricingInterface(dMaturity, dTenor, dStrike, iNPaths);
    }
    else if (iChoice == 76)
    {
        std::vector<double> dVect(1,0.0);
        Test::TestVectorDaughter sTestDaughter(dVect);
        sTestDaughter.DoSomething();
    }
    else if (iChoice == 77)
    {
        //  Martingality of Bond Price B(t,T1,T2)
        std::size_t iNPaths = 1000000;
        /*std::cout << "Enter the number of paths : "<< std::endl;
        std::cin >> iNPaths;*/
        double dT2 = 1,  dT1 = 0.75;
        std::cout << "Start Date of Forward discout factor : " << std::endl;
        std::cin >> dT1;
        std::cout << "Maturity of discount factor : "<< std::endl;
        std::cin >> dT2;
        Utilities::require(dT2 > dT1,"Maturity of discount factor is before start date");
        
        //  Initialization of LGM parameters 
        std::pair<std::vector<double>, std::vector<double> > dSigma;
        std::vector<std::pair<double, double> > dInitialYC;
        for (std::size_t i = 0 ; i < 4 * dT2 ; ++i)
        {
            dInitialYC.push_back(std::make_pair(0.25 * (i + 1), 0.03)); // YieldCurve = 3%
        }
        dSigma.first.push_back(0);
        dSigma.second.push_back(0.01); // volatility = 1%
        double dLambda = 0.05; // Mean reversion = 5%
        
        //  Initialization of classes
        Finance::TermStructure<double, double> sSigmaTS(dSigma.first, dSigma.second);
        Finance::YieldCurve sInitialYC("", "", dInitialYC, Utilities::Interp::LIN); // to change to SPLINE_CUBIC
        Processes::LinearGaussianMarkov sLGM(sInitialYC, dLambda, sSigmaTS);
        Finance::SimulationData sSimulationData;
        std::vector<double> dSimulationTenors;
        dSimulationTenors.push_back(dT1);
        
        clock_t start = clock();
        sLGM.Simulate(iNPaths, dSimulationTenors, sSimulationData, /*Step by Step MC ?*/ true);
        std::cout << "Simulation Time : "<< (double)(clock()-start)/CLOCKS_PER_SEC <<" sec" <<std::endl;
        
        //  change of probability to T forward neutral
        Finance::SimulationData sSimulationDataTForward;
        sLGM.ChangeOfProbability(dT1, sSimulationData, sSimulationDataTForward);
        
        //  compute forward bond prices
        std::vector<double> dForwardBondPrice;
        Finance::SimulationData::Cube sDataTForwardCube = sSimulationDataTForward.GetData().second;
        Finance::SimulationData::Cube sDataCube = sSimulationData.GetData().second;
        
        std::size_t iDate = 0;
        
        std::vector<double> dFactorRiskNeutral;
        std::size_t iNPaths0 = sDataCube[iDate].size();
        std::cout << iNPaths0 << std::endl;
        for (std::size_t iPath = 0 ; iPath < sDataCube[iDate].size() ; ++iPath)
        {
            dFactorRiskNeutral.push_back(sDataCube[iDate][iPath][0]);
        }
        Stats::Statistics sStats;
        std::vector<std::pair<double, std::size_t> > dEmpiricalDistribution = sStats.EmpiricalDistribution(dFactorRiskNeutral, 1000);
        Utilities::PrintInFile sEmpiricalDistributionFile("/Users/alexhum49/Desktop/FactorRiskNeutral.txt", false, 7);
        sEmpiricalDistributionFile.PrintDataInFile(dEmpiricalDistribution);
        
        std::vector<double> dDFT1FwdNeutral;
        for (std::size_t iPath = 0; iPath < iNPaths0 ; ++iPath)
        {
            double dFactorT1FwdNeutral = sDataTForwardCube[iDate][iPath][0];
            dDFT1FwdNeutral.push_back(sLGM.BondPrice(dT1, dT2, dFactorT1FwdNeutral));
        }
         
        std::cout << "Forward bond price by simulation (T1 Forward Neutral) : " << sStats.Mean(dDFT1FwdNeutral) << std::endl;

        std::cout << "Bond Price value : " << exp(-sInitialYC.YC(dT2) * dT2) / exp(-sInitialYC.YC(dT1) * dT1) << std::endl;
        double dMCPrice = 0;
        std::size_t iPath = 100;
        std::vector<double> dMeanPrice;
        for (std::size_t iLoop = 0 ; iLoop < iNPaths0 ; ++iLoop)
        {
            dMCPrice += dDFT1FwdNeutral[iLoop];
            if (iLoop % iPath == 0 && iLoop != 0)
            {
                dMeanPrice.push_back(dMCPrice / (iLoop + 1));
            }
        }
        Utilities::PrintInFile sPrint("/Users/alexhum49/Desktop/TestMartingality.txt", false, 6);
        sPrint.PrintDataInFile(dMeanPrice);
        std::cout << "Print in file : done ! "<< std::endl;
        
        //  End of test of martingality of Bond price
    }
    
    Stats::Statistics sStats;
    iNRealisations = dRealisations.size();
    if (iNRealisations > 0)
    {
        std::cout << "Mean : " << sStats.Mean(dRealisations) << std::endl;
    }
    
    std::cout << "GoodBye ! " << std::endl;
    return 0;
}


