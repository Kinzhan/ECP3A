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

void CapletPricingInterface(const double dMaturity, const double dTenor, const double dStrike, const std::size_t iNPaths);

void CapletPricingInterface(const double dMaturity, const double dTenor, const double dStrike, const std::size_t iNPaths)
{
    //  Initialization of LGM parameters 
    std::pair<std::vector<double>, std::vector<double> > dSigma;
    std::vector<std::pair<double, double> > dInitialYC;
    for (std::size_t i = 0 ; i < 4 * dMaturity + dTenor ; ++i)
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
    dSimulationTenors.push_back(dMaturity /*- dTenor*/);
    
    clock_t start = clock();
    sLGM.Simulate(iNPaths, dSimulationTenors, sSimulationData, /*Step by Step MC ?*/ /*iStepbyStepMC*/ true);
    std::cout << "Simulation Time : "<< (double)(clock()-start)/CLOCKS_PER_SEC <<" sec" <<std::endl;
    
    //sSimulationData.PrintInFile("/Users/alexhum49/Desktop/MyfileRiskNeutral.txt", false);
    //sSimulationData.PrintInFile("/Users/Kinz/Desktop/MyFilecopy.txt", false);
    
    // PRINT FACTORS AFTER CHANGE OF PROBA
    Finance::SimulationData sSimulationDataTForward;
    sLGM.ChangeOfProbability(dMaturity + dTenor, sSimulationData, sSimulationDataTForward);    
    //sSimulationDataTForward.PrintInFile("/Users/alexhum49/Desktop/MyfileTForward.txt", false);
    //sSimulationDataTForward.PrintInFile("/Users/Kinz/Desktop/MyfileRiskNeutral.txt", false);
    
    //Products::ProductsLGM sCaplet(sLGM);
    //std::vector<double> dPrice = sCaplet.Caplet(dMaturity - dTenor, dTenor, dStrike, sSimulationDataTForward);
    
    /*
     std::cout << "Do you want to test normality of simulated variables ? (0/1)"<<std::endl;
     std::size_t iTestOfNormality = 0;
     std::cin >> iTestOfNormality;
     
     if (iTestOfNormality)
     {
     std::cout << "Test of normality" << std::endl;
     
     std::cout << "Real variance is :" << std::endl;
     for (std::size_t iTenor = 0 ; iTenor < dSimulationTenors.size() ; ++iTenor)
     {
     std::cout << sSigmaTS.GetValues()[0] * sSigmaTS.GetValues()[0] * MathFunctions::Beta_OU(-2.0 * dLambda, dSimulationTenors[iTenor])<<std::endl;
     }
     
     std::cout << "Estimated variance is :" << std::endl;
     for (std::size_t iTenor = 0 ; iTenor < dSimulationTenors.size() ; ++iTenor)
     {
     double dMean = 0, dVariance = 0;
     for (std::size_t iPath = 0 ; iPath < 2 * iNPaths ; ++iPath)
     {
     double dLocalValue = sSimulationData.GetData().second[iTenor][iPath][0];
     dMean += dLocalValue;
     dVariance += dLocalValue * dLocalValue;
     }
     dMean /= (2.0 * iNPaths);
     dVariance /= (2.0 * iNPaths);
     dVariance -= dMean * dMean;
     std::cout << dVariance << std::endl;
     }
     }
     
     //Change of probability T-forward probability (T = maturity)
     std::cout << "Do you want to test the change of probability :" << std::endl;
     std::size_t iTestOfChangeOfProbability = 0;
     std::cin >> iTestOfChangeOfProbability;
     
     if (iTestOfChangeOfProbability)
     {
     std::cout << "Test of change of probability" << std::endl;
     Finance::SimulationData sSimulationDataTForward;
     sLGM.ChangeOfProbability(dMaturity, sSimulationData, sSimulationDataTForward);
     
     //sSimulationDataTForward.PrintInFile("/Users/alexhum49/Desktop/MyfileTForward.txt", false);
     sSimulationDataTForward.PrintInFile("/Users/Kinz/Desktop/MyfileRiskNeutral.txt", false);
     }
     */
    
    //COMPUTATION ON THE PRICE
    Products::ProductsLGM sProductLGM(sLGM);
    //double dPrice = 0.0, dStdDevPrice = 0.0;
    //std::vector<double> dPayoff = sProductLGM.Caplet(dMaturity - dTenor, dMaturity, dStrike, sSimulationDataTForward);
    std::vector<double> dPayoff = sProductLGM.Caplet(dMaturity, dMaturity + dTenor, dStrike, sSimulationDataTForward);
    
    double dMCPrice = 0;
    std::size_t iPath = 100;
    std::vector<double> dMeanPrice;
    for (std::size_t iLoop = 0 ; iLoop < iNPaths ; ++iLoop)
    {
        dMCPrice += dPayoff[iLoop];
        if (iLoop % iPath == 0 && iLoop != 0)
        {
            dMeanPrice.push_back(dMCPrice * exp(-sInitialYC.YC(dMaturity + dTenor) * dMaturity + dTenor) / iLoop);
        }
    }
    Utilities::PrintInFile sPrint("/Users/alexhum49/Desktop/TextCaplet.txt", false, 6);
    sPrint.PrintDataInFile(dMeanPrice);
    std::cout << "Print in file : done ! "<< std::endl;
    
    std::cout << "Final PV : " << dMeanPrice.back() << std::endl;
    
    /*for (std::size_t iPath = 0 ; iPath < iNPaths ; ++iPath)
     {
     //  Payment Date discount factor
     //std::vector<double> dDFs = sProductLGM.RiskNeutralDiscountFactor(iPath, sSimulationData);
     //double dDF = dDFs.back();
     double dDF = 1.0;
     dPrice += dDF * dPayoff[iPath];
     dStdDevPrice += dDF * dPayoff[iPath] * dDF * dPayoff[iPath];
     }
     dPrice /= iNPaths;
     dStdDevPrice /= iNPaths;
     dStdDevPrice -= dPrice * dPrice;
     
     std::cout << "Price = " << dPrice << std::endl;
     std::cout << "StdDev = " << dStdDevPrice << std::endl;*/
    
    //  Black-Scholes Price 
    if (!sSigmaTS.IsTermStructure())
    {
        double dVolSquareModel = (MathFunctions::Beta_OU(dLambda, dMaturity + dTenor) - MathFunctions::Beta_OU(dLambda, dMaturity)) * (MathFunctions::Beta_OU(dLambda, dMaturity + dTenor) - MathFunctions::Beta_OU(dLambda, dMaturity)) * dSigma.second[0] * dSigma.second[0] * (exp(2.0 * dLambda * (dMaturity)) - 1.0) / (2.0 * dLambda);
        std::cout << "Vol : " << sqrt(dVolSquareModel) << std::endl;
        double dForward = exp(sInitialYC.YC(dMaturity + dTenor) * dMaturity + dTenor) / exp(sInitialYC.YC(dMaturity) * (dMaturity));
        std::cout << "Black-Scholes Price : " << (1.0 + dTenor * dStrike) * MathFunctions::BlackScholes(dForward, 1.0 / (1.0 + dTenor * dStrike), sqrt(dVolSquareModel * (dMaturity)), -1) << std::endl;
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
        
        //double dVariable = 0;
        //std::cout<<"Enter value : ";
        //std::cin>>dVariable;
        
        //std::cout << "Result Value is : "<< sInterp1d.Interp1D(dVariable) <<std::endl;
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
        
        //std::cout << "Currency : " << sYieldCurve.GetCurrency() << std::endl;
        //std::cout << "Name : " << sYieldCurve.GetName() << std::endl;
        
        double dt = 0;
        std::cout<<"Enter value : ";
        std::cin>>dt;
        std::cout << "Yield Curve value is : " << sYieldCurve.YC(dt);

    }
    else if (iChoice == 5)
    {
        //  Test for Simulation Data
        //std::string cFile = "/Users/alexhum49/Desktop/MyFile.txt";
		std::string cFile = "/Users/Kinz/Desktop/MyFile.txt";
        Finance::SimulationData sSimulatedData;
        sSimulatedData.ReadFromFile(cFile);
        
        //sSimulatedData.PrintInFile("/Users/alexhum49/Desktop/MyFilecopy.txt", false);
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
        double dMaturity = 2.0, dTenor = 0.5, dStrike = 0.0;
        //std::size_t iStepbyStepMC = false;
        
        //  Input some variables
        /*std::cout << "Enter the number of paths : ";
        std::cin >> iNPaths;
        std::cout << "Enter the maturity of the caplet (in years): ";
        std::cin >> dMaturity;
        Utilities::require(dMaturity > 0 , "Maturity is negative");
		std::cout << "Enter Strike of caplet : ";
        std::cin >> dStrike;
        Utilities::require(dStrike > 0, "Strike is negative");
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
        sLGM.Simulate(iNPaths, dSimulationTenors, sSimulationData, /*Step by Step MC ?*/ /*iStepbyStepMC*/ true);
        std::cout << "Simulation Time : "<< (double)(clock()-start)/CLOCKS_PER_SEC <<" sec" <<std::endl;
        
        //  change of probability to T forward neutral
        Finance::SimulationData sSimulationDataTForward;
        sLGM.ChangeOfProbability(dT1, sSimulationData, sSimulationDataTForward);
        //sSimulationData.PrintInFile("/Users/alexhum49/Desktop/SimulationData.txt", false);
        //sSimulationDataTForward.PrintInFile("/Users/alexhum49/Desktop/SimulationDataTForward.txt", false);
        
        //  compute forward bond prices
        std::vector<double> dForwardBondPrice;
        Finance::SimulationData::Cube sDataTForwardCube = sSimulationDataTForward.GetData().second;
        Finance::SimulationData::Cube sDataCube = sSimulationData.GetData().second;
        std::vector<double> dFactorRiskNeutral, dFactorTForwardNeutral;
        
        std::size_t iDate = 0;
        
        Stats::Statistics sStats;
        for (std::size_t i = 0 ; i < iNPaths ; ++i)
        {
            dFactorRiskNeutral.push_back(sDataCube[iDate][i][0]);
            dFactorTForwardNeutral.push_back(sDataTForwardCube[iDate][i][0]);
        }
        
        std::vector<double> dDF;
        for (std::size_t iPath = 0; iPath < iNPaths ; ++iPath)
        {
            double dFactorTNeutral = sDataTForwardCube[iDate][iPath][0];
            dDF.push_back(sLGM.BondPrice(dT1, dT2, dFactorTNeutral));
        }
        Utilities::PrintInFile sPrintInFile("/Users/alexhum49/Desktop/DF.txt", false, 7);
        sPrintInFile.PrintDataInFile(dDF);
        
        std::cout << "Forward bond price by simulation : " << sStats.Mean(dDF) << std::endl;
        std::cout << "Bond Price value : " << exp(-sInitialYC.YC(dT2) * dT2) / exp(-sInitialYC.YC(dT1) * dT1) << std::endl;
        double dMCPrice = 0;
        std::size_t iPath = 100;
        std::vector<double> dMeanPrice;
        for (std::size_t iLoop = 0 ; iLoop < iNPaths ; ++iLoop)
        {
            dMCPrice += dDF[iLoop];
            if (iLoop % iPath == 0 && iLoop != 0)
            {
                dMeanPrice.push_back(dMCPrice / (iLoop + 1));
            }
        }
        Utilities::PrintInFile sPrint("/Users/alexhum49/Desktop/TestMartingality.txt", false, 6);
        sPrint.PrintDataInFile(dMeanPrice);
        std::cout << "Print in file : done ! "<< std::endl;
        
        //  Loop Test ? 
        /*std::size_t iLoopTest = 0;
        std::cout << "Loop test over maturities (0/1)? "<< std::endl;
        std::cin >> iLoopTest;
        if (iLoopTest == 1)
        {
            unsigned int iMaxMaturity = 30; // 30 year discount factor
            for (unsigned int iMaturity = 1 ; iMaturity < iMaxMaturity ; ++iMaturity)
            {
                double dT2 = static_cast<double>(iMaturity);
                //  Initialization of LGM parameters 
                std::vector<std::pair<double, double> > dInitialYC;
                for (std::size_t i = 0 ; i < 4 * dT2 ; ++i)
                {
                    dInitialYC.push_back(std::make_pair(0.25 * (i + 1), 0.03)); // YieldCurve = 3%
                }
                
                //  Initialization of classes
                Finance::YieldCurve sInitialYC("", "", dInitialYC, Utilities::Interp::LIN); // to change to SPLINE_CUBIC
                Processes::LinearGaussianMarkov sLGM(sInitialYC, dLambda, sSigmaTS);
                Finance::SimulationData sSimulationData;
                std::vector<double> dSimulationTenors;
                dSimulationTenors.push_back(dT1);
                
            sLGM.Simulate(iNPaths, dSimulationTenors, sSimulationData, true);
                
                //  change of probability to T forward neutral
                Finance::SimulationData sSimulationDataTForward;
                sLGM.ChangeOfProbability(dT2, sSimulationData, sSimulationDataTForward);    
                
                //  compute forward bond prices
                std::vector<double> dForwardBondPrice;
                Finance::SimulationData::Cube sDataTForwardCube = sSimulationDataTForward.GetData().second;
                std::size_t iDate = 0;
                std::vector<double> dDF;
                for (std::size_t iPath = 0; iPath < iNPaths ; ++iPath)
                {
                    double dFactor = sDataTForwardCube[iDate][iPath][0];
                    dDF.push_back(sLGM.BondPrice(dT1, dT2, dFactor, Processes::T_FORWARD_NEUTRAL));
                }
                
                Stats::Statistics sStats;
                std::cout << "Maturity : " << dT2 << std::endl;
                std::cout << "Forward bond price by simulation : " << sStats.Mean(dDF) << std::endl;
                std::cout << "Bond Price value : " << exp(-sInitialYC.YC(dT2) * dT2) / exp(-sInitialYC.YC(dT1) * dT1) << std::endl;;
                
            }
        }*/
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


