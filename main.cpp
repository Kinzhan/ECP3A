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

#include "ProductsLGM.h"

#include "Statistics.h"
#include "PrintInFile.h"
#include "FXMultiCurve.h"
#include "Date.h"
#include "StochasticBasisSpread.h"
#include "Coverage.h"
#include "ForwardRate.h"
#include <stdlib.h>
#include "Annuity.h"

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
    //Finance::YieldCurve sInitialYC("", "", dInitialYC, Utilities::Interp::LIN); // to change to SPLINE_CUBIC
    Finance::YieldCurve sDiscountCurve, sSpreadCurve, sForwardCurve; 
    sDiscountCurve = 0.03;
    sSpreadCurve = 0;
    sForwardCurve = sDiscountCurve + sSpreadCurve;
    Processes::LinearGaussianMarkov sLGM(sDiscountCurve, sSpreadCurve, dLambda, sSigmaTS);
    Finance::SimulationData sSimulationData;
    std::vector<double> dSimulationTenors;
    dSimulationTenors.push_back(dMaturity);
    
    clock_t start = clock();
    sLGM.Simulate(iNPaths, dSimulationTenors, sSimulationData, /*Step by Step MC ?*/ true);
    std::cout << "Simulation Time : "<< (double)(clock()-start)/CLOCKS_PER_SEC <<" sec" <<std::endl;
    start = clock();
    
    // PRINT FACTORS AFTER CHANGE OF PROBA
    Finance::SimulationData sSimulationDataTForward;
    sLGM.ChangeOfProbability(dMaturity, sSimulationData, sSimulationDataTForward); 
    
    std::cout << "Change of probability time : " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << std::endl;
    
    //COMPUTATION ON THE PRICE
    Products::ProductsLGM sProductLGM(sLGM);
    Processes::CurveName eCurveName = Processes::FORWARD;
    std::vector<double> dPayoff = sProductLGM.Caplet(dMaturity, dMaturity + dTenor, dMaturity + dTenor, dStrike, sSimulationDataTForward, eCurveName);
    
    double dMCPrice = 0;
    std::size_t iPath = 100;
    std::vector<double> dMeanPrice;
    iNPaths = dPayoff.size();
    double dDFPaymentDate = exp(-sDiscountCurve.YC(dMaturity) * dMaturity);
    
    for (std::size_t iLoop = 0 ; iLoop < iNPaths ; ++iLoop)
    {
        dMCPrice += dPayoff[iLoop];
        if (iLoop % iPath == 0 && iLoop != 0)
        {
            dMeanPrice.push_back(dMCPrice * dDFPaymentDate / (iLoop + 1));
        }
    }
    Utilities::PrintInFile sPrint("/Users/alexhum49/Desktop/TextCaplet.txt", false, 6);
    sPrint.PrintDataInFile(dMeanPrice);
    std::cout << "Print in file : done ! "<< std::endl;
    std::cout << "Final PV : ";
    if (dMeanPrice.size())
    {
        std::cout << dMeanPrice.back() << std::endl;
    }
    else
    {
        std::cout << dMCPrice * dDFPaymentDate / iNPaths << std::endl;
    }
    
    //  Black-Scholes Price 
    if (!sSigmaTS.IsTermStructure())
    {
        //  Sigma of HW1F model
        double dSigmaVol = dSigma.second[0];
        
        //  Integrated variance of fwd Zero-coupon bond
        double dVolSquareModel = (MathFunctions::Beta_OU(dLambda, dMaturity + dTenor) - MathFunctions::Beta_OU(dLambda, dMaturity)) * (MathFunctions::Beta_OU(dLambda, dMaturity + dTenor) - MathFunctions::Beta_OU(dLambda, dMaturity)) * dSigmaVol * dSigmaVol * (exp(2.0 * dLambda * (dMaturity)) - 1.0) / (2.0 * dLambda);
        
        //  B(t,T,T+\delta) = B(t,T+\delta) / B(t,T) --> forward discount factor
        double  dForward = exp(-sForwardCurve.YC(dMaturity + dTenor) * (dMaturity + dTenor)) / exp(-sForwardCurve.YC(dMaturity) * (dMaturity)),
                dStrikeZC = 1.0 / (1.0 + dTenor * dStrike);
        
        //  Output Black-Scholes result
        std::cout << "Black-Scholes Price : " ;
        std::cout << dDFPaymentDate * 1. / dStrikeZC * MathFunctions::BlackScholes(dForward, dStrikeZC, sqrt(dVolSquareModel), Finance::PUT) << std::endl;
    }
    
    std::cout<<"Total Time elapsed : " << (double)(clock()-start)/CLOCKS_PER_SEC <<" sec"<< std::endl;
}



int main()
{
    //  Initialization of Today Date 
    std::time_t lToday;
    std::tm *stm;
    time(&lToday);
    
    stm = localtime(&lToday);
    static Utilities::Date::MyDate sTodayDate(*stm);
    
    //  End of Initialization of Today Date
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
    std::cout << "7-  Black-Scholes" << std::endl;
    std::cout << "8-  Date" << std::endl;
    std::cout << "9-  Coverage" << std::endl;
    std::cout << "10- Annuity" << std::endl;
    std::cout << "75- Caplet Pricer HW1F" << std::endl;
    std::cout << "76- Test" << std::endl;
    std::cout << "77- Martingality of Bond Price" << std::endl;
    std::cout << "78- Martingality of Forward Libor" << std::endl;
    std::cout << "79- Test of MergeTermStructure" << std::endl;
    std::cout << "80- Quanto Adjustment (Obsolete)" << std::endl;
    std::cout << "81- Stochastic Basis Spread Parameters" << std::endl;
    std::cout << "82- Quanto Adjustment (Libor)" << std::endl;
    std::cout << "83- Caplet Price with Stochastic Basis Spread" << std::endl;
    std::cout << "84- Weight Calculation" << std::endl;
	std::cout << "85- Quanto Adjustment (Swap)" << std::endl;
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
    else if (iChoice == 7)
    {
        //  Test of Black-Scholes formula
        double dForward = 1.0, dStrike = 1.0, dVolatility = 0.20, dMaturity = 1, dFlatYCValue = 0.03;
        
        std::cout << "Enter Forward : " << std::endl;
        std::cin >> dForward;
        std::cout << "Enter Strike : " << std::endl;
        std::cin >> dStrike;
        std::cout << "Enter Maturity : " << std::endl;
        std::cin >> dMaturity;
        std::cout << "Enter Volatility :" << std::endl;
        std::cin >> dVolatility;
        std::cout << "Enter Flat Yield Curve Value : " << std::endl;
        std::cin >> dFlatYCValue;
        
        std::vector<std::pair<double, double> > dYieldCurve;
        for (std::size_t i = 0 ; i < 4 * dMaturity ; ++i)
        {
            dYieldCurve.push_back(std::make_pair(0.25 * (i + 1), dFlatYCValue));        
        }
        //  Initialization of classes
        Finance::YieldCurve sInitialYC("", "", dYieldCurve, Utilities::Interp::LIN); // to change to SPLINE_CUBIC
        
        double dDFPaymentDate = exp(-sInitialYC.YC(dMaturity) * dMaturity);
        std::cout << "Black-Scholes Price : " << dDFPaymentDate * MathFunctions::BlackScholes(dForward / dDFPaymentDate, dStrike, dVolatility * sqrt(dMaturity), Finance::CALL) << std::endl;
        
        //  End of Black-Scholes test
    }
    else if (iChoice == 8)
    {
        //  Beginning of date test
        std::cout << Utilities::Date::GetDate(sTodayDate) << std::endl;
        sTodayDate.Add(-1, Utilities::Date::DAY);
        sTodayDate.Print();
        
        sTodayDate.Add(-1, Utilities::Date::WEEK);
        sTodayDate.Print();
        
        sTodayDate.Add(-1, Utilities::Date::MONTH);
        sTodayDate.Print();
        
        sTodayDate.Add(-1, Utilities::Date::YEAR);
        sTodayDate.Print();
        
        //  End of date test
    }
    else if (iChoice == 9)
    {
        std::cout << "Beginning of test of coverage" << std::endl;
        Utilities::Date::MyDate sToday, sInOneYear(3.0);
        Finance::MyBasis eBasis = Finance::BONDBASIS;
        Finance::Coverage sCoverage(eBasis, sToday, sInOneYear);
        
        std::cout << "Coverage : " << sCoverage.ComputeCoverage() << std::endl;
    }
    else if (iChoice == 10)
    {
        std::cout << "Beginning of test of annuity : " << std::endl;
        Utilities::Date::MyDate sToday, sIn10Years(10);
        Finance::MyFrequency eFrequency = Finance::MyFrequencyAnnual;
        Finance::MyBasis eBasis = Finance::BONDBASIS;
        
        std::cout << "YieldCurve Value : " << std::endl;
        double dYCValue = 0.02;
        std::cin >> dYCValue;
        Finance::YieldCurve sYC;
        sYC = dYCValue;
        
        Finance::Annuity sAnnuity(sToday, sIn10Years, eBasis, eFrequency, sYC);
        std::cout << "Annuity : " << sAnnuity.ComputeAnnuity() << std::endl;
    }
    else if (iChoice == 11)
    {}
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
        std::cout << "Enter the maturity of the caplet (in years): ";
        std::cin >> dMaturity;
        Utilities::require(dMaturity > 0 , "Maturity is negative");
		std::cout << "Enter Strike of caplet : ";
        std::cin >> dStrike;
        std::cout << "Enter Tenor (in years) : ";
        std::cin >> dTenor;
        Utilities::require(dTenor > 0, "Tenor is negative");
        Utilities::require(dTenor < dMaturity, "Tenor is higher than caplet maturity");
        /*std::cout << "Step by Step MC : ";
        std::cin >> iStepbyStepMC;*/

        //  Loop for test of multiples strikes
        /*for (std::size_t iStrike = 0 ; iStrike < 11 ; ++iStrike)
        {
            std::cout << "Strike : " << iStrike * 0.01 << ";";
            CapletPricingInterface(dMaturity, dTenor, iStrike * 0.01, iNPaths);
        }*/
        //  Simple strike test
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
        
        std::vector<double> dDFT1FwdNeutral;
        for (std::size_t iPath = 0; iPath < iNPaths0 ; ++iPath)
        {
            double dFactorT1FwdNeutral = sDataTForwardCube[iDate][iPath][0];
            dDFT1FwdNeutral.push_back(sLGM.BondPrice(dT1, dT2, dFactorT1FwdNeutral, Processes::FORWARD));
        }
         
        std::cout << "Forward bond price by simulation (T1 Forward Neutral) : " << sStats.Mean(dDFT1FwdNeutral) << std::endl;

        std::cout << "Bond Price value : " << exp(-sInitialYC.YC(dT2) * dT2) / exp(-sInitialYC.YC(dT1) * dT1) << std::endl;
        
        //  End of test of martingality of Bond price
    }
    else if (iChoice == 78)
    {
        //  Beginning of test of forward libor
        std::size_t iNPaths = 1000000;
        
        double dT2 = 1,  dT1 = 0.75;
        std::cout << "Start Date of Forward Libor : " << std::endl;
        std::cin >> dT1;
        std::cout << "Tenor of Libor : "<< std::endl;
        std::cin >> dT2;
        Utilities::require(dT2 > 0,"Maturity of discount factor is before start date");
        //  End date of libor
        dT2 += dT1;
        
        //  Initialization of LGM parameters 
        std::pair<std::vector<double>, std::vector<double> > dSigma;
        dSigma.first.push_back(0);
        dSigma.second.push_back(0.01); // volatility = 1%
        double dLambda = 0.05; // Mean reversion = 5%
        
        //  Initialization of classes
        Finance::TermStructure<double, double> sSigmaTS(dSigma.first, dSigma.second);
        Finance::YieldCurve sDiscountCurve, sSpreadCurve, sForwardCurve; // to change to SPLINE_CUBIC
        sSpreadCurve = 0.01;
        sDiscountCurve = 0.03;
        sForwardCurve = sSpreadCurve + sDiscountCurve;
        Processes::LinearGaussianMarkov sLGM(sDiscountCurve, sSpreadCurve, dLambda, sSigmaTS);

        Finance::SimulationData sSimulationData;
        std::vector<double> dSimulationTenors;
        dSimulationTenors.push_back(dT1);
        
        clock_t start = clock();
        //  Risk-neutral simulation
        sLGM.Simulate(iNPaths, dSimulationTenors, sSimulationData, /*Step by Step MC ?*/ true);
        std::cout << "Simulation Time : "<< (double)(clock()-start)/CLOCKS_PER_SEC <<" sec" <<std::endl;
        
        //  change of probability to T forward neutral
        Finance::SimulationData sSimulationDataTForward;
        sLGM.ChangeOfProbability(dT2, sSimulationData, sSimulationDataTForward);
        
        //  compute forward libor values
        std::vector<double> dForwardBondPrice;
        Finance::SimulationData::Cube sDataT2ForwardCube = sSimulationDataTForward.GetData().second;
        Finance::SimulationData::Cube sDataCube = sSimulationData.GetData().second;
        
        std::size_t iDate = 0;
        
        std::vector<double> dFactorRiskNeutral;
        std::size_t iNPaths0 = sDataCube[iDate].size();
        for (std::size_t iPath = 0 ; iPath < sDataCube[iDate].size() ; ++iPath)
        {
            dFactorRiskNeutral.push_back(sDataCube[iDate][iPath][0]);
        }
        //  Empirical distribution of factors
        Stats::Statistics sStats;
        
        std::vector<double> dLiborFwdT2Neutral;
        for (std::size_t iPath = 0; iPath < iNPaths0 ; ++iPath)
        {
            double dFactorT2FwdNeutral = sDataT2ForwardCube[iDate][iPath][0];
            dLiborFwdT2Neutral.push_back(sLGM.Libor(dT1, dT1, dT2, dFactorT2FwdNeutral, Processes::FORWARD));
        }
        
        std::cout << "Forward libor price by simulation (T2 Forward Neutral) : " << sStats.Mean(dLiborFwdT2Neutral) << std::endl;
        
        std::cout << "Bond Price value : " << 1.0 / (dT2 - dT1) * (exp(-sForwardCurve.YC(dT1) * dT1) / exp(-sForwardCurve.YC(dT2) * dT2) - 1.0) << std::endl;

        //  End of test of forward libor
    }
    else if (iChoice == 79)
    {
		std::vector<double> dFixingsA;
		dFixingsA.push_back(1.0);
		dFixingsA.push_back(8.0);
		std::vector<double> dValuesA;
		dValuesA.push_back(0.01);
		dValuesA.push_back(0.08);
		std::vector<double> dFixingsB;
		dFixingsB.push_back(0.0);
		dFixingsB.push_back(5.0);
		dFixingsB.push_back(6.0);
		dFixingsB.push_back(10.0);
		std::vector<double> dValuesB;
		dValuesB.push_back(0.0);
		dValuesB.push_back(0.05);
		dValuesB.push_back(0.06);
		dValuesB.push_back(0.1);
		
		Finance::TermStructure<double,double> TermStructureA(dFixingsA, dValuesA), TermStructureB(dFixingsB, dValuesB);
		
		TermStructureA.MergeTermStructure(TermStructureB);
		
		std::vector<double> dFixingsANew = TermStructureA.GetVariables();
		std::vector<double> dValuesANew = TermStructureA.GetValues();
		std::vector<double> dFixingsBNew = TermStructureB.GetVariables();
		std::vector<double> dValuesBNew = TermStructureB.GetValues();
		
		size_t iSize1 = dFixingsANew.size(), iSize2 = dValuesANew.size(), iSize3 = dFixingsBNew.size(), iSize4 = dValuesBNew.size();
		
		for (size_t i=0; i<iSize1; ++i) {
			std::cout << "FixingsANew : " << dFixingsANew[i] << std::endl;
		}
		for (size_t i=0; i<iSize2; ++i) {
			std::cout << "ValuesANew : " << dValuesANew[i] << std::endl;
		}
		for (size_t i=0; i<iSize3; ++i) {
			std::cout << "FixingsBNew : " << dFixingsBNew[i] << std::endl;
		}
		for (size_t i=0; i<iSize4; ++i) {
			std::cout << "ValuesBNew : " << dValuesBNew[i] << std::endl;
		}
		/*
         std::cout << "FixingsANew : " << dFixingsANew << std::endl;
         std::cout << "dValuesANew  : " << dValuesANew << std::endl;
         std::cout << "FixingsBNew : " << dFixingsBNew << std::endl;
         std::cout << "dValuesBNew : " << dValuesBNew << std::endl;
		 */
	}
    else if (iChoice == 80)
    {
        //  Beginning of test of quanto adjustment
        std::cout << "Do you want to try the Term-structure (0/1)? " << std::endl;
        std::size_t iTSChoice = 0;
        std::cin >> iTSChoice;
        
        double dFXVol = 0.10, dFXForCorrel = 1.0, dSigma = 0.01/*, dLambda = 0.05*/;
        //for (double dSigma = 0.001 ; dSigma < 0.041 ; dSigma += 0.001)
        //for (double dFXVol = 0.01 ; dFXVol < 0.41 ; dFXVol += 0.01)
        //for (double dFXForCorrel = -1 ; dFXForCorrel < 1.0 ; dFXForCorrel += 0.05)
        for (double dLambda = 0.005 ; dLambda < 0.05 ; dLambda += 0.001)
        {
            
            /*std::cout << "Start Date of Libor : " << std::endl;
             double dT1 = 1;
             std::cin >> dT1;
             Utilities::require(dT1 > 0, "Start Date of Libor is negative");
             
             std::cout << "Tenor of Libor : " << std::endl;
             double dT2 = 0.5;
             std::cin >> dT2;
             Utilities::require(dT2 > 0, "Tenor of libor is negative");*/
            double dT1 = 1, dT2 = 1;
            dT2 += dT1;
            
            std::pair<std::vector<double>, std::vector<double> >    dFXVolTS,
                                                                    dCorrelFXForTS, 
                                                                    dSigmaTS = std::make_pair(std::vector<double>(1, 0.0), std::vector<double>(1,dSigma));
            
            if (iTSChoice == 0)
            {
                dFXVolTS = std::make_pair(std::vector<double>(1,0.0), std::vector<double>(1,dFXVol));
                dCorrelFXForTS = std::make_pair(std::vector<double>(1,0.0), std::vector<double>(1,dFXForCorrel));
                
            }
            else if (iTSChoice == 1)
            {
                for (std::size_t i = 0 ; i < 10 ; ++i)
                {
                    dFXVolTS.first.push_back(i * 0.5);
                    dFXVolTS.second.push_back(dFXVol);
                    
                    dCorrelFXForTS.first.push_back(i * 0.5);
                    dCorrelFXForTS.second.push_back(dFXForCorrel);
                }
            }
            
            //  Initialization of LGM parameters 
            std::vector<std::pair<double, double> > dInitialYC;
            for (std::size_t i = 0 ; i < 4 * dT2 ; ++i)
            {
                dInitialYC.push_back(std::make_pair(0.25 * (i + 1), 0.03)); // YieldCurve = 3%
            }
            //double dLambda = 0.05; // Mean reversion = 5%
            
            //  Initialization of classes
            Finance::TermStructure<double, double> sFXVolTS(dFXVolTS.first, dFXVolTS.second);
            Finance::TermStructure<double, double> sSigmaTS(dSigmaTS.first, dSigmaTS.second);
            Finance::TermStructure<double, double> sFXForCorrelTS(dCorrelFXForTS.first, dCorrelFXForTS.second);
            Finance::YieldCurve sInitialYC("", "", dInitialYC, Utilities::Interp::LIN); // to change to SPLINE_CUBIC
            Processes::LinearGaussianMarkov sLGMDomestic(sInitialYC, dLambda, sSigmaTS), sLGMForeign(sInitialYC, dLambda, sSigmaTS);
            
            Processes::FXMultiCurve sFXMultiCurve(sFXVolTS, sLGMForeign, sLGMDomestic, sFXForCorrelTS);
            
            //std::cout << /*"Multiplicative Quanto Adj : " << */sFXMultiCurve.QuantoAdjustmentMultiplicative(dT1, dT2) << std::endl;
            printf("%.8lf\n", sFXMultiCurve.QuantoAdjustmentAdditive(dT1, dT2));
            //std::cout << /*"Additive Quanto Adj : " <<*/ sFXMultiCurve.QuantoAdjustmentAdditive(dT1, dT2) << std::endl;
        }
        //  End of test of quanto adjustment
    }
    else if (iChoice == 81)
    {
        Finance::TermStructure<double, double> sSigmaCollatTS, sSigmaOISTS;
        double dSigmaCollat = 0.01, dSigmaOIS;
        std::cout << "Volatility of Discounting IFR : " << std::endl;
        std::cin >> dSigmaOIS;
        
        std::cout << "Volatility of Forwarding IFR : " << std::endl;
        std::cin >> dSigmaCollat;
        sSigmaCollatTS = dSigmaCollat;
        sSigmaOISTS = dSigmaOIS;
        
        double dLambdaCollat, dLambdaOIS;
        std::cout << "Mean reversion Discounting IFR : " << std::endl;
        std::cin >> dLambdaOIS;
        
        std::cout << "Mean reversion Forwarding IFR : " << std::endl;
        std::cin >> dLambdaCollat;
        
        /*double dRhoCollatOIS;
        std::cout << "Correlation OIS Collat : " << std::endl;
        std::cin >> dRhoCollatOIS;*/
        
        double dT = 1, dt = 0;
        
        if (dLambdaOIS != 0 || dLambdaCollat != 0)
        {
            std::cout << "Maturity : " << std::endl;
            std::cin >> dT;
        }
        
        Processes::StochasticBasisSpread sStochasticBasisSpread;
        
        /*std::cout << "blablablablabla" << std::endl;
        std::cout << "Choose the output : " << std::endl;
        std::cout << "1- Correlation Spread-Discounting " << std::endl;
        std::cout << "2- Vol of spread" << std::endl;*/
        std::size_t iOutput = 2;
        //std::cin >> iOutput;
        
        if (iOutput == 1)
        {
            for (double dRhoCollatOIS = 0 ; dRhoCollatOIS <= 1.0 ; dRhoCollatOIS += 0.01)
            {
                std::cout << dRhoCollatOIS << ";" << sStochasticBasisSpread.CorrelationSpreadOIS(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dRhoCollatOIS, dt, dT) << std::endl;
            }
            /*for (double dRhoCollatOIS = 0 ; dRhoCollatOIS <= 1.0001 ; dRhoCollatOIS += 0.01)
            {
                std::cout << dRhoCollatOIS;
                for (dT = 1 ; dT <= 20.001 ; dT += 5)
                {
                    std::cout << ";" << sStochasticBasisSpread.CorrelationSpreadOIS(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dRhoCollatOIS, dt, dT) ;
                }
                std::cout << std::endl;
            }*/
        }
        else if (iOutput == 2)
        {
            /*for (double dRhoCollatOIS = 0 ; dRhoCollatOIS <= 1.0 ; dRhoCollatOIS += 0.01)
            {
                std::cout << dRhoCollatOIS << ";" << sStochasticBasisSpread.VolSpread(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dRhoCollatOIS, dt, dT) << std::endl;
             }*/
            for (double dRhoCollatOIS = 0 ; dRhoCollatOIS <= 1.00 ; dRhoCollatOIS += 0.01)
            {
                std::cout << dRhoCollatOIS;
                for (dT = 1 ; dT <= 20.001 ; dT += 5)
                {
                    std::cout << ";" << sStochasticBasisSpread.VolSpread(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dRhoCollatOIS, dt, dT) ;
                }
                std::cout << std::endl;
            }
        }
    }
    else if (iChoice == 82)
    {
        Finance::TermStructure<double, double> sSigmaCollatTS, sSigmaOISTS;
        double dSigmaCollat, dSigmaOIS;
        std::cout << "Volatility of OIS IFR : " << std::endl;
        std::cin >> dSigmaOIS;
        
        std::cout << "Volatility of Collat IFR : " << std::endl;
        std::cin >> dSigmaCollat;
        sSigmaCollatTS = dSigmaCollat;
        sSigmaOISTS = dSigmaOIS;
        
        double dLambdaCollat = 0.01, dLambdaOIS = 0.01;
        std::cout << "Mean reversion OIS IFR : " << std::endl;
        std::cin >> dLambdaOIS;
        
        std::cout << "Mean reversion Collat IFR : " << std::endl;
        std::cin >> dLambdaCollat;
        
        double dRhoCollatOIS;
        std::cout << "Correlation OIS Collat : " << std::endl;
        std::cin >> dRhoCollatOIS;
        
        double dT1 = 1, dT2 = 1.25, dt = 0;
        std::cout << "Reset Date of Libor (in Years) : " << std::endl;
        std::cin >> dT1;
        
        std::cout << "Tenor of Libor (in years) : " << std::endl;
        std::cin >> dT2;
        dT2 += dT1;
        
        Processes::StochasticBasisSpread sStochasticBasisSpread;
        
        for (std::size_t iIntervals = 100 ; iIntervals < 1100 ; iIntervals += 100)
        {
            std::cout << iIntervals << ";" << sStochasticBasisSpread.LiborQuantoAdjustmentMultiplicative(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dRhoCollatOIS, dt, dT1, dT2, iIntervals) << std::endl;
        }
    }
    else if (iChoice == 83)
    {
        std::cout << "Caplet Pricing with Stochastic Basis Spreads" << std::endl;
        
        double dT1 = 1, dT2 = 1.50, dt = 0;
        /*std::cout << "Reset Date of Libor (in Years) : " << std::endl;
        std::cin >> dT1;*/
        
        /*std::cout << "Tenor of Libor (in years) : " << std::endl;
        std::cin >> dT2;
        dT2 += dT1;*/
        
        //double dStrike = 0;
        /*std::cout << "Strike of Caplet : " << std::endl;
        std::cin >> dStrike;*/
        
        Finance::YieldCurve sDiscountCurve, sForwardingCurve;
        double dFlatDiscountCurveValue = 0.03, dFlatForwardingCurveValue = 0.03;
        /*std::cout << "Flat discount curve value : " << std::endl;
        std::cin >> dFlatDiscountCurveValue;*/
        sDiscountCurve = dFlatDiscountCurveValue;
        
        /*std::cout << "Flat Forwarding curve value : " << std::endl;
        std::cin >> dFlatForwardingCurveValue;*/
        sForwardingCurve = dFlatForwardingCurveValue;
        
        Finance::TermStructure<double, double> sSigmaCollatTS, sSigmaOISTS;
        double dSigmaCollat = 0.01, dSigmaOIS = 0.01;
        /*std::cout << "Volatility of OIS IFR : " << std::endl;
        std::cin >> dSigmaOIS;*/
        
        /*std::cout << "Volatility of Collat IFR : " << std::endl;
        std::cin >> dSigmaCollat;*/
        sSigmaCollatTS = dSigmaCollat;
        sSigmaOISTS = dSigmaOIS;
        
        double dLambdaCollat = 0.05, dLambdaOIS = 0.05;
        /*std::cout << "Mean reversion OIS IFR : " << std::endl;
        std::cin >> dLambdaOIS;
        
        std::cout << "Mean reversion Collat IFR : " << std::endl;
        std::cin >> dLambdaCollat;*/
        
        double dRhoCollatOIS = 0.89;
        /*std::cout << "Correlation OIS Collat : " << std::endl;
        std::cin >> dRhoCollatOIS;*/
        
        Finance::ForwardRate sForwardRate(sForwardingCurve);
        Finance::DF sDiscountDF(sDiscountCurve);
        Processes::StochasticBasisSpread sStochasticBasisSpread;
        std::size_t iNIntervals =  100; // computation of numerical integral for quanto adjustment
        
        double  dForwardRate = sForwardRate.FwdRate(dT1, dT2), dVolatilitySquare;
    
        if (dLambdaCollat < BETAOUTHRESHOLD)
        {
            dVolatilitySquare = 0;
        }
        else 
        {
            double dFirstTerm = MathFunctions::Beta_OU(2.0 * dLambdaCollat, dT1 - dt);
            double dSecondTerm = MathFunctions::Beta_OU(2.0 * dLambdaCollat, dT2 - dt) - MathFunctions::Beta_OU(2.0 * dLambdaCollat, dT2 - dT1);
            double dThirdTerm = MathFunctions::Beta_OU(dLambdaCollat, dT1 + dT2 - 2.0 * dt) - MathFunctions::Beta_OU(dLambdaCollat, dT2 - dT1);
            dVolatilitySquare = dSigmaCollat * dSigmaCollat / (dLambdaCollat * dLambdaCollat) * (- dFirstTerm - dSecondTerm + 2.0 * dThirdTerm);
        }
        double dVolatility = sqrt(dVolatilitySquare);
        std::cout << "Volatility : " << dVolatility / sqrt(dT1 - dt) << std::endl;
        std::cout << "Strike ; Difference PVStoch - PVNoStoch" << std::endl;
        for (double dStrike = 0.001 ; dStrike < 0.1 ; dStrike += 0.001)
        {
            double dQuantoAdj = sStochasticBasisSpread.LiborQuantoAdjustmentMultiplicative(sSigmaOISTS,
                                                                                           sSigmaCollatTS,
                                                                                           dLambdaOIS, 
                                                                                           dLambdaCollat,
                                                                                           dRhoCollatOIS, 
                                                                                           dt,
                                                                                           dT1,
                                                                                           dT2,
                                                                                           iNIntervals);
            
            double dPVStochBasisSpread = MathFunctions::BlackScholes(dForwardRate * dQuantoAdj,
                                                                     dStrike, 
                                                                     sqrt(dVolatilitySquare), 
                                                                     Finance::CALL)* sDiscountDF.DiscountFactor(dT2) * (dT2 - dT1),
            dPVNoStochBasisSpread = MathFunctions::BlackScholes(dForwardRate, 
                                                            dStrike, 
                                                            dVolatility, 
                                                            Finance::CALL) * sDiscountDF.DiscountFactor(dT2) * (dT2 - dT1);


            std::cout << dStrike << ";" << (dPVStochBasisSpread - dPVNoStochBasisSpread)/* / (dQuantoAdj - 1.0)*/ << /*";" << dQuantoAdj - 1.0 <<*/ std::endl;
        
            //std::cout << ";" << dPVNoStochBasisSpread << std::endl;
        
            //std::cout << "Relative difference : " << (dPVStochBasisSpread - dPVNoStochBasisSpread) / dPVNoStochBasisSpread << std::endl;
        }
    }
    else if (iChoice == 84)
    {
        std::cout << "Calculation Type: " << std::endl;
        std::cout << "1- Manual" << std::endl;
        std::cout << "2- Weight class" << std::endl;
        std::size_t iCalculationType = 1;
        std::cin >> iCalculationType;
        
        if (iCalculationType ==1)
        {
            // We will compute the weight in a particular case 
            //  5 year annuity
            
            //std::cout << "Yield curve value " << std::endl;
            //double dYCValue = 0.01;
            //std::cin >> dYCValue;
            //Finance::YieldCurve sYieldCurve;
            //sYieldCurve = dYCValue;
            std::vector<std::pair<double, double> > dYC;
            for (std::size_t i = 0 ; i < 101 ; ++i)
            {
                if (i < 20)
                {
                    dYC.push_back(std::make_pair(i * 0.1, 0.01));
                }
                else if (i < 40)
                {
                    dYC.push_back(std::make_pair(i * 0.1, 0.013));
                }
                else if (i < 60)
                {
                    dYC.push_back(std::make_pair(i * 0.1, 0.016));
                }
                else
                {
                    dYC.push_back(std::make_pair(i * 0.1, 0.02));
                }
            }
            Finance::YieldCurve sYC("", "", dYC, Utilities::Interp::LIN);
            
            Finance::DF sDF(sYC);
            
            double dT = 10;
            std::vector<double> dS;
            for (double i = 0 ; i < dT ; ++i)
            {
                dS.push_back(i + 5);
            }
            double dEnd = 1.0;
            for (double dDay = 0 ; dDay < dEnd ; dDay += 1.0 / 252)
            {
                double dLevel = 0;
                for (std::size_t i = 0 ; i < dS.size() ; ++i)
                {
                    dLevel += sDF.DiscountFactor(dS[i] - dDay);
                }
                std::cout << (int)(dDay * 252) << ";" ;
                std::cout << sDF.DiscountFactor(dS[0] - dDay) / dLevel << ";" ;
                std::cout << sDF.DiscountFactor(dS[1] - dDay) / dLevel << ";" ;
                std::cout << sDF.DiscountFactor(dS[2] - dDay) / dLevel << ";" ;
                std::cout << sDF.DiscountFactor(dS[3] - dDay) / dLevel << ";" ;
                std::cout << sDF.DiscountFactor(dS[4] - dDay) / dLevel << std::endl ;
                /*printf("%d", (int)(dDay * 252));
                 printf(";%.7lf", sDF.DiscountFactor(dS[0] - dDay) / dLevel);
                 
                 printf(";%.7lf", sDF.DiscountFactor(dS[1] - dDay) / dLevel);
                 
                 printf(";%.7lf", sDF.DiscountFactor(dS[2] - dDay) / dLevel);
                 
                 printf(";%.7lf", sDF.DiscountFactor(dS[3] - dDay) / dLevel);
                 
                 printf(";%.7lf\n", sDF.DiscountFactor(dS[4] - dDay) / dLevel);
                 */
            }
        }
        else if (iCalculationType == 2)
        {
        
        }
    }
    else if (iChoice == 86)
    {
        std::cout << "Swaption pricing with stochastic basis spread : " << std::endl;
        //  Pricing of a 5Y-10Y swaption
        std::cout << "Maturity of swaption (in years) : " << std::endl;
        double dSwaptionMaturity = 5;
        std::cin >> dSwaptionMaturity;
        
        std::cout << "Swap Length (in years) : " << std::endl;
        double dSwapLength = 10;
        std::cin >> dSwapLength;
        
        Finance::YieldCurve sForwardingCurve, sDiscountingCurve;
        sForwardingCurve = 0.03;
        sDiscountingCurve = 0.03;
        
        Finance::TermStructure<double, double> sSigmaf, sSigmad;
        double dVol = 0.01;
        sSigmad = dVol;
        sSigmaf = dVol;
        
        double dLambdad = 0.05, dLambdaf = 0.05;
        double dRhofd = 0.9;
        
        //  3M frequency on the floating leg
        std::vector<double> dS;
        for (std::size_t i = 0 ; i < dSwapLength * 4 ; ++i)
        {
            dS.push_back(0.25 * i + dSwaptionMaturity);
        }
        
        //  Annual frequency on the fixed leg
        std::vector<double> dT;
        for (std::size_t i = 0 ; i < dSwapLength; ++i) {
            dT.push_back(i + dSwaptionMaturity);
        }
        
        Processes::StochasticBasisSpread sStochasticBasisSpread;
        std::cout << "QuantoAdj : " << sStochasticBasisSpread.SwapQuantoAdjustmentMultiplicative(sSigmad, sSigmaf, dLambdad, dLambdaf, dRhofd, sDiscountingCurve, sForwardingCurve, 0, dS, dT) << std::endl;
    }
	else if (iChoice == 85) {
		double dFirstFixing = 1.0 ;
		double dFloatLegTenor = 0.25, dFixedLegTenor = 1.0 ;
		double dSwapMaturity = 10.0 ;
		
		double dLambdaOIS = 0.05, dLambdaCollat = 0.05 ;
		double dSigmaOIS = 0.05, dSigmaCollat = 0.05 ;
		double dRiskRateOIS = 0.05, dRiskRateCollat = 0.05 ;
		double dCorrelOISCollat = 0.8 ;
		
		// Calculation of dS, dT
		std::vector <double> dFloatLegFixing, dFixedLegFixing ;
		
		for (std::size_t iFixing = 0; iFixing <= dSwapMaturity / dFloatLegTenor; ++iFixing) {
			dFloatLegFixing.push_back(dFirstFixing + iFixing * dFloatLegTenor) ;
		}
		
		for (std::size_t iFixing = 0; iFixing <= dSwapMaturity / dFixedLegTenor; ++iFixing) {
			dFixedLegFixing.push_back(dFirstFixing + iFixing * dFixedLegTenor) ;
		}
		
		// Initialization of the yield curves
		Finance::YieldCurve sYieldCurveCollat ;
		sYieldCurveCollat = dRiskRateCollat ;
		
		Finance::YieldCurve sYieldCurveOIS ;
		sYieldCurveOIS = dRiskRateOIS ;
		
		// Initialization of the sigma term structures
		std::vector <double> dSigmaOISTS, dSigmaCollatTS, dSigmaFixing ;
		for (std::size_t iFixing = 0; iFixing <= dSwapMaturity; ++iFixing) {
			dSigmaFixing.push_back(iFixing) ;
			dSigmaOISTS.push_back(dSigmaOIS) ;
			dSigmaCollatTS.push_back(dSigmaCollat) ;
		}
		Finance::TermStructure <double, double> sSigmaOISTS(dSigmaFixing, dSigmaOISTS) ;
		Finance::TermStructure <double, double> sSigmaCollatTS(dSigmaFixing, dSigmaCollatTS) ;
		
		Processes::StochasticBasisSpread sStochasticBasisSpread ;
		std::cout << sStochasticBasisSpread.SwapQuantoAdjustmentMultiplicative(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dCorrelOISCollat, sYieldCurveOIS, sYieldCurveCollat, 0, dFloatLegFixing, dFixedLegFixing) << std::endl ;
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


