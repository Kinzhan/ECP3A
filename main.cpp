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
#include "Weights.h"
#include "SwapMonoCurve.h"
#include "CalibrationPms.h"

void CapletPricingInterface(const double dMaturity, const double dTenor, const double dStrike, std::size_t iNPaths, const double dLambda, const double dSigmaValue, const double dDiscountValue);

void CapletPricingInterface(const double dMaturity, const double dTenor, const double dStrike, std::size_t iNPaths, const double dLambda = 0.05, const double dSigmaValue = 0.01, const double dDiscountValue = 0.03)
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
    dSigma.second.push_back(dSigmaValue); // volatility = 1%
    //double dLambda = 0.05; // Mean reversion = 5%
    
    //  Initialization of classes
    Finance::TermStructure<double, double> sSigmaTS(dSigma.first, dSigma.second);
    //Finance::YieldCurve sInitialYC("", "", dInitialYC, Utilities::Interp::LIN); // to change to SPLINE_CUBIC
    Finance::YieldCurve sDiscountCurve, sSpreadCurve, sForwardCurve; 
    sDiscountCurve = dDiscountValue;
    sSpreadCurve = 0;
    sForwardCurve = sDiscountCurve + sSpreadCurve;
    Processes::LinearGaussianMarkov sLGM(sDiscountCurve, sSpreadCurve, dLambda, sSigmaTS);
    Finance::SimulationData sSimulationData;
    std::vector<double> dSimulationTenors;
    dSimulationTenors.push_back(dMaturity);
    
    clock_t start = clock();
    sLGM.Simulate(iNPaths, dSimulationTenors, sSimulationData, /*Step by Step MC ?*/ true);
    //std::cout << "Simulation Time : "<< (double)(clock()-start)/CLOCKS_PER_SEC <<" sec" <<std::endl;
    start = clock();
    
    // PRINT FACTORS AFTER CHANGE OF PROBA
    Finance::SimulationData sSimulationDataTForward;
    sLGM.ChangeOfProbability(dMaturity, sSimulationData, sSimulationDataTForward); 
    
    //std::cout << "Change of probability time : " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << std::endl;
    
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
    //Utilities::PrintInFile sPrint("/Users/alexhum49/Desktop/TextCaplet.txt", false, 6);
	/*Utilities::PrintInFile sPrint("/Users/kinzhan/Desktop/TextCaplet.txt", false, 6);
    sPrint.PrintDataInFile(dMeanPrice);
    std::cout << "Print in file : done ! "<< std::endl;*/
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
    
    //std::cout<<"Total Time elapsed : " << (double)(clock()-start)/CLOCKS_PER_SEC <<" sec"<< std::endl;
}

//  Volatility of zero coupon bond P(t,T)
double Gamma(double dLambda, double dSigma, double dt, double dT);
double Gamma(double dLambda, double dSigma, double dt, double dT)
{
    if (dLambda > 1.0e-07)
    {
        return dSigma / dLambda * (1. - exp(-dLambda * (dT - dt)));
    }
    else 
    {
        return dSigma * (dT - dt);
    }
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
    std::cout << "11- Calibration of 1F Model" << std::endl;
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
    std::cout << "86- Swaption Pricing with Stochastic Basis Spread" << std::endl;
    std::cout << "87- Volatility of weights" << std::endl;
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
        dVariablesAndValues.push_back(std::make_pair(1.0, 1.0));
        dVariablesAndValues.push_back(std::make_pair(2.0, 6.0));
		dVariablesAndValues.push_back(std::make_pair(9.0, 8.0));
		dVariablesAndValues.push_back(std::make_pair(5.0, 7.0));
        dVariablesAndValues.push_back(std::make_pair(3.5, 3.0));
        //dVariablesAndValues.push_back(std::make_pair(5.0, 7.0));
        //dVariablesAndValues.push_back(std::make_pair(9.0, 8.0));
        //Finance::YieldCurve sYieldCurve("EUR", "EUROIS", dVariablesAndValues, Utilities::Interp::SPLINE_CUBIC);
        Finance::YieldCurve sYieldCurve("EUR", "EUROIS", dVariablesAndValues, Utilities::Interp::LIN);
		
        double dt = 0;
        std::cout<<"Enter value : ";
        std::cin>>dt;
        std::cout << "Yield Curve value is : " << sYieldCurve.YC(dt) << std::endl;

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
    {
        //  Beginning of calibration of 1F model parameters
        try 
        {
            Calibration::CalibrationPms sCalib;
            std::string cFileName = "/Users/alexhum49/Desktop/Libor6M.txt";
            std::vector<double> dData = sCalib.LoadDataFromFile(cFileName);
            Calibration::NewtonPms sNewtonPms(0.0001, 100);
            double dDeltaT = 1/12.;
            sCalib.NewtonRaphsonAlgorithm(dData, dDeltaT, sNewtonPms);
        } 
        catch (const std::string & cError) 
        {
            std::cout << cError << std::endl;
        }
        catch (const char * cError)
        {
            std::cout << cError << std::endl;
        }
        //  End of calibration of 1F model parameters
    }
    else if (iChoice == 75)
    {
        //  Test for Hull-White model
        
        std::cout << "Caplet Pricing by simulation" << std::endl;
        std::size_t iNPaths = 1000000;
        double dMaturity = 4.0, dTenor = 2, dStrike = 0.02;
        //std::size_t iStepbyStepMC = false;
        
        //  Input some variables
        /*std::cout << "Enter the number of paths : ";
        std::cin >> iNPaths;
        std::cout << "Enter the maturity of the caplet (in years): ";
        std::cin >> dMaturity;
        Utilities::require(dMaturity > 0 , "Maturity is negative");
		std::cout << "Enter Strike of caplet : ";
        std::cin >> dStrike;
        std::cout << "Enter Tenor (in years) : ";
        std::cin >> dTenor;*/
        Utilities::require(dTenor > 0, "Tenor is negative");
        //Utilities::require(dTenor < dMaturity, "Tenor is higher than caplet maturity");
        /*std::cout << "Step by Step MC : ";
        std::cin >> iStepbyStepMC;*/

        //  Loop for test of multiples strikes
        /*for (std::size_t iStrike = 0 ; iStrike < 11 ; ++iStrike)
        {
            std::cout << "Strike : " << iStrike * 0.01 << ";";
            CapletPricingInterface(dMaturity, dTenor, iStrike * 0.01, iNPaths);
        }*/
        /*//  Simple strike test
		for (double dRiskValue = 0.002; dRiskValue <= 0.05; dRiskValue+=0.002) {
			CapletPricingInterface(dMaturity, dTenor, dStrike, iNPaths, 0.25, 0.01, dRiskValue);
		}*/
        CapletPricingInterface(dMaturity, dTenor, dStrike, iNPaths, 0.25, 0.01, 0.03);
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
        double dT2 = 4,  dT1 = 2;
        /*std::cout << "Start Date of Forward discout factor : " << std::endl;
        std::cin >> dT1;
        std::cout << "Maturity of discount factor : "<< std::endl;
        std::cin >> dT2;*/
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
        double dLambda = 1; // Mean reversion = 5%
        
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
        //sLGM.ChangeOfProbability(0, sSimulationData, sSimulationDataTForward);
        
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
		
		// For the report charts, comparison of the ZC forward under risk neutral \ forward probabilities
		std::vector<double> dBuckets ;
		{
			dBuckets.push_back(0.835342) ;
			dBuckets.push_back(0.837599) ;
			dBuckets.push_back(0.839856) ;
			dBuckets.push_back(0.842113) ;
			dBuckets.push_back(0.84437) ;
			dBuckets.push_back(0.846627) ;
			dBuckets.push_back(0.848884) ;
			dBuckets.push_back(0.851141) ;
			dBuckets.push_back(0.853398) ;
			dBuckets.push_back(0.855655) ;
			dBuckets.push_back(0.857912) ;
			dBuckets.push_back(0.860169) ;
			dBuckets.push_back(0.862426) ;
			dBuckets.push_back(0.864684) ;
			dBuckets.push_back(0.866941) ;
			dBuckets.push_back(0.869198) ;
			dBuckets.push_back(0.871455) ;
			dBuckets.push_back(0.873712) ;
			dBuckets.push_back(0.875969) ;
			dBuckets.push_back(0.878226) ;
			dBuckets.push_back(0.880483) ;
			dBuckets.push_back(0.88274) ;
			dBuckets.push_back(0.884997) ;
			dBuckets.push_back(0.887254) ;
			dBuckets.push_back(0.889511) ;
			dBuckets.push_back(0.891768) ;
			dBuckets.push_back(0.894025) ;
			dBuckets.push_back(0.896283) ;
			dBuckets.push_back(0.89854) ;
			dBuckets.push_back(0.900797) ;
			dBuckets.push_back(0.903054) ;
			dBuckets.push_back(0.905311) ;
			dBuckets.push_back(0.907568) ;
			dBuckets.push_back(0.909825) ;
			dBuckets.push_back(0.912082) ;
			dBuckets.push_back(0.914339) ;
			dBuckets.push_back(0.916596) ;
			dBuckets.push_back(0.918853) ;
			dBuckets.push_back(0.92111) ;
			dBuckets.push_back(0.923367) ;
			dBuckets.push_back(0.925625) ;
			dBuckets.push_back(0.927882) ;
			dBuckets.push_back(0.930139) ;
			dBuckets.push_back(0.932396) ;
			dBuckets.push_back(0.934653) ;
			dBuckets.push_back(0.93691) ;
			dBuckets.push_back(0.939167) ;
			dBuckets.push_back(0.941424) ;
			dBuckets.push_back(0.943681) ;
			dBuckets.push_back(0.945938) ;
			dBuckets.push_back(0.948195) ;
			dBuckets.push_back(0.950452) ;
			dBuckets.push_back(0.952709) ;
			dBuckets.push_back(0.954966) ;
			dBuckets.push_back(0.957224) ;
			dBuckets.push_back(0.959481) ;
			dBuckets.push_back(0.961738) ;
			dBuckets.push_back(0.963995) ;
			dBuckets.push_back(0.966252) ;
			dBuckets.push_back(0.968509) ;
			dBuckets.push_back(0.970766) ;
			dBuckets.push_back(0.973023) ;
			dBuckets.push_back(0.97528) ;
			dBuckets.push_back(0.977537) ;
			dBuckets.push_back(0.979794) ;
			dBuckets.push_back(0.982051) ;
			dBuckets.push_back(0.984308) ;
			dBuckets.push_back(0.986565) ;
			dBuckets.push_back(0.988823) ;
			dBuckets.push_back(0.99108) ;
			dBuckets.push_back(0.993337) ;
			dBuckets.push_back(0.995594) ;
			dBuckets.push_back(0.997851) ;
			dBuckets.push_back(1.000108) ;
			dBuckets.push_back(1.002365) ;
			dBuckets.push_back(1.004622) ;
			dBuckets.push_back(1.006879) ;
			dBuckets.push_back(1.009136) ;
			dBuckets.push_back(1.011393) ;
			dBuckets.push_back(1.01365) ;
			dBuckets.push_back(1.015907) ;
			dBuckets.push_back(1.018164) ;
			dBuckets.push_back(1.020422) ;
			dBuckets.push_back(1.022679) ;
			dBuckets.push_back(1.024936) ;
			dBuckets.push_back(1.027193) ;
			dBuckets.push_back(1.02945) ;
			dBuckets.push_back(1.031707) ;
			dBuckets.push_back(1.033964) ;
			dBuckets.push_back(1.036221) ;
			dBuckets.push_back(1.038478) ;
			dBuckets.push_back(1.040735) ;
			dBuckets.push_back(1.042992) ;
			dBuckets.push_back(1.045249) ;
			dBuckets.push_back(1.047506) ;
			dBuckets.push_back(1.049763) ;
			dBuckets.push_back(1.052021) ;
			dBuckets.push_back(1.054278) ;
			dBuckets.push_back(1.056535) ;
			dBuckets.push_back(1.058792) ;
		}
		
		std::vector<std::pair<double, std::size_t> > sDistribution = sStats.EmpiricalDistribution(dDFT1FwdNeutral,100);
		//std::vector<std::pair<double, std::size_t> > sDistribution = sStats.EmpiricalDistribution(dDFT1FwdNeutral,dBuckets);
		
		Utilities::PrintInFile sPrintDistribution("/Users/kinzhan/Desktop/ZCEmpiricalDistribution.txt", false, 6);
		sPrintDistribution.PrintDataInFile(sDistribution);
        
        //  End of test of martingality of Bond price
    }
    else if (iChoice == 78)
    {
        //  Beginning of test of forward libor
        std::size_t iNPaths = 1000000;
        
        double dT2 = 1,  dT1 = 2;
        //std::cout << "Start Date of Forward Libor : " << std::endl;
        //std::cin >> dT1;
        //std::cout << "Tenor of Libor : "<< std::endl;
        //std::cin >> dT2;
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
        sSpreadCurve = 0.00;
        sDiscountCurve = 0.03;
        //sForwardCurve = sSpreadCurve + sDiscountCurve;
		sForwardCurve = sDiscountCurve;
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
        //sLGM.ChangeOfProbability(0, sSimulationData, sSimulationDataTForward);
		
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
		
		// For the report charts, comparison of the ZC forward under risk neutral \ forward probabilities
		std::vector<double> dBuckets ;
		{
			dBuckets.push_back(-0.037035);
			dBuckets.push_back(-0.03558);
			dBuckets.push_back(-0.034125);
			dBuckets.push_back(-0.03267);
			dBuckets.push_back(-0.031215);
			dBuckets.push_back(-0.029761);
			dBuckets.push_back(-0.028306);
			dBuckets.push_back(-0.026851);
			dBuckets.push_back(-0.025396);
			dBuckets.push_back(-0.023941);
			dBuckets.push_back(-0.022487);
			dBuckets.push_back(-0.021032);
			dBuckets.push_back(-0.019577);
			dBuckets.push_back(-0.018122);
			dBuckets.push_back(-0.016667);
			dBuckets.push_back(-0.015213);
			dBuckets.push_back(-0.013758);
			dBuckets.push_back(-0.012303);
			dBuckets.push_back(-0.010848);
			dBuckets.push_back(-0.009393);
			dBuckets.push_back(-0.007939);
			dBuckets.push_back(-0.006484);
			dBuckets.push_back(-0.005029);
			dBuckets.push_back(-0.003574);
			dBuckets.push_back(-0.002119);
			dBuckets.push_back(-0.000665);
			dBuckets.push_back(0.00079);
			dBuckets.push_back(0.002245);
			dBuckets.push_back(0.0037);
			dBuckets.push_back(0.005155);
			dBuckets.push_back(0.006609);
			dBuckets.push_back(0.008064);
			dBuckets.push_back(0.009519);
			dBuckets.push_back(0.010974);
			dBuckets.push_back(0.012429);
			dBuckets.push_back(0.013883);
			dBuckets.push_back(0.015338);
			dBuckets.push_back(0.016793);
			dBuckets.push_back(0.018248);
			dBuckets.push_back(0.019703);
			dBuckets.push_back(0.021157);
			dBuckets.push_back(0.022612);
			dBuckets.push_back(0.024067);
			dBuckets.push_back(0.025522);
			dBuckets.push_back(0.026977);
			dBuckets.push_back(0.028431);
			dBuckets.push_back(0.029886);
			dBuckets.push_back(0.031341);
			dBuckets.push_back(0.032796);
			dBuckets.push_back(0.034251);
			dBuckets.push_back(0.035705);
			dBuckets.push_back(0.03716);
			dBuckets.push_back(0.038615);
			dBuckets.push_back(0.04007);
			dBuckets.push_back(0.041525);
			dBuckets.push_back(0.042979);
			dBuckets.push_back(0.044434);
			dBuckets.push_back(0.045889);
			dBuckets.push_back(0.047344);
			dBuckets.push_back(0.048799);
			dBuckets.push_back(0.050253);
			dBuckets.push_back(0.051708);
			dBuckets.push_back(0.053163);
			dBuckets.push_back(0.054618);
			dBuckets.push_back(0.056073);
			dBuckets.push_back(0.057527);
			dBuckets.push_back(0.058982);
			dBuckets.push_back(0.060437);
			dBuckets.push_back(0.061892);
			dBuckets.push_back(0.063347);
			dBuckets.push_back(0.064801);
			dBuckets.push_back(0.066256);
			dBuckets.push_back(0.067711);
			dBuckets.push_back(0.069166);
			dBuckets.push_back(0.070621);
			dBuckets.push_back(0.072075);
			dBuckets.push_back(0.07353);
			dBuckets.push_back(0.074985);
			dBuckets.push_back(0.07644);
			dBuckets.push_back(0.077895);
			dBuckets.push_back(0.079349);
			dBuckets.push_back(0.080804);
			dBuckets.push_back(0.082259);
			dBuckets.push_back(0.083714);
			dBuckets.push_back(0.085169);
			dBuckets.push_back(0.086623);
			dBuckets.push_back(0.088078);
			dBuckets.push_back(0.089533);
			dBuckets.push_back(0.090988);
			dBuckets.push_back(0.092443);
			dBuckets.push_back(0.093897);
			dBuckets.push_back(0.095352);
			dBuckets.push_back(0.096807);
			dBuckets.push_back(0.098262);
			dBuckets.push_back(0.099717);
			dBuckets.push_back(0.101171);
			dBuckets.push_back(0.102626);
			dBuckets.push_back(0.104081);
			dBuckets.push_back(0.105536);
			dBuckets.push_back(0.106991);
		}
		
		//std::vector<std::pair<double, std::size_t> > sDistribution = sStats.EmpiricalDistribution(dLiborFwdT2Neutral,100);
		std::vector<std::pair<double, std::size_t> > sDistribution = sStats.EmpiricalDistribution(dLiborFwdT2Neutral,dBuckets);
		
		Utilities::PrintInFile sPrintDistribution("/Users/kinzhan/Desktop/LiborEmpiricalDistribution.txt", false, 6);
		sPrintDistribution.PrintDataInFile(sDistribution);
		
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
        double dSigmaCollat = 0.01, dSigmaOIS = 0.01;
        /*std::cout << "Volatility of Discounting IFR : " << std::endl;
        std::cin >> dSigmaOIS;
        
        std::cout << "Volatility of Forwarding IFR : " << std::endl;
        std::cin >> dSigmaCollat;*/
        sSigmaCollatTS = dSigmaCollat;
        sSigmaOISTS = dSigmaOIS;
        
        double dLambdaCollat = 0.01, dLambdaOIS = 0.01;
        /*std::cout << "Mean reversion Discounting IFR : " << std::endl;
        std::cin >> dLambdaOIS;
        
        std::cout << "Mean reversion Forwarding IFR : " << std::endl;
        std::cin >> dLambdaCollat;*/
        
        double dRhoCollatOIS = 0.5;
        /*std::cout << "Correlation OIS Collat : " << std::endl;
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
        double dSigmaCollat = 0.1, dSigmaOIS = 0.01;
        /*std::cout << "Volatility of OIS IFR : " << std::endl;
        std::cin >> dSigmaOIS;
        
        std::cout << "Volatility of Collat IFR : " << std::endl;
        std::cin >> dSigmaCollat;*/
        sSigmaCollatTS = dSigmaCollat;
        sSigmaOISTS = dSigmaOIS;
        
        double dLambdaCollat = 0.1, dLambdaOIS = 0.05;
        /*std::cout << "Mean reversion OIS IFR : " << std::endl;
        std::cin >> dLambdaOIS;
        
        std::cout << "Mean reversion Collat IFR : " << std::endl;
        std::cin >> dLambdaCollat;*/
        
        double dRhoCollatOIS = 0.95;
        /*std::cout << "Correlation OIS Collat : " << std::endl;
        std::cin >> dRhoCollatOIS;*/
        
        double dT1 = 1, dT2 = 2., dt = 0.;

        /*std::cout << "Reset Date of Libor (in Years) : " << std::endl;
        std::cin >> dT1;
        
        std::cout << "Tenor of Libor (in years) : " << std::endl;
        std::cin >> dT2;*/
        //dT2 += dT1;
        
        Processes::StochasticBasisSpread sStochasticBasisSpread;
        
        /*for (std::size_t iIntervals = 100 ; iIntervals < 1100 ; iIntervals += 100)
        {
            std::cout << iIntervals << ";" << sStochasticBasisSpread.LiborQuantoAdjustmentMultiplicative(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dRhoCollatOIS, dt, dT1, dT2, iIntervals) << std::endl;
        }*/
		
		std::size_t iIntervals = 300 ;
		
		for (double dSigmaCollatValue = 0.01 ; dSigmaCollatValue <= 0.101 ; dSigmaCollatValue += 0.001)
        {
			sSigmaCollatTS = dSigmaCollatValue;
            //for (double dLambdaCollatValue = 0.01 ; dLambdaCollatValue <= 0.101 ; dLambdaCollatValue += 0.01)
			{
            double dLambdaCollatValue = 0.05;
				dLambdaCollat = dLambdaCollatValue;
            double dRhoValue = 0.8;
				//for (double dRhoValue = 0.05 ; dRhoValue <= 1 ; dRhoValue += 0.10)
				{
					dRhoCollatOIS = dRhoValue;
					std::cout << dSigmaCollatValue << ";" << dLambdaCollatValue << ";" << dRhoValue << ";" << sStochasticBasisSpread.LiborQuantoAdjustmentMultiplicative(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dRhoCollatOIS, dt, dT1, dT2, iIntervals) - 1 << std::endl;
				}
			}
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
        double dT = 10;
        std::vector<double> dS;
        for (double i = 0 ; i < dT ; ++i)
        {
            dS.push_back(i + 5);
        }
        
        if (iCalculationType ==1)
        {
            // We will compute the weight in a particular case 
            //  5 year annuity
            
            Finance::DF sDF(sYC);
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
            Finance::Weights sWeights(sYC, dS);
            std::vector<double> dWeights = sWeights.GetWeights();
            for (std::size_t i = 0; i < dWeights.size() ; ++i)
            {
                std::cout << i << ";" << dWeights[i] << std::endl;
            }
        }
    }
    else if (iChoice == 85) 
    {
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
    else if (iChoice == 86)
    {
        std::cout << "Swaption pricing with stochastic basis spread : " << std::endl;
        //  Pricing of a 5Y-10Y swaption
        //std::cout << "Maturity of swaption (in years) : " << std::endl;
        //double dSwaptionMaturity = 5;
        //std::cin >> dSwaptionMaturity;
        
        //std::cout << "Swap Length (in years) : " << std::endl;
        //double dSwapLength = 10;
        //std::cin >> dSwapLength;
        
        std::cout << "Choose the chart " << std::endl;
        std::cout << "1- Additive QA (Length And Maturity)" << std::endl;
        std::cout << "2- Swap values" << std::endl;
        std::cout << "3- Swaption Pricing" << std::endl;
        std::cout << "4- Additive QA (Model Parameters)" << std::endl;
        std::size_t iTest = 1;
        std::cin >> iTest;
        
        Finance::YieldCurve sForwardingCurve, sDiscountingCurve;
        sDiscountingCurve = 0.03;
        sForwardingCurve = 0.03;
        
        std::cout << "Flat Forwarding curve ? (0/1)" << std::endl;
        bool bFlatForwardingCurve = true;
        std::cin >> bFlatForwardingCurve;
        
        if (!bFlatForwardingCurve)
        {
            std::vector<std::pair<double, double> > dVectOfPair;
            for (std::size_t i = 0 ; i < 31 ; ++i)
            {
                dVectOfPair.push_back(std::make_pair(i, -0.005 + 0.0001 * i / 30));
            }
            Finance::YieldCurve sYC("", "", dVectOfPair, Utilities::Interp::LIN);
            sForwardingCurve = sForwardingCurve + sYC;
        }
        
        Finance::TermStructure<double, double> sSigmaf, sSigmad;
        double dVol = 0.01;
        sSigmad = dVol;
        sSigmaf = dVol;
        
        double dLambdad = 0.05, dLambdaf = 0.05;
        double dRhofd = 0.9;
        
        std::cout << "SwapLength ; SwapBegin ; ";
        
        //  3M frequency on the floating leg
        if (iTest == 1)
        {
            std::cout << "Additive QA " << std::endl;
            double dEpsilon = 0.001;
            for (double dSwapLength = 1 ; dSwapLength < 21 ; ++dSwapLength)
            {
                for (double dSwaptionMaturity = 1 ; dSwaptionMaturity < 10 ; ++dSwaptionMaturity)
                {
                    std::vector<double> dS;
                    for (std::size_t i = 0 ; i <= dSwapLength * 4 + dEpsilon; ++i)
                    {
                        dS.push_back(0.25 * i + dSwaptionMaturity);
                    }
                    
                    //  Annual frequency on the fixed leg
                    std::vector<double> dT;
                    for (std::size_t i = 0 ; i < dSwapLength + dEpsilon ; ++i)
                    {
                        dT.push_back(i + dSwaptionMaturity);
                    }
                    
                    Processes::StochasticBasisSpread sStochasticBasisSpread;
                    std::cout << /*"QuantoAdj : "*/ dSwapLength << ";" << dSwaptionMaturity << ";" << sStochasticBasisSpread.SwapQuantoAdjustmentMultiplicative(sSigmad, sSigmaf, dLambdad, dLambdaf, dRhofd, sDiscountingCurve, sForwardingCurve, 0, dS, dT) - 1. << std::endl;
                }
            }
        }
        else if (iTest == 2)
        {
            std::cout << "Swap Values " << std::endl;
            for (double dSwapLength = 1 ; dSwapLength < 21 ; ++dSwapLength)
            {
                for (double dSwaptionMaturity = 1 ; dSwaptionMaturity < 11 ; ++dSwaptionMaturity)
                {
                    Utilities::Date::MyDate sDateBegin(dSwaptionMaturity), sDateEnd(dSwaptionMaturity + dSwapLength);
                    Finance::Annuity sAnnuity(sDateBegin, sDateEnd, Finance::BONDBASIS, Finance::MyFrequencyAnnual, sForwardingCurve);
                    Finance::SwapMonoCurve sSwapMonoCurve(sDateBegin, sDateEnd, Finance::MyFrequencyAnnual, Finance::BONDBASIS, sForwardingCurve);
                    
                    
                    std::cout << dSwapLength << ";" << dSwaptionMaturity << ";" << sSwapMonoCurve.ComputeSwap() << std::endl;
                }
            }
        }
        else if (iTest == 3)
        {
            std::cout << "Swaption Premia" << std::endl;
            std::cout << "Not yet implemented" << std::endl;
        }
        else if (iTest == 4)
        {
            // We take a 10Y swap beginning in 5Y
            //  Frequency Annual
            //
            std::cout << std::endl;
            double dBegin = 5., dLength = 10., dEpsilon = 0.00001;
            std::vector<double> dT, dS;
            for (std::size_t i = 0 ; i < 4 * dLength + dEpsilon ; ++i)
            {
                dT.push_back(dBegin + i * 0.25);
            }
            for (std::size_t i = 0 ; i < dLength + dEpsilon ; ++i)
            {
                dS.push_back(dBegin + i);
            }
            
            double dLambdaCollat, dRhoCollatOIS, dLambdaOIS = 0.05;
            Finance::TermStructure<double, double> sSigmaCollatTS, sSigmaOISTS;
            double dSigmaOIS = 0.02;
            sSigmaOISTS = dSigmaOIS;
            Processes::StochasticBasisSpread sStochasticBasisSpread;
            double dt = 0;
            
            for (double dSigmaCollatValue = 0.01 ; dSigmaCollatValue <= 0.03 ; dSigmaCollatValue += 0.002)
            {
                //double dSigmaCollatValue = 0.009;
                sSigmaCollatTS = dSigmaCollatValue;
                for (double dLambdaCollatValue = 0.01 ; dLambdaCollatValue <= 0.101 ; dLambdaCollatValue += 0.01)
                {
                    //double dLambdaCollatValue = 0.05;
                    dLambdaCollat = dLambdaCollatValue;
                    //double dRhoValue = 0.8;
                    for (double dRhoValue = 0.01 ; dRhoValue <= 1.01 ; dRhoValue += 0.10)
                    {
                        dRhoCollatOIS = dRhoValue;
                        std::cout << dSigmaCollatValue << ";" << dLambdaCollatValue << ";" << dRhoValue << ";" << sStochasticBasisSpread.SwapQuantoAdjustmentMultiplicative(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dRhoCollatOIS, sDiscountingCurve, sForwardingCurve, dt, dS, dT) - 1 << std::endl;
                    }
                }
            }
        }
    }
    else if (iChoice == 87)
    {
        //  Swap Length = 10Y, Frequency 3M (Floating Leg), Annual on the fixed leg --> 40 dates in the schedule
        std::size_t iNDatesFixedLeg = 10, iNDatesFloatingLeg = 40;
        std::vector<double> dDatesFixedLeg, dDatesFloatingLeg;
        for (std::size_t i = 0 ; i <= iNDatesFixedLeg ; ++i)
        {
            dDatesFixedLeg.push_back(i);
        }
        for (std::size_t i = 1 ; i <= iNDatesFloatingLeg ; ++i)
        {
            dDatesFloatingLeg.push_back(i * 0.25);
        }
        
        double dYCValue = 0.03;
        Finance::YieldCurve sYieldCurve;
        sYieldCurve = dYCValue;
        
        Finance::Weights sWeights(sYieldCurve, dDatesFixedLeg);
        
        double dLambda = 0.05;
        double dSigma = 0.01;
        
        double dAnnuityVol = 0;
        for (std::size_t i = 0 ; i < iNDatesFixedLeg ; ++i)
        {
            double dGamma_i = Gamma(dLambda, dSigma, 0, dDatesFixedLeg[i+1]);
            dAnnuityVol += sWeights.GetWeight(i) * dGamma_i ;
        }
        
        std::cout << "Date ; FwdLiborVolatility ; WeightVolatility" << std::endl;
        for (std::size_t i = 0 ; i < iNDatesFloatingLeg ; ++i)
        {
            double dFwdLiborVolatility = Gamma(dLambda, dSigma, 0, dDatesFloatingLeg[i+1]) - Gamma(dLambda, dSigma, 0, dDatesFloatingLeg[i]);
            double dWeightVolatility = - Gamma(dLambda, dSigma, 0, dDatesFloatingLeg[i+1]) + dAnnuityVol;
            std::cout << i << ";" << dFwdLiborVolatility << ";" << dWeightVolatility << std::endl;
        }
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


