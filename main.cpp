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

int main()
{
    std::vector<double> dRealisations;
    
    std::cout << "Hello" << std::endl;
    std::size_t iNRealisations = 10, iAntitheticVariables = 1, iChoice;
    
    std::cout<< "Which variables do you want to simulate (1: Uniform, 2: Gaussian, 3: Interpolation, 4: YieldCurve, 5: SimulationData, 75: HullWhite, 76: Test) ? ";
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
    else if (iChoice ==5)
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
    else if (iChoice == 75)
    {
        //  Test for Hull-White model
        
        std::cout << "Caplet Pricing by simulation" << std::endl;
        std::size_t iNPaths = 10;
        double dMaturity = 1.0, dStrike = 0.01, dTenor = 0.25;
        //std::size_t iStepbyStepMC = false;
        
        //  Input some variables
        std::cout << "Enter the number of paths : ";
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
        Utilities::require(dTenor < dMaturity, "Tenor is higher than caplet maturity");
        /*std::cout << "Step by Step MC : ";
        std::cin >> iStepbyStepMC;*/

        //  Initialization of LGM parameters 
        std::pair<std::vector<double>, std::vector<double> > dSigma;
        std::vector<std::pair<double, double> > dInitialYC;
        for (std::size_t i = 0 ; i < 4 * dMaturity ; ++i)
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
        for (std::size_t i = 0 ; i < dMaturity * 10 ; ++i)
        {
            dSimulationTenors.push_back(0.1 * (i + 1));
        }
        
        clock_t start = clock();
        sLGM.Simulate(iNPaths, dSimulationTenors, sSimulationData, /*Step by Step MC ?*/ /*iStepbyStepMC*/ true);
        std::cout << "Simulation Time : "<< (double)(clock()-start)/CLOCKS_PER_SEC <<" sec" <<std::endl;
        
        //sSimulationData.PrintInFile("/Users/alexhum49/Desktop/MyfileRiskNeutral.txt", false);
		sSimulationData.PrintInFile("/Users/Kinz/Desktop/MyFilecopy.txt", false);
		
		// PRINT FACTORS AFTER CHANGE OF PROBA
		Finance::SimulationData sSimulationDataTForward;
        sLGM.ChangeOfProbability(dMaturity, sSimulationData, sSimulationDataTForward);    
        //sSimulationDataTForward.PrintInFile("/Users/alexhum49/Desktop/MyfileTForward.txt", false);
		sSimulationDataTForward.PrintInFile("/Users/Kinz/Desktop/MyfileRiskNeutral.txt", false);
		

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
        double dPrice = 0.0, dStdDevPrice = 0.0;
        std::vector<double> dPayoff = sProductLGM.Caplet(dMaturity - dTenor, dMaturity, dStrike, sSimulationData);
        /*for (std::size_t iPath = 0 ; iPath < iNPaths ; ++iPath)
        {
            //  Payment Date discount factor
            //std::vector<double> dDFs = sProductLGM.RiskNeutralDiscountFactor(iPath, sSimulationData);
            //double dDF = dDFs.back();
            double dDF = 1.0;
			dPrice += dDF * dPayoff[iPath];
            dStdDevPrice += dDF * dPayoff[iPath] * dDF * dPayoff[iPath];
        }
        */
        dPrice /= iNPaths;
        dStdDevPrice /= iNPaths;
        dStdDevPrice -= dPrice * dPrice;
        
        std::cout << "Price = " << dPrice << std::endl;
        std::cout << "StdDev = " << dStdDevPrice << std::endl;

        std::cout<<"Total Time elapsed : " << (double)(clock()-start)/CLOCKS_PER_SEC <<" sec"<< std::endl;
    }
    else if (iChoice == 76)
    {
        std::vector<double> dVect(1,0.0);
        Test::TestVectorDaughter sTestDaughter(dVect);
        sTestDaughter.DoSomething();
    }
    double dMean = 0.0, dStdDev = 0.0;
    
    iNRealisations = dRealisations.size();
    if (iNRealisations > 0)
    {
        for (std::size_t iRealisation = 0 ; iRealisation < iNRealisations ; ++iRealisation)
        {
            dMean += dRealisations[iRealisation];
            dStdDev += dRealisations[iRealisation] * dRealisations[iRealisation];
        }
        dMean /= iNRealisations;
        dStdDev /= iNRealisations;
        
        dStdDev -= dMean * dMean;
        dStdDev = sqrt(dStdDev);
        
        std::cout << "Mean : " << dMean << std::endl;
        std::cout << "StdDev : " << dStdDev << std::endl;
    }
    // To print something in a file

    /*FILE *pFile;
    pFile = fopen("/Users/alexhum49/Desktop/Myfile.txt", "w"); // open in writing mode
    if (pFile == NULL)
    {
        std::cout << "Error in opening the file" << std::endl;
    }
    try {
        fprintf(pFile, "%.5lf %.5lf", dMean, dStdDev);
    } catch (std::exception & ) {
        
    }
    fclose(pFile);*/
    
    return 0;
}


