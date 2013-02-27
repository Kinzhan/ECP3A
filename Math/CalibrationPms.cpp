//
//  CalibrationPms.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 23/02/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>
#include "CalibrationPms.h"

namespace Calibration {
    
    NewtonPms::NewtonPms(double dTolerance, std::size_t iNIterMax) : dTolerance_(dTolerance), iNIterMax_(iNIterMax)
    {}
    
    NewtonFunction::NewtonFunction(double dDeltaT, const std::vector<double> & dData, double dEpsilonDerivative) : dDeltaT_(dDeltaT), dData_(dData), dEpsilonDerivative_(dEpsilonDerivative)
    {}
    
    double NewtonFunction::func(double dLambda)
    {
        double dSum = 0., dSumSquare = 0.;
        for (std::size_t i = 0 ; i < dData_.size() ; ++i)
        {
            double dDiff = (dData_[i + 1] - dData_[i] * exp(-dLambda * dDeltaT_));
            dSum += dData_[i] * dDiff;
            dSumSquare += dDiff * dDiff;
        }
        
        double d1MinusExp = 1 - exp(-2. * dLambda * dDeltaT_);
        
        return 2. * dDeltaT_ * exp(-dLambda * dDeltaT_) / d1MinusExp - 1 / dLambda + d1MinusExp / (dLambda * dSumSquare) * (2 * dLambda * dDeltaT_ * exp(-dLambda * dDeltaT_) / d1MinusExp * dSum + (d1MinusExp - dLambda * 0.5 / dDeltaT_) / (d1MinusExp * d1MinusExp) * dSumSquare);
    }
    
    double NewtonFunction::dfdlambda(double dLambda)
    {
        //  Centered approximation of the derivative
        return (func(dLambda + dEpsilonDerivative_) - func(dLambda - dEpsilonDerivative_)) * 0.5 / dEpsilonDerivative_;
    }
    
    CalibrationPms::CalibrationPms()
    {
        //  initialize value for Newton-Raphson Algorithm
        dLambda_ = 0.05;
        dSigma_ = 0.01;
    }
    
    CalibrationPms::~CalibrationPms()
    {}
    
    double CalibrationPms::GetLambda() const
    {
        return dLambda_;
    }
    
    double CalibrationPms::GetSigma() const
    {
        return dSigma_;
    }
    
    std::vector<double> CalibrationPms::LoadDataFromFile(const std::string &cFileName) const
    {
        std::ifstream sFile;
        sFile.open(cFileName.c_str());
        if (sFile.is_open())
        {
            std::vector<double> dResults;
            std::string cLine;
            while (std::getline(sFile, cLine)) 
            {
                std::istringstream iss(cLine);
                double dValue;
                while (iss >> dValue)
                {
                    dResults.push_back(dValue);
                }
            }
            return dResults;
        }
        else 
        {
            throw "Could not open file : " + cFileName;
        }
    }
    
    void CalibrationPms::ComputeSigma(const std::vector<double> & dData, double dDeltaT)
    {
        double dSum = 0;
        std::size_t iNDataSize = dData.size();
        for (std::size_t i = 0 ; i < iNDataSize - 1 ; ++i)
        {
            dSum += (dData[i+1] - dData[i] * exp(-dLambda_ * dDeltaT)) * (dData[i+1] - dData[i] * exp(-dLambda_ * dDeltaT));
        }
        dSigma_ = 2 * dLambda_ / (iNDataSize * (1 - exp(-2. * dLambda_ * dDeltaT))) * dSum;
    }
    
    void CalibrationPms::NewtonRaphsonAlgorithm(const std::vector<double> &dData, double dDeltaT, const NewtonPms & sNewtonPms)
    {
        NewtonFunction sNewtonFunction(dDeltaT, dData);
        std::size_t iIter = 0;
        double dLambdaOld = dLambda_ / 2;
        std::cout << "Iter ; Lambda ; Sigma" << std::endl;
        while (iIter <= sNewtonPms.iNIterMax_ && std::abs(dLambda_ - dLambdaOld) > sNewtonPms.dTolerance_)
        {
            dLambdaOld = dLambda_;
            double dDerivativeValue = sNewtonFunction.dfdlambda(dLambdaOld),
                    dfValue = sNewtonFunction.func(dLambdaOld);
            dLambda_ -= dfValue / dDerivativeValue ;
            
            //  compute sigma at each step of the algorithm
            ComputeSigma(dData, dDeltaT);
            
            std::cout << iIter << ";" << dLambda_ << ";" << dSigma_ << std::endl;
            iIter ++;
        }
    }
}