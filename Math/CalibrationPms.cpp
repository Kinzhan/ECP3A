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
    
    NewtonFunctionMu0::NewtonFunctionMu0(double dDeltaT, const std::vector<double> & dData, double dEpsilonDerivative) : dDeltaT_(dDeltaT), dData_(dData), dEpsilonDerivative_(dEpsilonDerivative)
    {}
    
    double NewtonFunctionMu0::func(double dLambda)
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
    
    double NewtonFunctionMu0::dfdlambda(double dLambda)
    {
        //  Centered approximation of the derivative
        return (func(dLambda + dEpsilonDerivative_) - func(dLambda - dEpsilonDerivative_)) * 0.5 / dEpsilonDerivative_;
    }
    
    NewtonFunction::NewtonFunction(double dDeltaT, const std::vector<double> & dData, double dEpsilonDerivative) : dDeltaT_(dDeltaT), dData_(dData), dEpsilonDerivative_(dEpsilonDerivative)
    {}
    
    double NewtonFunction::VolatilitySquare(double dLambda)
    {
        double dMu = 0;
        double dSum = 0;
        std::size_t iN = dData_.size();
        
        for (std::size_t i = 0 ; i < iN - 1 ; ++i)
        {
            dMu += dData_[i + 1] - dData_[i] * exp(-dLambda * dDeltaT_);
        }
        dMu /= iN;
        for (std::size_t i = 0 ; i < iN - 1 ; ++i)
        {
            double dDiff = dData_[i + 1] - dData_[i] * exp(-dLambda * dDeltaT_) - dMu ;
            dSum += dDiff * dDiff;
        }
        return 2. * dLambda / ((1 - exp(-2. * dLambda * dDeltaT_)) * iN) * dSum;
    }
    
    double NewtonFunction::Mu(double dLambda)
    {
        double dMu = 0 ; 
        std::size_t iN = dData_.size();
        for (std::size_t i = 0 ; i < iN ; ++i)
        {
            dMu += dData_[i + 1] - dData_[i] * exp(-dLambda * dDeltaT_);
        }
        return 1. / (iN * (1 - exp(-2 * dLambda * dDeltaT_))) * dMu;
    }
    
    double NewtonFunction::func(double dLambda)
    {
        std::size_t iN = dData_.size();
        double dSum = 0, dMu = Mu(dLambda);
        for (std::size_t i = 0 ; i < iN - 1 ; ++i)
        {
            dSum += (dData_[i] - dMu) * (dData_[i + 1] - dData_[i] * exp(-dLambda * dDeltaT_) - dMu * (1 - exp(-dLambda * dDeltaT_)));
        }
        double dResult = 0;
        dResult -= iN * dDeltaT_ * exp(-2 * dLambda * dDeltaT_) / (1 - exp(-2 * dLambda * dDeltaT_));
        dResult -= iN / (2. * dDeltaT_);
        dResult += (4. * dLambda * dDeltaT_ * exp(-2. * dLambda * dDeltaT_) - 2 * (1 - exp(-2. * dLambda * dDeltaT_))) / ((1 - exp(-2. * dLambda * dDeltaT_)) * (1 - exp(-2. * dLambda * dDeltaT_))) * iN / (2. * dLambda);
        
        dResult += 4. * dLambda * dDeltaT_ * exp(-dLambda * dDeltaT_) / ((1 - exp(-2. * dLambda * dDeltaT_)) * (1 - exp(-2. * dLambda * dDeltaT_)) * VolatilitySquare(dLambda)) * dSum;
        
        return dResult;
    }
    
    double NewtonFunction::dfdlambda(double dlambda)
    {
        //  Centered approximation of the derivative
        return (func(dlambda + dEpsilonDerivative_) - func(dlambda - dEpsilonDerivative_)) * 0.5 / dEpsilonDerivative_;
    }
    
    CalibrationPms::CalibrationPms(double dLambda)
    {
        //  initialize value for Newton-Raphson Algorithm
        dLambda_ = dLambda;
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
    
    double CalibrationPms::GetMu() const
    {
        return dMu_;
    }
    
    std::vector<double> CalibrationPms::LoadDataFromFile(const std::string &cFileName, bool bIsLibor, double dDeltaT) const
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
                    if (bIsLibor)
                    {
                        dResults.push_back(1 / dDeltaT * log(1 + dDeltaT * dValue));
                    }
                    else
                    {
                        dResults.push_back(dValue);
                    }
                }
            }
            return dResults;
        }
        else 
        {
            throw "Could not open file : " + cFileName;
        }
    }
    
    void CalibrationPms::ComputeSigmaMu0(const std::vector<double> & dData, double dDeltaT)
    {
        double dSum = 0;
        std::size_t iNDataSize = dData.size();
        for (std::size_t i = 0 ; i < iNDataSize - 1 ; ++i)
        {
            dSum += (dData[i+1] - dData[i] * exp(-dLambda_ * dDeltaT)) * (dData[i+1] - dData[i] * exp(-dLambda_ * dDeltaT));
        }
        dSigma_ = 2 * dLambda_ / (iNDataSize * (1 - exp(-2. * dLambda_ * dDeltaT))) * dSum;
        dSigma_ = sqrt(dSigma_);
    }
    
    void CalibrationPms::ComputeMu(const std::vector<double> &dData, double dDeltaT)
    {
        double dSum = 0; 
        std::size_t iN = dData.size();
        for (std::size_t i = 0 ; i < iN - 1 ; ++i)
        {
            dSum += (dData[i + 1] - dData[i] * exp(-dLambda_ * dDeltaT));
        }
        dMu_ = dSum / (iN * (1 - exp(-dLambda_ * dDeltaT)));
    }
    
    void CalibrationPms::ComputeSigma(const std::vector<double> &dData, double dDeltaT)
    {
        ComputeMu(dData, dDeltaT);
        double dSum = 0;
        std::size_t iN = dData.size();
        for (std::size_t i = 0 ; i < iN - 1 ; ++i)
        {
            double dDiff = dData[i + 1] - dData[i] * exp(-dLambda_ * dDeltaT) - dMu_ * (1 - exp(-dLambda_ * dDeltaT));
            dSum += dDiff * dDiff;
        }
        dSigma_ = 2. * dLambda_ / ((1 - exp(-2. * dLambda_ * dDeltaT)) * iN) * dSum;
        dSigma_ = sqrt(dSigma_);
    }
    
    void CalibrationPms::NewtonRaphsonAlgorithmMu0(const std::vector<double> &dData, double dDeltaT, const NewtonPms & sNewtonPms)
    {
        NewtonFunctionMu0 sNewtonFunction(dDeltaT, dData);
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
            ComputeSigmaMu0(dData, dDeltaT);
            
            std::cout << iIter << ";" << dLambda_ << ";" << dSigma_ << std::endl;
            iIter ++;
        }
    }
    
    void CalibrationPms::NewtonRaphsonAlgorithm(const std::vector<double> &dData, double dDeltaT, const Calibration::NewtonPms &sNewtonPms)
    {
        NewtonFunction sNewtonFunction(dDeltaT, dData);
        std::size_t iIter = 0;
        double dLambdaOld = dLambda_ / 2;
        std::cout << "Iter ; Lambda ; Sigma ; Mu ; f(Lambda)" << std::endl;
        while (iIter <= sNewtonPms.iNIterMax_ && std::abs(sNewtonFunction.func(dLambda_)) > sNewtonPms.dTolerance_)
        {
            dLambdaOld = dLambda_;
            double dDerivativeValue = sNewtonFunction.dfdlambda(dLambdaOld),
            dfValue = sNewtonFunction.func(dLambdaOld);
            dLambda_ -= dfValue / dDerivativeValue ;
            
            //  compute sigma at each step of the algorithm
            ComputeSigma(dData, dDeltaT);
            ComputeMu(dData, dDeltaT);
            
            std::cout << iIter << ";" << dLambda_ << ";" << dSigma_ << ";" << dMu_ << ";" << sNewtonFunction.func(dLambda_) << std::endl;
            iIter ++;
        }
    }
}