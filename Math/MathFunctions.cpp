//
//  MathFunctions.cpp
//  FinanceTools
//
//  Created by Alexandre HUMEAU on 04/03/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Constants.h"
#include "MathFunctions.h"
#include "Require.h"
#include <cmath>

namespace MathFunctions {

    double NormalCumulativeDistribution(const double x)
    {
        // return the cumulative normal distribution taken in x.
        if (x < 0.0)
        {
            return 1 - NormalCumulativeDistribution(-x);
        }
        else
        {
            // take the approximation from Abramovitz and Stegun
            // N(x) = 1 - exp(.5(x*x))*(b1 * t+ b2 * t^2 + b3 * t^3 + b4 * t^4 + b5 * t^5) / sqrt(2* PI)
            // b0 = 0.2316419,
            // b1 = 0.319381530,
            // b2 = −0.356563782,
            // b3 = 1.781477937,
            // b4 = −1.821255978,
            // b5 = 1.330274429
            // t = 1 / (1 + b0 * x)
            double b0 = 0.2316419, b1 = 0.319381530, b2 = -0.356563782, b3 = 1.781477937, b4 = -1.821255978, b5 = 1.330274429, t = 1 / (1 + b0 * x);
            return 1 - exp(0.5 * x * x) * (b1 * t + b2 * t*t + b3 * t * t * t + b4 * t * t * t * t + b5 * t * t * t * t * t) / sqrtpi(2.0);
        }
    }

    double sqrtpi(const double x)
    {
        return sqrt(PI * x);
    }

    // More accurate cumulative normal function
    // Implementation in Better Approximations for cumulatives normal functions by Graeme West

    double AccCumNorm(const double x)
    {
        double dRes;
        double dAbsx = fabs(x);
        double dBuild = 0.0;
        if (dAbsx > 37.0)
        {
            dRes = 0.0;
        }
        else
        {
            double dExp = exp(-dAbsx * dAbsx * 0.5);
            if (dAbsx < 7.07106781186547)
            {
                dBuild = 3.52624965998911E-02 * dAbsx + 0.700383064443688;
                dBuild = dBuild * dAbsx + 6.37396220353165;
                dBuild = dBuild * dAbsx + 33.912866078383;
                dBuild = dBuild * dAbsx + 112.079291497871;
                dBuild = dBuild * dAbsx + 221.213596169931;
                dBuild = dBuild * dAbsx + 220.206867912376;

                dRes = dExp * dBuild;

                dBuild = 8.83883476483184E-02 * dAbsx + 1.75566716318264;
                dBuild = dBuild * dAbsx + 16.064177579207;
                dBuild = dBuild * dAbsx + 86.7807322029461;
                dBuild = dBuild * dAbsx + 296.564248779674;
                dBuild = dBuild * dAbsx + 637.333633378831;
                dBuild = dBuild * dAbsx + 793.826512519948;
                dBuild = dBuild * dAbsx + 440.413735824752;

                dRes /= dBuild;
            }
            else
            {
                dBuild = dAbsx + 0.65;
                dBuild = dAbsx + 4.0 / dBuild;
                dBuild = dAbsx + 3.0 / dBuild;
                dBuild = dAbsx + 2.0 / dBuild;
                dBuild = dAbsx + 1.0 / dBuild;
                dRes = dExp / dBuild / 2.506628274631;
            }
        }
        if (x > 0)
        {
            return 1.0 - dRes;
        }
        else
        {
            return dRes;
        }
    }

    // More accurate bivariate cumulative normal function
    // Implementation in Better Approximations for cumulatives normal functions by Graeme West

    double AccBivarCumNorm(const double x, const double  y, const double rho)
    {
        Utilities::require(-1.0 < rho < 1.0);
        int i;
        double dRes;
        double a[] = {0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992};
        double W[] ={0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042};
        double h1, h2, LH = 0.0, h12, h3, h5, h6, h7, r1, r2, r3, rr, AA, ab;
        h1 = x;
        h2 = y;
        h12 = 0.5 * (h1 * h1 + h2 * h2);
        if (fabs(rho) > 0.7)
        {
            r2 = 1 - rho * rho;
            r3 = sqrt(r2);

            if (rho < 0.0)
            {
                h2 = -h2;
            }
            h3 = h1 * h2;
            h7 = exp(-0.5 * h3);
            if (fabs(rho) < 1.0)
            {
                h6 = fabs(h1 - h2);
                h5 = 0.5 * h6 * h6;
                h6 /= r3;
                AA = 0.5 - 0.125 * h3;
                ab = 3.0 - 2.0 * AA * h5;
                LH = 0.13298076 * h6 * ab * (1 - AccCumNorm(h6)) - exp(-h5 / r2) * (ab + AA * r2) * 0.053051647;

                for (i = 0 ; i < 5 ; ++ i)
                {
                    r1 = r3 * a[i];
                    rr = r1 * r1;
                    r2 = sqrt(1 - rr);
                    LH -= W[i] * exp(-h5 / rr) * (exp(-h3 / (1.0 + r2)) / r2 / h7 - 1.0 - AA * rr);
                }
            }
            dRes = LH * r3 * h7 + AccCumNorm(std::min(h1,h2));
            if (rho < 0.0)
            {
                dRes = AccCumNorm(h1) - dRes;
            }
        }
        else
        {
            h3 = h1 * h2;
            if (rho != 0.0)
            {
                for (i = 0 ; i < 5 ; ++ i)
                {
                    r1 = rho * a[i];
                    r2 = 1 - r1 * r1;
                    LH += W[i] * exp((r1 * h3 - h12) / r2) / sqrt(2.0);
                }
            }
            dRes = AccCumNorm(h1) * AccCumNorm(h2) - LH * rho;
        }
        return dRes;
    }

    double InvCumNorm(const double p)
    {
        Utilities::require(0.0 <= p <= 1.0);
        // Implementation comes from http://home.online.no/~pjacklam/notes/invnorm/

        //Coefficients in rational approximations.

        double a[] = {0.0, -3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00},
        b[] = {0.0, -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01},
        c[] = {0.0, -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00,  4.374664141464968e+00, 2.938163982698783e+00},
        d[] = {0.0, 7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};

        //Define break-points.

        double p_low  = 0.02425, p_high = 1 - p_low;
        double q,x = 0,r;

        //Rational approximation for lower region.

        if (0 < p && p < p_low)
        {
            q = sqrt(-2*log(p));
            x = (((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6]) / ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1);
        }

        // Rational approximation for central region.

        else if (p_low <= p && p <= p_high)
        {
            q= p - 0.5;
            r = q*q;
            x = (((((a[1]*r+a[2])*r+a[3])*r+a[4])*r+a[5])*r+a[6])*q / (((((b[1]*r+b[2])*r+b[3])*r+b[4])*r+b[5])*r+1);
        }
        // Rational approximation for upper region.

        else if (p_high < p && p < 1)
        {
            q = sqrt(-2*log(1-p));
            x = -(((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6]) / ((((d[1]*q+d[2]*q+d[3])*q+d[4])*q+1));
        }

        return x;
    }

    double Beta_OU(const double dLambda, const double dt)
    {
        double x = dLambda * dt;
        if (fabs(x) > BETAOUTHRESHOLD)
        {
            return (1.0 - exp(-x)) / dLambda;
        }
        else
        {
            return x - 0.5 * x * dt + 1 / 6.0 * x * x * dt - 1 / 24.0 * x * x * x * dt;
        }
    }
    
    // returns sum(exp(dLambda * u)du, u=dt1..dt2)
	double SumExp(const double dLambda, const double dt1, const double dt2)
    {
		return Beta_OU(-dLambda, dt2) - Beta_OU(-dLambda, dt1);
    }
	
	// Determinist part for the HW1F model, before factorization with beta(T) - beta(t)
	// = sum[{a^2*(beta(T) + betadt) - 2*beta(s))}ds, s=0..t] / 2
	double DeterministPartZCBHW1F(const double dLambda, const double dt, const double dT, const Finance::TermStructure <double, double> & dTermStructure)
	{
		// Reading the TermStructure of sigma
		std::vector<double> dVariables = dTermStructure.GetVariables() ;
		std::vector<double> dValues = dTermStructure.GetValues() ;
        
        if (dTermStructure.IsTermStructure())
        {
            std::cout << "Not yet implemented" << std::endl;
        }
        else
        {
            double dSigma = dValues[0];
            if (std::abs(dLambda) > BETAOUTHRESHOLD)
            {
                // If lambda is not too small
                return 0.5 * dSigma * dSigma / (dLambda * dLambda) * (4.0 * (exp(dLambda * dt) - 1.0) - (exp(-dLambda * dt) + exp(-dLambda * dT)) *(exp(2.0 * dLambda * dt) - 1.0));
            }
            else
            {
                // Small values of lambda
                std::cout << "To be implemented";
            }
        }
        //  Error value
        return std::numeric_limits<double>::infinity();
        
        //  Old code which was bugged
		/*if (std::abs(dLambda) > BETAOUTHRESHOLD)
        {			
			// if dLambda isnt too small
			while (dTValues[k] <= dt && k+1 < SizeS) {
				s2 = dSValues[k] * dSValues[k] ;
				sum1 = exp(dLambda * dTValues[k]) + exp(dLambda * std::min(dTValues[k+1],dt)) - 2.0 ;
				sum1 /= dLambda ;
				sum1 *= SumExp(-2.0*dLambda, dTValues[k], std::min(dTValues[k+1],dt)) ;
				sum1 -= (2.0 / dLambda * (SumExp(-dLambda, dTValues[k], std::min(dTValues[k+1],dt)) - (std::min(dTValues[k+1],dt) - dTValues[k])))  ;
				sum1 *= s2 ;
				sum2 += sum1 ;
				++k ;
			}
            
			if (dTValues[k] <= dt && k+1 == SizeS) {
				s2 = dSValues[k] * dSValues[k] ;
				sum1 = exp(dLambda * dTValues[k]) + exp(dLambda * dt) - 2.0 ;
				sum1 /= dLambda ;
				sum1 *= SumExp(-2.0*dLambda, dTValues[k], dt) ;
				sum1 -= (2.0 / dLambda * (SumExp(-dLambda, dTValues[k], dt) - (dt - dTValues[k])))  ;
				sum1 *= s2 ;
				sum2 += sum1 ;
            }
			return -0.5 * sum2 ;
		}
		
		else
        {
            // For small values of dLambda
			while (dTValues[k] <= dt && k+1 < SizeS) {
				s2 = dSValues[k] * dSValues[k] ;
				sum1 = (std::min(dTValues[k+1],dt) - dTValues[k]) * (-(dT + dt) * (std::min(dTValues[k+1],dt) + dTValues[k])
														+ (dT * dT + dt * dt) * 0.5
														- (dTValues[k]* dTValues[k] + std::min(dTValues[k+1],dt) * std::min(dTValues[k+1],dt) - std::min(dTValues[k+1],dt) * dTValues[k]) / 3) ;
				sum1 *= dLambda ;
				sum1 += (std::min(dTValues[k+1],dt) - dTValues[k]) * (dT + dt + dTValues[k] + std::min(dTValues[k+1],dt)) ;
				sum1 *= s2 ;
				sum2 += sum1 ;
				++k ;
			}
            
			if (dTValues[k] <= dt && k+1 == SizeS) {
				s2 = dSValues[k] * dSValues[k] ;
				sum1 = (dT - dTValues[k]) * ((dT + dt) * (dT + dTValues[k])
											 + (dT * dT + dt * dt) * 0.5
											 - (dTValues[k] * dTValues[k] + dT * dT - dT * dTValues[k]) / 3) ;
				sum1 *= dLambda ;
				sum1 += (dT - dTValues[k]) * (dT + dt + dTValues[k] + dT) ;
				sum1 *= s2 ;
				sum2 += sum1 ;
            }

			return -0.5 * sum2 ;
        }*/
	}

    double BlackScholes(const double dForward, const double dStrike, const double dStdDev, const int iPhi)
    {
        // iPhi = 1 iff it is a call price else it is a put price
        Utilities::require((iPhi == 1) || (iPhi ==-1));
        Utilities::require(dForward > 0.0);
        Utilities::require(dStrike > 0.0);
        Utilities::require(dStrike > 0.0);
        
        double d1 = log(dForward / dStrike) / dStdDev + 0.5 * dStdDev, d2 = d1 - dStdDev;
        
        return iPhi * (dForward * AccCumNorm(iPhi * d1) - dStrike * AccCumNorm(iPhi * d2));
    }

	// 
	double QuasiMC2(const std::size_t iNum, const std::size_t iBase) {
		// Calculation of the max power of iBase in the decomposition, it is an 
		std::size_t iSize = log((double)iNum) / log((double)iBase) ;

		// Initialization of iBase^0 = 1
		std::size_t iPowBase = 1 ;

		// iCoeff will contain the iBase^i * (int)( iNum / iBase^i), the coefficient iCoeffNumBase of the decomposition are the iCoeff(i) - iCoeff(i+1)
		// It is initialized for i = 0
		std::size_t iCoeff = iNum ;

		// Initialization of the coefficients divided by iPowBase si it is a DOUBLE
		double iCoeffNumBase = 0;

		double dQuasiMC2 = 0 ;

		for (std::size_t i = 0; i < iSize + 1; i++) {
			iCoeffNumBase = iCoeff ;
			// Increase the power of iBase
			iPowBase *= iBase ;
			// Calculation of the next iCoeff
			iCoeff = iPowBase * (std::size_t)(iNum / iPowBase) ;
			iCoeffNumBase -= iCoeff ;
			iCoeffNumBase /= iPowBase ;
			dQuasiMC2 += iCoeffNumBase ;
		}

		return dQuasiMC2 ;
	}

}