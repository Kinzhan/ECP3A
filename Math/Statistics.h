//
//  Statistics.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 22/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_Statistics_h
#define Seminaire_Statistics_h

#include <vector>

namespace Stats {
    class Statistics
    {
    public:
        Statistics();
        virtual ~Statistics();
        
        //  Method to compute the mean of data
        virtual double Mean(const std::vector<double> & dData) const;
        
        //  Method to compute the median of data
        virtual double Median(const std::vector<double> & dData) const;
        
        //  Method to compute the standard deviation of data
        virtual double StandardDeviation(const std::vector<double> & dData) const;
        
        //  Method to compute the variance of data
        virtual double Variance(const std::vector<double> & dData) const;
        
        //  Method to compute the Quantile of data
        virtual double Quantile(const double dQuantile, const std::vector<double> & dData) const;
        
        //  Method to compute the empirical distribution of the Data given a fixed number of buckets
        virtual std::vector<std::pair<double,std::size_t> > EmpiricalDistribution(const std::vector<double> & dData, const std::size_t iNBuckets) const;
        
        //  Method to compute the empirical distribution of the data given fixed points
        virtual std::vector<std::pair<double,std::size_t> > EmpiricalDistribution(const std::vector<double> & dData, const std::vector<double> & dPoints) const;
    };
}

#endif
