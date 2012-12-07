//
//  SimulationData.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 24/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_Header_h
#define Seminaire_Header_h

#include <iostream>
#include <map>
#include <cmath> // for floor and pow
#include "VectorUtilities.h"
#include "StringUtilities.h"

namespace Finance{
    
    class SimulationData
    {
    public:
        //           Dates       Paths
        typedef std::vector<std::vector<std::vector<double> > > Cube;
        typedef std::pair<std::vector<double>, Cube> Data;
        
        SimulationData()
        {};
        
        virtual ~SimulationData()
        {
            DateList_.clear();
            Data_.first.clear();
            Data_.first.~vector();
            
            for (Cube::iterator itDates = Data_.second.begin() ; itDates != Data_.second.end() ; ++itDates)
            {
                for (std::vector<std::vector<double> >::iterator itPaths = (*itDates).begin() ; itPaths != (*itDates).end() ; ++itPaths)
                {
                    (*itPaths).clear();
                }
            }
        };
        
        void ReadFromFile(const char * cFile)
        {
            std::ifstream stream;
            stream.open(cFile);
            std::string cDate, cPath, cValue, cValues;
            stream >> cDate >> cPath >> cValue;
            long lDate, lPath;
            
            while (!stream.eof())
            {
                //  Get a string and then split the string
                stream >> lDate >> lPath >> cValues;
                std::size_t iDate;
                if (!Utilities::IsFound(DateList_, lDate, &iDate))
                {
                    DateList_.push_back(lDate);
                }
                
                //  Split vectorValues
                std::vector<std::string> cVectorValues = Utilities::Split(cValues, " ");
                
                Data_.first = std::vector<double>(1,0.0);
                std::vector<double> dVectorValues;
                for (std::size_t i = 0 ; i < cVectorValues.size() ; ++i)
                {
                    //  Convert the string values that represent doubles to doubles
                    dVectorValues.push_back(atof(cVectorValues[i].c_str()));
                }
                //  Put the vector in the Data_ map
                Put(iDate, lPath, dVectorValues);
            }
        };
        
        void ReadFromFile(const std::string & cFile)
        {
            ReadFromFile(cFile.c_str());
        };
        
        void PrintInFile(const char * cFile, bool bAppend) const
        {
            FILE* sFile = fopen(cFile, bAppend ? "a" : "w");
            if (sFile)
            {
                fprintf(sFile, "Date Path Value\n");
                std::size_t iDate = 0;
                for (Cube::const_iterator itDates = Data_.second.begin() ; itDates != Data_.second.end() ; ++itDates)
                {
                    std::size_t iPath = 0;
                    for (std::vector<std::vector<double> >::const_iterator itPaths = (*itDates).begin() ; itPaths != (*itDates).end() ; ++itPaths)
                    {
                        fprintf(sFile, "%ld %lu", DateList_[iDate], iPath);
                        for (std::size_t i = 0 ; i < (*itPaths).size() ; ++i)
                        {
                            fprintf(sFile, " %.10lf", (*itPaths)[i]);
                        }
                        fprintf(sFile, "\n");
                        //increment iPath
                        iPath++;
                    }
                    //  increment iDate
                    iDate++;
                }
            }
            fclose(sFile);
        }
        
        //  Method to put values for a date at a specific path in Data_
        void Put(const std::size_t iDate, const std::size_t iPath, const std::vector<double> & dValues)
        {
            if (iDate < Data_.second.size())
            {
                if (iPath < Data_.second[iDate].size())
                {
                    Data_.second[iDate][iPath] = dValues;
                }
                else
                {
                    std::size_t iNPaths = Data_.second[iDate].size();
                    for (std::size_t i = iNPaths ; i < iPath + 1 ; ++i)
                    {
                        std::vector<double> dValues0(1,0.0);
                        Data_.second[iDate].push_back(dValues0);
                    }
                    Data_.second[iDate].back() = dValues;
                }
            }
            else
            {
                std::size_t iNDates = Data_.second.size();
                for (std::size_t i = iNDates ; i < iDate + 1 ; ++i)
                {
                    std::vector<std::vector<double> > dVectorValues(iPath + 1, std::vector<double>(1,0.0));
                    Data_.second.push_back(dVectorValues);
                    Data_.second[iDate].back() = dValues;
                }
            }
        }
        
        //  Method to set all the dates
        void SetDates(const std::vector<double> & dDates)
        {
            DateList_.resize(dDates.size());
            for (std::size_t i = 0 ; i < dDates.size() ; ++i)
            {
                DateList_[i] = (int)floor(dDates[i] * 365);
            }
        }
        
        void SetDates(const std::vector<long> & lDates)
        {
            DateList_ = lDates;
        }
        
        //  Get an element from the Data
        std::vector<double> Get(const std::size_t iDate, const std::size_t iPath)
        {
            //if (Data_.second.count(iDate) != 0)
            if (iDate < Data_.second.size())
            {
                if (iPath < Data_.second[iDate].size())
                {
                    return Data_.second[iDate][iPath];
                }
            }
            else
            {
                std::cout<<"Date "<< DateList_[iDate] << " not found"<<std::endl;
            }
            //  return an error value
            return std::vector<double>(1, pow(10.0,10));
        }
        
        Data GetData() const
        {
            return Data_;
        }
        
        std::vector<long> GetDateList() const
        {
            return DateList_;
        }
        
    private:
        Data Data_;
        std::vector<long> DateList_;
    };
}

#endif
