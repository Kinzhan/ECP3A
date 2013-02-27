//
//  PrintInFile.cpp
//  TP
//
//  Created by Alexandre HUMEAU on 22/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <sstream>
#include "PrintInFile.h"

namespace Utilities {
    PrintInFile::PrintInFile() : bAppend_(false)
    {}
    
    PrintInFile::PrintInFile(const std::string & cFileName, const bool bAppend, const std::size_t iPrecision) : cFileName_(cFileName), bAppend_(bAppend), iPrecision_(iPrecision)
    {
        std::stringstream out ;
        out << iPrecision_;
        cPrecision_ = out.str();
    }
    
    PrintInFile::~PrintInFile()
    {
        cFileName_.clear();
    }
    
    void PrintInFile::PrintDataInFile(const std::vector<double> &dData)
    {
        FILE * sFile;
        sFile = fopen(cFileName_.c_str(), bAppend_ ? "a" : "w");
        //  Check if the file is well opened
        if (sFile)
        {
            for (std::size_t i = 0 ; i < dData.size() ; ++i)
            {
                //  Print Data in the file with the given precision
                fprintf(sFile, ("%."+ cPrecision_ + "lf").c_str(), dData[i]);
                fprintf(sFile, "\n");
            }
            // Close the file
            fclose(sFile);
        }
    }
    
    void PrintInFile::PrintDataInFile(const std::vector<std::pair<double, std::size_t> > &dData)
    {
        FILE * sFile;
        sFile = fopen(cFileName_.c_str(), bAppend_ ? "a" : "w");
        //  Check if the file is well opened
        if (sFile)
        {
            for (std::size_t i = 0 ; i < dData.size() ; ++i)
            {
                //  Print Data in the file with the given precision
                fprintf(sFile, ("%."+ cPrecision_ + "lf,%lu").c_str(), dData[i].first, dData[i].second);
                fprintf(sFile, "\n");
            }
            // Close the file
            fclose(sFile);
        }
    }
    
    void PrintInFile::PrintDataInFile(const std::vector<std::pair<double, double> > &dData)
    {
        FILE * sFile;
        sFile = fopen(cFileName_.c_str(), bAppend_ ? "a" : "w");
        //  Check if the file is well opened
        if (sFile)
        {
            for (std::size_t i = 0 ; i < dData.size() ; ++i)
            {
                //  Print Data in the file with the given precision
                fprintf(sFile, ("%."+ cPrecision_ + "lf,%."+ cPrecision_ + "lf").c_str(), dData[i].first, dData[i].second);
                fprintf(sFile, "\n");
            }
            // Close the file
            fclose(sFile);
        }
    }
    
    void PrintInFile::PrintDataInFile(const std::vector<std::vector<double> > & dData)
    {
        FILE * sFile;
        sFile = fopen(cFileName_.c_str(), bAppend_ ? "a" : "w");
        //  Check if the file is well opened
        if (sFile)
        {
            for (std::size_t i = 0 ; i < dData.size() ; ++i)
            {
                //  Print Data in the file with the given precision
                //fprintf(sFile, ("%."+ cPrecision_ + "lf,%."+ cPrecision_ + "lf").c_str(), dData[i].first, dData[i].second);
                for (std::size_t j = 0 ; j < dData[i].size() ; ++j)
                {
                    fprintf(sFile, ("%."+ cPrecision_ + "lf,").c_str(), dData[i][j]);
                }
                fprintf(sFile, "\n");
            }
            // Close the file
            fclose(sFile);
        }
    }
    
    void PrintInFile::PrintDataInFile(const std::map<double, std::map<double, double> > &mData)
    {
        FILE * sFile;
        sFile = fopen(cFileName_.c_str(), bAppend_ ? "a" : "w");
        //  Check if the file is well opened
        if (sFile)
        {
            fprintf(sFile,",");
            std::map<double, std::map<double, double> >::const_iterator iter = mData.begin();
            for (std::map<double, double>::const_iterator it = iter->second.begin() ; it != iter->second.end() ; ++it)
            {
                fprintf(sFile, ("%."+ cPrecision_ + "lf,").c_str(), it->first);
            }
            fprintf(sFile, "\n");
            for ( ; iter != mData.end() ; ++iter)
            {
                fprintf(sFile, ("%."+ cPrecision_ + "lf,").c_str(), iter->first);
                for (std::map<double, double>::const_iterator it = iter->second.begin() ; it != iter->second.end() ; ++it)
                {
                    fprintf(sFile, ("%."+ cPrecision_ + "lf,").c_str(), it->second);
                }
                fprintf(sFile, "\n");
            }
            fclose(sFile);
        }
    }
}