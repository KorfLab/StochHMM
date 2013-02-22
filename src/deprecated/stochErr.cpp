//
//  StochError.cpp

//Copyright (c) 2007-2012 Paul C Lott 
//University of California, Davis
//Genome and Biomedical Sciences Facility
//UC Davis Genome Center
//Ian Korf Lab
//Website: www.korflab.ucdavis.edu
//Email: lottpaul@gmail.com
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of
//this software and associated documentation files (the "Software"), to deal in
//the Software without restriction, including without limitation the rights to
//use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
//the Software, and to permit persons to whom the Software is furnished to do so,
//subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
//FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
//COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
//IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "stochErr.h"

namespace StochHMM{
    
    
    //!Set the time of the error
    void stochErrInfo::setTime(){
        time_t rawtime;
        time(&rawtime);
        sTime=ctime(&rawtime);
        sTime=sTime.substr(0,24);
    }
    
    //!Get Error information as string
    std::string stochErrInfo::stringify(){
        std::string output;
        std::stringstream out;
        
        output+="Time:\t" + sTime + "\n";
        output+="File:\t" + sFile + "\n";
        output+="Function:\t" + sFunction + "\n";
        
        out << sLine;
        
        output+="Line:\t" + out.str() + "\n";
        if (!sInfo.empty()){
            output+="Additional Error Information:\t" + sInfo + "\n";
        }
        
        
        return output;
    }
    
    
    //! Print error information to stderr
    void stochErrInfo::print(){
        std::cerr << stringify() << std::endl;
    }
    
    //!Create a stochErr
    //!Error doesn't contain any additional information
    stochErr::stochErr(const errStoch err, const char* func, const char* file, int ln) : sErr(err)
    {
        stochErrInfo* sinfo = new stochErrInfo(func,file,ln);
        sInfo.push_back(sinfo);
    }
    
    
    //!Create a stochErr
    //!Error: contains additional information
    //!\param err Type of error
    //!\param infor Error information
    //!\param func Function where error occurred
    //!\param file File where error occurred
    //!\param ln Line where error occurred
    stochErr::stochErr(const errStoch err, const char* info, const char* func, const char* file, int ln) : sErr(err)
    {
        stochErrInfo* sinfo = new stochErrInfo(info,func,file,ln);
        sInfo.push_back(sinfo);
    }
    
    //    stochErr::stochErr(stochErr& err){
    //        sErr = err.sErr;
    //        sPriority = err.sPriority;
    //        sInfo = err.sInfo;
    //    }
    
    //!Destroy the stochErr
    stochErr::~stochErr() throw(){
        if (sInfo.size()!=0){
            for(size_t i=0;i<sInfo.size();i++){
                delete sInfo[i];
                sInfo[i] = NULL;
            }
        }
        sErr = sNoError;
    }
    
    //!Get readable error priority
    //! \return const char* Error Priority
    const char* stochErr::getPriority() const{
        return priority[sPriority];
    }
    
    //!Get readable Error call
    //! \return const char* Error that occured
    const char* stochErr::getError() const{
        return ErrInfo[sErr];
    }
    
    
    //!Get readable error call and information
    //!Will return the whole stack throw information
    const char* stochErr::what() const throw(){
        std::string output;
        output= "Thrown Error:\t";
        output+= getError();
        output+="\nError Information:\n";
        
        for(size_t i=0;i<sInfo.size();i++){
            output+=sInfo[i]->stringify();
        }
        return output.c_str();
    }
    
}