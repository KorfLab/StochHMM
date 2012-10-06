//
//  stochMath.cpp

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

#include "stochMath.h"
namespace StochHMM{

    
    double addLog(double& first, double& second){
        //cout << first <<endl;
        //cout << second <<endl;
        if (first==-INFINITY){
            return second;
        }
        else if (second==-INFINITY){
            return first;
        }
        else if (first>second){
            return first+log(1+exp(second-first));
        }
        else{
            return second+log(1+exp(first-second));
        }
    }
    
    //----------------------------------------------------------------------------//
    // Description: addVectorCombinatorial                                                              
    // Adds the values of vectors combinatorial
    // Example [ 1 4 ] + [ 2 7 ] = [ 2 6 8 11 ] 
    // 
    //----------------------------------------------------------------------------//
    void addVectorCombinatorial(std::vector<int>& result, std::vector<int>& lhs, std::vector<int>& rhs){
        
        if (result.size()>0){
            result.clear();
        }
        
        //If either vector is empty then return copy of the other
        if (lhs.size()==0){
            result.assign(rhs.begin(),rhs.end());
            return;
        }
        else if (rhs.size()==0){
            result.assign(lhs.begin(),lhs.end());
            return;
        }
        
        for(size_t i=0;i<lhs.size();i++){
            for(size_t j=0;j<rhs.size();j++){
                result.push_back(lhs[i]+rhs[j]);
            }
        }
        return;
    }
    
    
    void addVectorCombinatorial(std::vector<double>& result, std::vector<double>& lhs, std::vector<double>& rhs){
        if (result.size()>0){
            result.clear();  //Make sure result is empty
        }
        
        //If either vector is empty then return copy of the other
        if (lhs.size()==0){
            result.assign(rhs.begin(),rhs.end());
            return;
        }
        else if (rhs.size()==0){
            result.assign(lhs.begin(),lhs.end());
            return;
        }
        
        for(size_t i=0;i<lhs.size();i++){
            for(size_t j=0;j<rhs.size();j++){
                result.push_back(lhs[i]+rhs[j]);
            }
        }
        return;
    }
    
    
    void multiplyVectorCombinatorial(std::vector<double>& result, std::vector<double>&lhs, std::vector<double>&rhs){
        if (result.size()>0){
            result.clear();
        }
        
        //If either vector is empty then return copy of the other
        if (lhs.size()==0){
            result.assign(rhs.begin(),rhs.end());
            return;
        }
        else if (rhs.size()==0){
            result.assign(lhs.begin(),lhs.end());
            return;
        }
        
        for(size_t i=0;i<lhs.size();i++){
            for(size_t j=0;j<rhs.size();j++){
                result.push_back(lhs[i]*rhs[j]);
            }
        }
        return;
    }
    
    void multiplyVectorCombinatorial(std::vector<int>& result, std::vector<int>&lhs, std::vector<int>&rhs){
        if (result.size()>0){
            result.clear();
        }
        
        //If either vector is empty then return copy of the other
        if (lhs.size()==0){
            result.assign(rhs.begin(),rhs.end());
            return;
        }
        else if (rhs.size()==0){
            result.assign(lhs.begin(),lhs.end());
            return;
        }
        
        
        for(size_t i=0;i<lhs.size();i++){
            for(size_t j=0;j<rhs.size();j++){
                result.push_back(lhs[i]*rhs[j]);
            }
        }
        return;
    }
    
    
    void addValueToVector(std::vector<int>& vec, int value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]+=value;
        }
        return;
    }
    
    void addValueToVector(std::vector<double>& vec, double value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]+=value;
        }
        return;
    }
    
    
    void multiplyValueToVector(std::vector<double>& vec, double value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]*=value;
        }
        return;
    }
    
    void multiplyValueToVector(std::vector<int>& vec, int value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]*=value;
        }
        return;
    }
    
    
    void divideValueToVector(std::vector<int>& vec, int value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]/=value;
        }
        return;
    }
    
    
    void divideValueToVector(std::vector<double>& vec, double value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]/=value;
        }
        return;
    }
    
    
    
    
    void generateUnknownIndices(std::vector<int>& results, int alphabetSize, int order ,int value){
        results.assign(alphabetSize,0);
        for (int i=0;i<alphabetSize;i++){
            results[i]= value + i * POWER[order-1][alphabetSize-1];
        }
        return;
    }
    
    
    
    double sumVector(std::vector<double>& data){
        double sum=0;
        for(size_t i=0;i<data.size();i++){
            sum+=data[i];
        }
        return sum;
    }
    
    
    
    double minVector(std::vector<double>& data){
        return *min_element(data.begin(), data.end()); 
    }
    
    
    double maxVector(std::vector<double>& data){
        return *max_element(data.begin(), data.end());
    }
    
    
    double avgVector(std::vector<double>& data){
        return sumVector(data) / double(data.size());
    }
    
    void logVector(std::vector<double>& data){
        for(size_t i=0;i<data.size();i++){
            data[i]=log(data[i]);
        }
        return;
    }
    
    void expVector(std::vector<double>& data){
        for(size_t i=0;i<data.size();i++){
            data[i]=exp(data[i]);
        }
        return;
    }
    
    void probVector(std::vector<double>& data){
        double sum=sumVector(data);
        for(size_t iter=0;iter<data.size();iter++){
            if (sum==0){
                data[iter]=0;
            }
            else{
                data[iter]/=sum;
            }
        }
        return;
    }
    
    
    //Linear interpolation of y value given two pair<x,y> and x value
    double interpolate(std::pair<double,double>& a, std::pair<double,double>& b, double& cx){
        //std::cout << a.first << "\t" << a.second <<std::endl;
        //std::cout << b.first << "\t" << b.second <<std::endl;
        
        return a.second+(b.second-a.second)*((cx-a.first)/(b.first-a.first));
    }
    
    //Linear extrapolation of y value given two pair<x,y> and x value
    double extrapolate(std::pair<double,double>& a, std::pair<double,double>& b, double& cx){
        return a.second+((cx-a.first)/(b.first-a.first))*(b.second-a.second);
    }
    
}
