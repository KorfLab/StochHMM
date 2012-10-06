//
//  weight.cpp
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

#include "weight.h"
namespace StochHMM{


    weight::weight(){
        absolute=false;
        absoluteValue=-INFINITY;
        maxValue=NULL;
        minValue=NULL;
        vals=NULL;
    }

    weight::~weight(){
        delete maxValue;
        delete minValue;
        delete vals;
        maxValue=NULL;
        minValue=NULL;
        vals=NULL;
        return;
    }

    void weight::setMaxWeight(double& x, double& y){
        if (maxValue!=NULL){
            delete maxValue;
        }
        
        maxValue=new(std::nothrow) std::pair<double,double>;
        if (maxValue==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        maxValue->first=x;
        maxValue->second=y;
        return;
    }

    void weight::setMinWeight(double& x, double& y){
        if (minValue!=NULL){
            delete minValue;
        }
        
        minValue=new(std::nothrow) std::pair<double,double>;
        if (minValue==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        minValue->first=x;
        minValue->second=y;
        return;
    }

    void weight::setWeights(std::vector<double>& xVals, std::vector<double>& yVals){
        if (maxValue==NULL && minValue==NULL){
            setWeights(xVals, yVals, NULL, NULL);
        }
        else if (maxValue==NULL){
            setWeights(xVals, yVals, minValue, NULL);
        }
        else if (minValue==NULL){
            setWeights(xVals, yVals, NULL, maxValue);
        }
        else{
            setWeights(xVals,yVals,minValue,maxValue);
        }
    }

    //Set weights given the coordinates
    void weight::setWeights( std::vector<double>& xVals, std::vector<double>& yVals, std::pair<double,double>* xMin, std::pair<double,double>* xMax ){
        
        if (vals!=NULL){
            delete vals;
        }
        
        vals=new(std::nothrow) std::vector<std::pair<double,double> >;
        if (vals==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        if (xVals.size()!=yVals.size()){
            std::cerr << "X-Y values are not the same size." <<std::endl;
            return;
        }
        
        for(size_t i=0;i<xVals.size();i++){
            std::pair<double,double> xy (xVals[i],yVals[i]);
            vals->push_back(xy);
        }
        return;
    }


    bool compXval(const std::pair<double,double>& a, const std::pair<double,double>& b){
        return (a.first<b.first);
    }


    //Given a score get the weighted score
    double weight::getWeightedScore(double score){
        //Find score that flanks or is equal to the score we supply
        //If there is one equal to then we return it, otherwise
        //we get the two values and interpolate the score linearly
        //unless we are on the ends then we extrapolate
        
        if (absolute){
            return score + absoluteValue;
        }
        else if (maxValue!=NULL && score > maxValue->first){
            return maxValue->second;
        }
        else if (maxValue!=NULL && score < minValue->first){
            return minValue->second;
        }
        else{
            std::vector<std::pair<double,double> >::iterator it;
            
            std::pair<double,double> temp (score,0);
            
            it=lower_bound(vals->begin(),vals->end(),temp ,compXval);
            
            //std::cout << (*it).first << "\t" << exp((*it).first) <<std::endl;
            
            if ((*it).first==score){
                return (*it).second;
            }
            else if (it==vals->begin()){ //extrapolate based on beginning
                return extrapolate((*it),(*(it+1)),score);
            }
            else if (it==vals->end()){  //extrapolate based on end
                return extrapolate((*(it-2)), (*(it-1)),score);
            }
            else{
                return interpolate((*(it-1)),(*it),score);
            }
        }
    }

    std::string weight::stringify(){
        std::string scaleString;
        
        //Process Values and scaled values;
        scaleString+= "\tVALUE:\tLOG\t[";
        for(size_t i=0;i<vals->size();i++){
            if (i>0){
                scaleString+= ",";
            }
            scaleString+= double_to_string((*vals)[i].first);
        }
        scaleString+= "]\n";
        
        scaleString+= "\tSCALED:\tLOG\t[";
        for(size_t i=0;i<vals->size();i++){
            if (i>0){
                scaleString+= ",";
            }
            scaleString+= double_to_string((*vals)[i].second);
        }
        scaleString+="]\n";
        
        if (minValue!=NULL){
            scaleString+= "\tMIN_VALUE:\tLOG\t" + double_to_string(minValue->first) + "\n";
            scaleString+= "\tMIN_SCALED:\tLOG\t" + double_to_string(minValue->second) + "\n";
        }
        
        if (maxValue!=NULL){
            scaleString+= "\tMAX_VALUE:\tLOG\t" + double_to_string(maxValue->first) + "\n";
            scaleString+= "\tMAX_SCALED:\tLOG\t" + double_to_string(maxValue->second) + "\n";
        }
        
        return scaleString;
    }
    
    bool weight::parse(const std::string& txt){
        stringList lst;
        lst.splitString(txt,"\n");
        
        bool table=false;
        std::vector<double> xVal;
        std::vector<double> yVal;
        
        bool min=false;
        double minXVal;
        double minYVal;
        
        bool max=false;
        double maxXVal;
        double maxYVal;
        
        for(size_t i=0;i<lst.size();i++){
            stringList line;
            clear_whitespace(lst[i],":[]\n ");
            line.splitString(lst[i],"\t,");
            if (line.size()==0){
                std::cerr << "Unable to parse weight line: " << lst[i] << std::endl;
                return false;
            }
                
            if (line[0].compare("SCALE")==0){
                if (line.size()==2){
                    name=line[1];
                }
                else{
                    std::cerr << "Weight Scale is missing a name" << std::endl;
                    return false;
                }
            }
            else if (line[0].compare("ABSOLUTE")==0){
                if (line.size()<2){
                    std::cerr << "Weight Absolute value is missing " << std::endl;
                    return false; 
                }
                
                
                double val;
                if (!stringToDouble(line[1], val)){
                    std::cerr << "Weighted Absoolute value not numeric: " << line[1] << std::endl;
                    return false; 
                }
                                
                setAbsolute(val);
                return true;
            }
            else if (line[0].compare("VALUE")==0){
                
                if (line.size()<3){
                    std::cerr << "Weight VALUE line is missing values" << std::endl;
                    return false; 
                }
                
                bool px=line[1].compare("P(X)")==0;
                for(int i=2;i<line.size();i++){
                    
                    double tempValue;
                    if (!stringToDouble(line[i], tempValue)){
                        std::cerr << "Weighted VALUE not numeric: " << line[i] << std::endl;
                        return false;
                    }
                    
                    if (px){
                        xVal.push_back(log(tempValue));
                    }
                    else{
                        xVal.push_back(tempValue);
                    }
                }
                table=true;
            }
            else if (line[0].compare("SCALED")==0){
                
                if (line.size()<3){
                    std::cerr << "Weight SCALED line is missing values" << std::endl;
                    return false; 
                }
                
                bool px = line[1].compare("P(X)")==0;
                for(int i=2;i<line.size();i++){
                    
                    double tempValue;
                    if (!stringToDouble(line[i], tempValue)){
                        std::cerr << "Weighted SCALED value not numeric: " << line[i] << std::endl;
                        return false;
                    }
                    
                    if (px){
                        yVal.push_back(log(tempValue));
                    }
                    else{
                        yVal.push_back(tempValue);
                    }
                }
                table=true;
            }
            else if (line[0].compare("MIN_VALUE")==0){
                
                if (line.size()<3){
                    std::cerr << "Weight MIN_VALUE line is missing values" << std::endl;
                    return false; 
                }
                
                double tempValue;
                if (!stringToDouble(line[2], tempValue)){
                    std::cerr << "Weighted VALUE not numeric: " << line[i] << std::endl;
                    return false;
                }
                
                
                if (line[1].compare("P(X)")==0){
                    minXVal=log(tempValue);
                }
                else{
                    minXVal=tempValue;
                }
                min=true;
            }
            else if (line[0].compare("MIN_SCALED")==0){
                
                if (line.size()<3){
                    std::cerr << "Weight MIN_SCALED line is missing values" << std::endl;
                    return false; 
                }
                
                double tempValue;
                if (!stringToDouble(line[2], tempValue)){
                    std::cerr << "Weighted MIN_SCALED value not numeric: " << line[i] << std::endl;
                    return false;
                }

                
                
                if (line[1].compare("P(X)")==0){
                    minYVal=log(tempValue);
                    
                }
                else{
                    minYVal=tempValue;
                    
                }
                min=true;
            }
            else if (line[0].compare("MAX_VALUE")==0){
                
                
                if (line.size()<3){
                    std::cerr << "Weight MAX_VALUE line is missing values" << std::endl;
                    return false; 
                }
                
                double tempValue;
                if (!stringToDouble(line[2], tempValue)){
                    std::cerr << "Weighted MAX_VALUE value not numeric: " << line[i] << std::endl;
                    return false;
                }
                

                
                if (line[1].compare("P(X)")==0){
                    maxXVal=log(tempValue);
                }
                else{
                    maxXVal=tempValue;
                }
                max=true;
            }
            else if (line[0].compare("MAX_SCALED")==0){
                
                if (line.size()<3){
                    std::cerr << "Weight MAX_SCALED line is missing values" << std::endl;
                    return false; 
                }
                
                double tempValue;
                if (!stringToDouble(line[2], tempValue)){
                    std::cerr << "Weighted MAX_SCALED value not numeric: " << line[i] << std::endl;
                    return false;
                }
                
                
                if (line[1].compare("P(X)")==0){
                    maxYVal=log(tempValue);
                    
                }
                else{
                    maxYVal=tempValue;
                }
                max=true;
            }
        }
        
        if (table && (yVal.size()==xVal.size())){
            setWeights(xVal, yVal);
        }
        else{
            std::cerr << "Scaling Values entries don't match up. They are of different sizes" <<std::endl;
            return false;
        }
        
        if (min){
            setMinWeight(minXVal, minYVal);
        }
        
        if (max){
            setMaxWeight(maxXVal, maxYVal);
        }
        
        return true;
    }
    
    
    
    void weights::addWeight(weight* wt){
        
        std::string& name=wt->getName();
        if (!wts.count(name)){
            wts[name]=wt;
            return;
        }
        else{
            std::cerr << "Weight of name: " << name << " already exists and is duplicated.  Please remove duplicates.\n";
            exit(1);
        }
    }
    
    weight* weights::operator[](std::string& name){
        if (wts.count(name)){
            return wts[name];
        }
        else{
            return NULL;
        }
    }
    
    void weights::print(){
        std::cout << stringify() << std::endl;
    }
    
    std::string weights::stringify(){
        std::string weightString;
        std::string lnSep(50,'=');
        
        weightString+= "SCALING VALUES\n" + lnSep + "\n";
        std::map<std::string,weight*>::iterator it;
        for(it=wts.begin();it!=wts.end();it++){
            weightString+="SCALE:\t" + (*it).first + "\n";
            weightString+=(*it).second->stringify();
        }
        
        weightString+="\n";
        return weightString;
        
    }
    
    
    
}
