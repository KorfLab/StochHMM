//
//  pwm.cpp
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

#include "pwm.h"
namespace StochHMM{


positionWeightMatrix::positionWeightMatrix(){
    size=0;
    valType=PROBABILITY;
}

double positionWeightMatrix::scoreMotif(sequence* seq, int startIter){
    
	int sequenceSize=seq.size();
    int stopIter=startIter+size;
    if (stopIter>sequenceSize-1){
        return -INFINITY;
    }
    
    double score(0.0);
    
	for (int i=0;i<size;i++){
        int letter=seq->seq[startIter+i]
            
        // TODO: Need to implement ambiguous scoring here
        score+=log(weightMatrix[i][letter]/background[letter]);
	}
    
	return 1-getPValue(score);
}

double positionWeightMatrix::getPValue(){
    
    
    
    
    return;
}

void positionWeightMatrix::setBackground(vector<double>& bg){
    background.clear();
    background.assign(bg.begin(),bg.end());
    return;
}


void positionWeightMatrix::initializePValue(){
    
    //Compute all possible scores for the matrix
    vector<double>* results = new(std::nothrow) vector<double>(weightMatrix[0]);
    
    if (results==NULL){
        std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
        exit(1);
    }
    
    for(size_t i=1;i<weightMatrix.size();i++){
        vector<double>* temp = new(std::nothrow) vector <double>;
        
        if (temp==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        addVectorCombinatorial(*temp, *results, weightMatrix[i]);
        delete results;
        results=temp;
    }
    
    size_t Nscores=results->size();
    
    //Sort the scores
    sort(results->begin(), results->end());
    
    
    //figure out p-value 
    
    return;
}

}