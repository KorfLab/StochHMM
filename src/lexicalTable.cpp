//
//  Lexical.cpp
//  StochHMM
//
//  Created by Paul Lott on 4/2/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#include "lexicalTable.h"



namespace StochHMM{
    
    lexicalTable::lexicalTable(){
        max_order=0;
        
        logProb=NULL;
        counts = NULL;
        prob = NULL;
        
        
        unknownScoreType=NO_SCORE;
        unknownDefinedScore=-INFINITY;
        
        return;
    }
    
    lexicalTable::~lexicalTable(){
        delete logProb;
        delete prob;
        delete counts;
        
        logProb=NULL;
        prob=NULL;
        counts=NULL;
    }
    
    void lexicalTable::createTable(int rows, int columns, int pseudocount, valueType typ){
        if (typ==COUNTS){
            if (counts!=NULL){
                delete counts;
            }
            std::vector<double> temp_columns(columns,pseudocount);
            counts=new(std::nothrow) std::vector<std::vector<double> > (rows,temp_columns);
            
            if (counts==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        return;
    }
    
    void lexicalTable::print() {
        std::cout << stringify() << std::endl;
    }
    
    
    std::vector<std::vector<double> >* lexicalTable::getCountsTable(){
        if (counts==NULL){
            counts = new(std::nothrow) std::vector<std::vector<double> >;
            
            if (counts==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
			
        }
        
        return counts;
    }
    
    
    std::vector<std::vector<double> >* lexicalTable::getProbabilityTable(){
        if (prob==NULL){
            prob = new(std::nothrow) std::vector<std::vector<double> >;
            
            if (prob==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
			
        }
        
        return prob;
    }
    
    
    std::vector<std::vector<double> >* lexicalTable::getLogProbabilityTable(){
        if (logProb==NULL){
            logProb = new(std::nothrow) std::vector<std::vector<double> >;
            
            if (logProb==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        
        return logProb;
    }
    
    
    double lexicalTable::getValue(sequences& seqs, size_t pos){
        double value;
        
		if (max_order>pos){
			value=_get_reduced_order(seqs, pos);
		}
		else{
			value=_getValue(seqs,pos);
		}
		
        
        return value;
    }
    
    
    //!Add a track to an emission
    //!\param trk Pointer to track
    //!\param orderValue order of emission from track
    void lexicalTable::addTrack(track* trk,int orderValue){
        trcks.push_back(trk);
        alphabets.push_back(trk->getAlphaSize());
        order.push_back(orderValue);
        if (orderValue>max_order){
            max_order=orderValue;
        }
        
    }
	
    //!Set the type of counts in the emission 2D table provided by the user
    //!\param temp vector of vectors of doubles
    //!\param emmType Type of value (COUNTS, PROBABILITY, LOG_PROB)
    void lexicalTable::assignTable(std::vector<std::vector<double > >* temp, valueType emmType){
        if (emmType==COUNTS){
            counts=temp;
        }
        else if (emmType == PROBABILITY){
            prob=temp;
        }
        else if (emmType == LOG_PROB){
            logProb=temp;
        }
    }
	
    
    
    
    // TODO:  Test this function and use to get the required order not just the zeroth order
    
    // TODO: Track index returning wrong value;
    // TODO: Check function for abilty to use ambiguous characters
    //! Calculate the emission probability for a reduced-order emission
    //! Reduced order is called when the order distribution defined is greater then the position in the sequence. For example, when we want to find the what came 5bp before but we're only at position 1 of the sequence.
    //! \param seqs pointer to sequences
    //! \param position Position in the sequence
    //! \return double Score of the transition
    double lexicalTable::_get_reduced_order(sequences& seqs, size_t iter){
        int x=0;
        size_t size=trcks.size();
        int total_size=1;
        
        
        for(size_t i=0;i<size;i++){
            size_t index=trcks[i]->getIndex();
            
            int letter=seqs.seqValue(index,iter);
            int x_subtotal=0;
            x_subtotal+=letter;
            for(size_t j=i+1;j<size;j++){
                x_subtotal*=alphabets[j];
            }
            total_size*=alphabets[i];
            x+=x_subtotal;
        }
        
        std::vector<double> temp (total_size,0);
        
        double sum=0;
        for(size_t i=0;i<counts->size();i++){
            for(size_t j=0;j<(*counts)[i].size();j++){
                temp[j]+=(*counts)[i][j];
                //cout << counts[i][j] <<endl;
                sum+=(*counts)[i][j];
            }
        }
        
        
        if (sum>0 && temp[x]>0){
            return log(temp[x]/sum);
        }
        else {
            return -INFINITY;
        }
    }
    
    
    //TODO: Check conversion size_t to int;
    
    double lexicalTable::_getScore(Index& xVal, Index& yVal){
        if (!xVal.isAmbiguous() && !yVal.isAmbiguous()){
            return (*logProb)[yVal[0]][xVal[0]];
        }
        
        std::vector<double> scores;
        for(size_t i=0;i<yVal.size();i++){
            for(size_t j=0;j<xVal.size();j++){
                scores.push_back((*logProb)[yVal[0]][xVal[0]]);
            }
        }
        
        if (unknownScoreType==AVERAGE_SCORE){
            return avgVector(scores);
        }
        else if (unknownScoreType==LOWEST_SCORE){
            return minVector(scores);
        }
        else if (unknownScoreType==HIGHEST_SCORE){
            return maxVector(scores);
        }
        
        return -INFINITY;
    }
    
    
    //Calculate the index for emissions.   Emissions can be compound tracks, so we need to calculate
    //index for those where we may have an emission A1 or ATG120 and the index for the conditional sequence
    double lexicalTable::_getValue(sequences& seqs, size_t iter){
        
        Index xValue;
        Index yValue;
        
        size_t size=trcks.size();
        
        //Determining Index for applicable sequences and each track
        for(size_t i=0;i<size;i++){
            int index=trcks[i]->getIndex();
            int letter=seqs.seqValue(index,iter);
            
            Index x_subtotal;
            
            if (letter<0){
                if (unknownScoreType==DEFINED_SCORE){
                    return unknownDefinedScore;
                }
                else if (unknownScoreType==NO_SCORE){
                    return -INFINITY;
                }
                else{
                    x_subtotal.setAmbiguous(trcks[i]->getAmbiguousSet(letter));
                }
                
            }
            else{
                x_subtotal+=letter;
            }
            
            for(size_t j=i+1;j<size;j++){
                x_subtotal*=alphabets[j];
            }
            
            xValue+=x_subtotal;
            
            Index y_subtotal;
            
            if (order[i]==0){
                continue;
            }
            else{
                
                for(int k=order[i];k>=1;k--){
                    int prev_letter=seqs.seqValue(index,iter-k);
                    Index letter;
                    if (prev_letter<0){
                        if (unknownScoreType==DEFINED_SCORE){
                            return unknownDefinedScore;
                        }
                        else if (unknownScoreType==NO_SCORE){
                            return -INFINITY;
                        }
                        else{
                            letter.setAmbiguous(trcks[i]->getAmbiguousSet(prev_letter));
                            y_subtotal+=letter*POWER[k-1][alphabets[i]-1];
                        }
                    }
                    else{
                        y_subtotal+=prev_letter*POWER[k-1][alphabets[i]-1];
                    }
                }
                
                for(size_t j=i+1;j<size;j++){
                    y_subtotal*=POWER[order[j]][alphabets[j]-1];  //Calculate offset based on order of next track and the alphabet size -1 of the next track
                }
            }
            yValue+=y_subtotal;
        }
        
        return _getScore(xValue,yValue);
    }
    
    
    
    
    
    std::string lexicalTable::stringify(){
        std::string tbl("");
        size_t tracks_size = trcks.size();
        
        if(tracks_size==0){
            std::cerr << "Can't print out table without track and order being set for lexicalTable\n";
            exit(1);
        }
        
        //Output Column Headers
        size_t columns(1);
        std::vector<size_t> alphaSizes;
        //alphaSizes.push_back(0);
        
        for(size_t i = 0;i<trcks.size();i++){
            size_t alphaSz = trcks[i]->getAlphaSize();
            columns*=alphaSz;
            alphaSizes.push_back(alphaSz);
        }
        
        
        reverse(alphaSizes.begin(),alphaSizes.end());
        
        std::string colHeader("@");
        
        for(size_t i = 0;i<columns;i++){
            size_t indexValue = i;
            for(size_t tr=0;tr<trcks.size();tr++){
                
                if (tr>0){
                    colHeader+= "|";
                }
                
                size_t val(0);
                if (tr<trcks.size()-1){
                    val= floor(indexValue/alphaSizes[tr]);
                    indexValue-=val*alphaSizes[tr];
                }
                else{
                    val = indexValue;
                }
                
                colHeader+=trcks[tr]->convertIndexToWord(val, 1);
            }
            colHeader+="\t";
        }
        
        tbl+=colHeader + "\n";
        
        std::vector<std::vector<double> >* temp;
        
        if (logProb!=NULL){
            temp=logProb;
        }
        else if (prob!=NULL){
            temp=prob;
        }
        else if (counts!=NULL){
            temp=counts;
        }
        else{
            std::cerr << "No table is defined\n";
            return "";
        }
        
        //TODO: Fix row header for other orders
        bool rowHeader = (temp->size()>1) ? true : false;
        for(size_t i=0;i<temp->size();i++){
            std::string header("");
            
            if (rowHeader){
                size_t indexValue = i;
                
                for(size_t tr=0;tr<trcks.size();tr++){
                    
                    if (tr>0 && order[tr]>0){
                        header+= "|";
                    }
                    
                    
                    size_t val(0);
                    if (tr<trcks.size()-1){
                        double pwr = POWER[order[tr+1]][trcks[tr+1]->getAlphaSize()-1];
                        val= floor(indexValue/pwr);
                        indexValue-=val*pwr;
                    }
                    else{
                        val = indexValue;
                    }
                    
                    header+=trcks[tr]->convertIndexToWord(val, order[tr]);
                }
                tbl+="@" + header + "\t";
                
            }
            
            for(size_t j=0;j<(*temp)[i].size();j++){
                if (j>0){
                    tbl+="\t";
                }
                tbl+=double_to_string((*temp)[i][j]);
            }
            tbl+="\n";
        }
        return tbl;
    }
	
    
    
    
}