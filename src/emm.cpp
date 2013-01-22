//emm.cpp
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

#include "emm.h"
namespace StochHMM{

    //! Create an emission
    emm::emm(){
        real_number= false;
        complement = false;
        continuous = false;
        
        function = false;
        lexFunc = NULL;
        
        tagFunc=NULL;
    }

        
    //! Destroy an emission    
    emm::~emm(){
        delete lexFunc;
        delete tagFunc;
        function = false;

        lexFunc=NULL;
        tagFunc=NULL;
    }
    
    //!Parse an emission from text
    //!\param txt  String representation of emission
    //!\param trks Tracks used by the model
    //!\param wts Weights used by the model
    //!\param funcs State functions used by the model
    bool emm::parse(std::string& txt,tracks& trks, weights* wts, StateFuncs* funcs){
        if (!_processTags(txt,trks, wts, funcs)){
            return false;
        }
        
        stringList ln;
        ln.splitString(txt,"\n");
        size_t idx;
        if (ln.contains("EMISSION")){
            idx = ln.indexOf("EMISSION");
        }
        else{
            std::cerr << "Missing EMISSION tag from emission. Please check the formatting.   This is what was handed to the emission class:\n " <<  txt << std::endl;
            return false;
        }
        
        stringList line;
        line.splitString(ln[idx], "\t,: ");
        
        size_t typeBegin;
        
        valueType  valtyp;
        if (line.contains("P(X)")){
            typeBegin = line.indexOf("P(X)");
            valtyp=PROBABILITY;
        }
        else if (line.contains("LOG")){
            typeBegin = line.indexOf("LOG");
            valtyp=LOG_PROB;
        }
        else if (line.contains("COUNTS")){
            typeBegin = line.indexOf("COUNTS");
            valtyp=COUNTS ;
        }
        else if (line.contains("REAL_NUMBER")){
            typeBegin = line.indexOf("REAL_NUMBER");
            real_number = true;
            if (line.contains("COMPLEMENT") || line.contains("1-P(X)")) {
                complement=true;
            }
        }
		else if (line.contains("CONTINUOUS")){
			typeBegin = line.indexOf("CONTINUOUS");
			continuous=true;
			if (line.contains("COMPLEMENT") || line.contains("1-P(X)")) {
                complement=true;
            }
		}
        else if (line.contains("FUNCTION")){
            typeBegin = line.indexOf("FUNCTION");
            function=true;
        }
        else {
            std::string info = "Couldn't parse Value type in the Emission: " + txt  + " Please check the formatting.   The allowed types are: P(X), LOG, COUNTS, or REAL_NUMBER. \n";
            std::cerr << info << std::endl;
            
            //errorInfo(sCantParseLine, info.c_str());
        }
        
        
        //remaining tracks and Orders then set Track
        std::vector<track*> tempTracks;
        for(size_t i=1;i<typeBegin;i++){
            track* tk = trks.getTrack(line[i]);
            if (tk==NULL){
                std::cerr << "Lexical Transition tried to add a track named: " << line[i] << " . However, there isn't a matching track in the model.  Please check to model formatting.\n";
                return false;
            }
            else{
                tempTracks.push_back(tk);
            }
        }
        
        
        if (real_number){
            
            if (tempTracks.size()>1){
                std::cerr << "Multiple tracks listed under Real Track Emission Definition\n";
                return false;
            }
            
            realTrack = tempTracks[0];
            return true;
        }
		else if (continuous){
			
			if (tempTracks.size()>1){
                std::cerr << "Multiple tracks listed under CONTINUOUS Track Emission Definition\n";
                return false;
            }
            
            realTrack = tempTracks[0];
			
			idx = ln.indexOf("PDF");
			line.splitString(ln[idx],"\t:, ");
			
			size_t function_idx = line.indexOf("PDF") + 1;
			pdfName = line[function_idx];
			
			size_t parameter_idx = line.indexOf("PARAMETERS");
			dist_parameters = new(std::nothrow) std::vector<double>;
			
			for(size_t i = parameter_idx+1 ; i< line.size() ; i++){
				double value;
				stringToDouble(line[i], value);
				dist_parameters->push_back(value);
			}
			
			
            pdf = funcs->getPDFFunction(pdfName);
            
            return true;
		}
        else if (function){
            std::string& functionName = line[typeBegin+1];
            lexFunc = new(std::nothrow) emissionFuncParam(functionName,funcs,tempTracks[0]);
            
            if (lexFunc==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            return true;
        }
        else{
            
            
            if (ln.contains("ORDER")){
                idx=ln.indexOf("ORDER");
            }
            else{
                std::cerr << "Couldn't find ORDER in non-Real_Number emission.  Please check the formatting" << std::endl;
                return false;
                //errorInfo(sCantParseLine, "Couldn't find ORDER in non-Real_Number emission.  Please check the formatting\n");
            }
            
            std::vector<int> tempOrder;
            line.splitString(ln[idx],"\t:,");
            
            size_t orderIdx = line.indexOf("ORDER");
            orderIdx++;
            
            size_t ambIdx; 
            bool containsAmbig = line.contains("AMBIGUOUS");
            
            if (containsAmbig){ambIdx=line.indexOf("AMBIGUOUS");}
            else{ ambIdx= line.size();}
            
            for(size_t i=orderIdx;i<ambIdx;i++){
                
                int tempValue;
                if (!stringToInt(line[i], tempValue)){
                    std::cerr << "Emission Order not numeric" << std::endl;
                    return false;
                }
                
                if (tempValue>32){
                    std::cerr << "Emission order is greater than 32.  Must be 32 or less" << std::endl;
                    return false;
                }
                
                tempOrder.push_back(tempValue);
            }
            
            if (tempOrder.size() == tempTracks.size()){
                for(size_t i=0;i<tempOrder.size();i++){
                    scores.addTrack(tempTracks[i], tempOrder[i]);
                }
            }
            else{
                std::cerr << "Different number of tracks and orders parsed in Emission: " << txt << " Check the formatting of the Emission" << std::endl;
                return false;
            }
            
            
            //Parse Ambiguous Tag Info
            if (containsAmbig){
                ambIdx++;
                if (line.size()<=ambIdx){
                    std::cerr << "No scoring type after AMBIGUOUS label\nAssuming AVG\n";
                    scores.setUnkScoreType(AVERAGE_SCORE);
                }
                else if (line[ambIdx].compare("AVG")==0){scores.setUnkScoreType(AVERAGE_SCORE);}
                else if (line[ambIdx].compare("MAX")==0){scores.setUnkScoreType(HIGHEST_SCORE);}
                else if (line[ambIdx].compare("MIN")==0){scores.setUnkScoreType(LOWEST_SCORE);}
                else if (line[ambIdx].compare("P(X)")==0)
                {
                    scores.setUnkScoreType(DEFINED_SCORE);
                    
                    ambIdx++;
                    if (ambIdx>=line.size()){
                        std::cerr << "Missing Ambiguous Value" << std::endl;
                        
                        return false;
                    }
                                    
                    double tempValue;
                    if (!stringToDouble(line[ambIdx], tempValue)){
                        std::cerr << "Ambiguous Value couldn't be parsed: "<< line[ambIdx] << std::endl;
                        return false;
                    }
                    
                    scores.setUnkScore(log(tempValue));
                }
                else if (line[ambIdx].compare("LOG")==0){
                    scores.setUnkScoreType(DEFINED_SCORE);
                    ambIdx++;
                    
                    if (ambIdx>=line.size()){
                        std::cerr << "Missing Ambiguous Value" << std::endl;
                        
                        return false;
                    }
                    
                    double tempValue;
                    if (!stringToDouble(line[ambIdx], tempValue)){
                        std::cerr << "Ambiguous Value couldn't be parsed: "<< line[ambIdx] << std::endl;
                        return false;
                    }
                    
                    scores.setUnkScore(tempValue);
                }
            }
        
            //Get Tables
            int expectedColumns = 1;
            int expectedRows = 1;
            for(size_t i = 0; i<scores.getNumberOfAlphabets(); i++){
                expectedColumns*=scores.getAlphaSize(i);
                expectedRows*=POWER[scores.getOrder(i)][scores.getAlphaSize(i)-1];
            }
            
            std::vector<std::vector<double> >* log_prob = scores.getLogProbabilityTable();
            std::vector<std::vector<double> >* prob = scores.getProbabilityTable();
            std::vector<std::vector<double> >* counts = scores.getCountsTable();
            
            
            for (size_t iter = 2; iter< ln.size();iter++){
                
                
                //If it's the first line check for a '#' indicating that the column header is present
                if (iter==2 && ln[iter][0]=='@'){
                    continue;
                }
                
                line.splitString(ln[iter],"\t ");
                
                //Check for Row header
                if (line[0][0]=='@'){
                    line.pop_ith(0);
                }
                
                
                std::vector<double> temp = line.toVecDouble();
                if (temp.size() != expectedColumns){
                    std::string info = "The following line couldn't be parsed into the required number of columns.   Expected Columns: " + int_to_string(expectedColumns) + "\n The line appears as: "  + ln[iter] ;
                    
                    std::cerr << info << std::endl;
                    return false;
                    //errorInfo(sCantParseLine, info.c_str());
                }
                else{
                    if (valtyp == PROBABILITY){
                        prob->push_back(temp);
                        logVector(temp);
                        log_prob->push_back(temp);
                    }
                    else if (valtyp == LOG_PROB){
                        log_prob->push_back(temp);
                        expVector(temp);
                        prob->push_back(temp);
                    }
                    else if (valtyp == COUNTS){
                        counts->push_back(temp);
                        probVector(temp);
                        prob->push_back(temp);
                        logVector(temp);
                        log_prob->push_back(temp);
                    }
                }
            }
            
            if (log_prob->size() != expectedRows){
                std::cerr << " The Emission table doesn't contain enough rows.  Expected Rows: " << expectedRows << " \n Please check the Emission Table and formatting for " <<  txt << std::endl;
                return false;
            }
			
			scores.initialize_emission_table();
			
        }
        
        return true;
    }
    
    bool emm::_processTags(std::string& txt, tracks& trks,weights* wts, StateFuncs* funcs){
        stringList lst = extractTag(txt);
                
        if (lst.size() == 0){
            return true;
        }
        
        tagFunc=new(std::nothrow) emissionFuncParam();
        
        
        if (tagFunc==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        
        if (!tagFunc->parse(lst, trks, wts,funcs)){
            std::cerr << "Couldn't parse Emission Tag: " << lst.stringify() << std::endl;
            return false;
        }
        std::string trackName = tagFunc->getTrackName();
        
        
        //Don't need if we check in the parsing of emissionFuncParam....
        if (!trks.isTrackDefined(trackName)){
            std::cerr << "No Track defined with name:\t" << trackName << "\nEmission Tag:\n" << txt << std::endl;
            return false;
        }
        
        
        return true;
    }


//TODO: Need to double check everything in this function
//TODO: It's too important to have a bug in it.
    
    //!Calculate the emission value given a position in the sequence
    //!If emission is a real number it will return the value from the real number track
    //!If emission is a sequence then it will get the value and return it
    //!\param seqs  Sequences to use
    //!\param iter Position within the sequences
    //!\return double log(prob) value of emission
    double emm::get_emission(sequences& seqs,size_t pos){
        double final_emission;
        
        if (real_number){
            
            final_emission=seqs.realValue(realTrack->getIndex(),pos);
            if (complement){
                final_emission=log(1-exp(final_emission));
            }
        }
        else if (function){
            final_emission=lexFunc->evaluate(seqs, pos);
        }
		else if (continuous){
			
			final_emission = (*pdf)(seqs.realValue(realTrack->getIndex(),pos),*dist_parameters);
			
			if (complement){
				final_emission=log(1-exp(final_emission));
			}
		}
        else{
            final_emission=scores.getValue(seqs, pos);
        }
        
        
        if (tagFunc!=NULL){
            final_emission+=tagFunc->evaluate(seqs, pos);
        }
        
        return final_emission;
    }


//    // TODO:  Test this function and use to get the required order not just the zeroth order
//
//    // TODO: Track index returning wrong value;
//    //! Calculate the emission probability for a reduced-order emission
//    //! Reduced order is called when the order distribution defined is greater then the position in the sequence. For example, when we want to find the what came 5bp before but we're only at position 1 of the sequence.
//    //! \param seqs pointer to sequences
//    //! \param position Position in the sequence
//    //! \return double Score of the transition
//    double emm::_get_reduced_order(sequences& seqs, int iter){
//        int x=0;
//        size_t size=trcks.size();
//        int total_size=1;
//        
//        
//        for(size_t i=0;i<size;i++){
//            int index=trcks[i]->getIndex();
//            
//            int letter=seqs.seqValue(index,iter);
//            int x_subtotal=0;
//            x_subtotal+=letter;
//            for(size_t j=i+1;j<size;j++){
//                x_subtotal*=alphabets[j];
//            }
//            total_size*=alphabets[i];
//            x+=x_subtotal;
//        }
//        
//        std::vector<double> temp (total_size,0);
//        
//        double sum=0;
//        for(size_t i=0;i<counts->size();i++){
//            for(size_t j=0;j<(*counts)[i].size();j++){
//                temp[j]+=(*counts)[i][j];
//                //cout << counts[i][j] <<endl;
//                sum+=(*counts)[i][j];
//            }
//        }
//        
//        
//        if (sum>0 && temp[x]>0){
//            return log(temp[x]/sum);
//        }
//        else {
//            return -INFINITY;
//        }
//    }
//
//
//    double emm::_getScore(Index& xVal, Index& yVal){
//        if (!xVal.isAmbiguous() && !yVal.isAmbiguous()){
//            return (*log_emm)[yVal[0]][xVal[0]];
//        }
//        
//        std::vector<double> scores;
//        for(size_t i=0;i<yVal.size();i++){
//            for(size_t j=0;j<xVal.size();j++){
//                scores.push_back((*log_emm)[yVal[0]][xVal[0]]);
//            }
//        }
//        
//        if (unknownScoreType==AVERAGE_SCORE){
//            return avgVector(scores);
//        }
//        else if (unknownScoreType==LOWEST_SCORE){
//            return minVector(scores);
//        }
//        else if (unknownScoreType==HIGHEST_SCORE){
//            return maxVector(scores);
//        }
//        
//        return -INFINITY;
//    }
//
//
//    //Calculate the index for emissions.   Emissions can be compound tracks, so we need to calculate
//    //index for those where we may have an emission A1 or ATG120 and the index for the conditional sequence
//    double emm::_getValue(sequences& seqs, int iter){
//        
//        Index xValue;
//        Index yValue;
//        
//        size_t size=trcks.size();
//        
//        //Determining Index for applicable sequences and each track
//        for(size_t i=0;i<size;i++){
//            int index=trcks[i]->getIndex();
//            int letter=seqs.seqValue(index,iter);
//            
//            Index x_subtotal;
//
//            if (letter<0){
//                if (unknownScoreType==DEFINED_SCORE){
//                    return unknownDefinedScore;
//                }
//                else if (unknownScoreType==NO_SCORE){
//                    
//                }
//                else{
//                    x_subtotal.setAmbiguous(trcks[i]->getAmbiguousSet(letter));
//                }
//                
//            }
//            else{
//                x_subtotal+=letter;
//            }
//            
//            for(size_t j=i+1;j<size;j++){
//                x_subtotal*=alphabets[j];
//            }
//            
//            xValue+=x_subtotal;
//            
//            Index y_subtotal;
//            
//            if (order[i]==0){
//                continue;
//            }
//            else{
//                
//                for(int k=order[i];k>=1;k--){
//                    int prev_letter=seqs.seqValue(index,iter-k);
//                    Index letter;
//                    if (prev_letter<0){
//                        if (unknownScoreType==DEFINED_SCORE){
//                            return unknownDefinedScore;
//                        }
//                        else if (unknownScoreType==NO_SCORE){
//                            return -INFINITY;
//                        }
//                        else{
//                            letter.setAmbiguous(trcks[i]->getAmbiguousSet(prev_letter));
//                            y_subtotal+=letter*POWER[k-1][alphabets[i]-1];
//                        }
//                    }
//                    else{
//                        y_subtotal+=prev_letter*POWER[k-1][alphabets[i]-1];
//                    }
//                }
//                
//                for(size_t j=i+1;j<size;j++){
//                    y_subtotal*=POWER[order[j]][alphabets[j]-1];  //Calculate offset based on order of next track and the alphabet size -1 of the next track
//                }
//            }
//            yValue+=y_subtotal;
//        }
//        
//        return _getScore(xValue,yValue);
//    }
    
//    //!Add a track to an emission
//    //!\param trk Pointer to track
//    //!\param orderValue order of emission from track
//    void emm::addTrack(track* trk,int orderValue){
//        trcks.push_back(trk);
//        alphabets.push_back(trk->getAlphaSize());
//        order.push_back(orderValue);
//        if (orderValue>max_order){
//            max_order=orderValue;
//        }
//        
//    }
    
    
    //!Is the emission from a real Number track
    //!\return true if track is real number track
    //!\return false if track is alphanumerical track
    bool emm::isReal(){
        if (real_number && realTrack->getAlphaType()==REAL){
            return true;
        }
        else{
            return false;
        }
    }

//    //!Set the type of counts in the emission 2D table provided by the user
//    //!\param temp vector of vectors of doubles 
//    //!\param emmType Type of value (COUNTS, PROBABILITY, LOG_PROB)
//    void emm::assignMatrix(std::vector<std::vector<double > >* temp, valueType emmType){
//        if (emmType==COUNTS){
//            counts=temp;
//        }
//        else if (emmType == PROBABILITY){
//            emmi=temp;
//        }
//        else if (emmType == LOG_PROB){
//            log_emm=temp;
//        }
//    }
        
    //Get the string representation of the emission
    //\return std::string 
    std::string emm::stringify(){
        std::string emissionString("EMISSION:\t");
        
        if (real_number){
            
            emissionString+=realTrack->getName();
            emissionString+=":\t";
            
            emissionString+="REAL_NUMBER";
            if (complement){
                emissionString+=":\tCOMPLEMENT\t";
            }
            else{
                emissionString+="\t";
            }
            
            
            if (tagFunc){
                emissionString+=tagFunc->stringify();
            }
            emissionString+="\n";
        }
		else if (continuous){
			emissionString+=realTrack->getName();
			emissionString+="\tCONTINUOUS";
			
			if (tagFunc){
                emissionString+=tagFunc->stringify();
            }
            emissionString+="\n\t";
			
			emissionString+="PDF:\t";
			emissionString+=pdfName + "\tPARAMETERS:\t";
			emissionString+=join(*dist_parameters, ',');
			emissionString+= "\n";
			
		}
        else if (function){
            emissionString+=lexFunc->getTrack()->getName();
            emissionString+="\tFUNCTION:\t";
            emissionString+=lexFunc->getName();
            emissionString+="\t";
            
            if (tagFunc){
                emissionString += tagFunc->stringify();
            }
            
            emissionString+="\n";
        }
        else{
            
            for(size_t i=0;i<scores.trackSize();i++){
                if (i>0){
                    emissionString+=",";
                }
                emissionString+=scores.getTrack(i)->getName();
            }
            
            emissionString+=":\t";
            
            emissionString+="LOG";
            
            if (tagFunc){
                
                emissionString += "\t";
                
                emissionString += tagFunc->stringify();
            }
            
            emissionString+="\n\tORDER:\t";
            
            for(size_t i=0;i<scores.trackSize();i++){
                if (i>0){
                    emissionString+=",";
                }
                emissionString+=int_to_string(scores.getOrder(i));
            }
            
            unknownCharScoringType  ambTemp = scores.getAmbScoringType();
            
            if (ambTemp!=NO_SCORE){
                emissionString+="\tAMBIGUOUS:\t";
                emissionString+=(ambTemp==HIGHEST_SCORE)? "MAX":
                                (ambTemp==LOWEST_SCORE)? "MIN":
                                (ambTemp==AVERAGE_SCORE)? "AVG": "LOG:" + double_to_string(scores.getAmbDefinedScore());
            }
            emissionString+="\n";

            emissionString+=scores.stringify();
        }
        emissionString+="\n";
        
        return emissionString;
    }

}
