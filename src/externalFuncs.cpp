//
//  externalFuncs.cpp
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

#include "externalFuncs.h"
namespace StochHMM{
    
    //TODO: Create non-model import way of creating and defining externalFuncs class
    
    //!Create transFuncParam from a stringList parsed from model for given track and applying a weight as defined in the model
    //! \param lst  stringList from parsing the line in the function
    //! \param trcks tracks defined in the model
    //! \param wts pointer to weights as defined in the model
    //! \param funcs pointer to the state functions defined by user and referenced in the model
    
    transitionFuncParam::transitionFuncParam(){
        transFunc=NULL;
    }
    
//    transitionFuncParam::transitionFuncParam(stringList& lst, tracks& trcks, weights* wts, StateFuncs* funcs){
//        transFunc=NULL;
//        parse(lst,trcks, wts,funcs);
//    }

    
    //!Parse the stringList from model import 
    bool transitionFuncParam::parse(stringList& lst, tracks& trcks, weights* wts, StateFuncs* funcs){
        
        size_t idx;
        
        //FUNCTION NAME (REQUIRED)
        if (lst.contains("FUNCTION")){
            idx=lst.indexOf("FUNCTION");
            idx++;
            transFuncName = lst[idx];
            if (funcs!=NULL){
                transFunc = funcs->getTransitionFunction(transFuncName);
            }
        }
        else{
            std::cerr << "Tag was parsed but contains no FUNCTION: . Please check the formatting of the tags\n" << std::endl;
            return false;
            
            //errorInfo(sCantParseLine, "Tag was parsed but contains no FUNCTION: . Please check the formatting of the tags\n");
        }
        
        //Implement Which track to pass to function
        
        if (lst.contains("TRACK")){
            idx=lst.indexOf("TRACK");
            idx++;
            trackName=lst[idx];
            transFuncTrack = trcks.getTrack(trackName);
        }
        else{
            std::cerr << "Tag was parsed but contains no TRACK: . Please check the formatting of the tags\n" << std::endl;
            return false;
            //errorInfo(sCantParseLine, "Tag was parsed but contains no TRACK: . Please check the formatting of the tags\n");
        }
            
        
        
        if (lst.contains("SCALE")){
            idx=lst.indexOf("SCALE");
            idx++;
            if (isNumeric(lst[idx])){
                transFuncScaling = new(std::nothrow) weight;
                
                if (transFuncScaling==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                double tempValue;
                if (!stringToDouble(lst[idx], tempValue)){
                    std::cerr << "Ambiguous Value couldn't be parsed: "<< lst[idx] << std::endl;
                    return false;
                }
                
                transFuncScaling->setAbsolute(tempValue);
            }
            else{
                std::string weightName=lst[idx];
                if (wts->count(weightName)){
                    transFuncScaling = (*wts)[weightName];
                }
            }
            
        }
        
        
        //Process Traceback Commands
        
        //Parse the TO Traceback Labels in the function tag
        const std::string tbLabels[]= {"TO_LABEL","TB->LABEL","TO_GFF","TB->GFF","TO_STATE","TB->STATE","TO_START","TB->START","DIFF_STATE"};
        const tracebackIdentifier typs[]= {STATE_LABEL,STATE_LABEL,STATE_GFF,STATE_GFF,STATE_NAME,STATE_NAME,START_INIT,START_INIT,DIFF_STATE};
        
        for(int i=0;i<9;i++){
            if (lst.contains(tbLabels[i])){
                
                transFuncTraceback=true;
                
                idx=lst.indexOf(tbLabels[i]);
                transFuncTracebackIdentifier=typs[i];
                
                if (typs[i]!=START_INIT || typs[i]!=DIFF_STATE){
                    transFuncTracebackString=lst[idx+1];
                }
            }
        }
        
        
        // Parse the Combine tags
        if (transFuncTraceback){
            //Process Traceback Combining Commands
            std::string combineLabels[] = {"COMBINE_LABEL","COMBINE_GFF","COMBINE_STATE", "NO_COMBINE"};
            const combineIdentifier comtyps[] = {STATELABEL,STATEGFF,STATENAME,FULL};
            bool combineTypeProvided=false;
            for(int i=0;i<4;i++){
                if (lst.contains(combineLabels[i])){
                    idx=lst.indexOf(combineLabels[i]);
                    transFuncCombineIdentifier=comtyps[i];
                    transFuncCombineString=lst[idx+1];
                    combineTypeProvided=true;
                }
            }
            
            if (!combineTypeProvided){
                std::cerr << "Transition Function tag with a traceback was called, but no CombineType was provided.  If no traceback is needed then please remove traceback command. \n" << std::endl;
                return false;
                
                //errorInfo(sTagIncorrect, "Transition Function tag with a traceback was called, but no CombineType was provided.  If no traceback is needed then please remove traceback command. \n");
            }
        }
        
        
        return true;
    }
    
    //!Print the externalFuncts to stdout
    void transitionFuncParam::print(){
        std::cout<<stringify()<<std::endl;
    }
    
    //!Get string representation of the externalFuncs
    //! \return std::string Representation of externalFuncs
    std::string transitionFuncParam::stringify(){
        std::string exFuncString;
        exFuncString+="[ FUNCTION:\t" + transFuncName + "\t";
        exFuncString+="TRACK:\t" + trackName + "\t";
        if (transFuncTraceback){
            exFuncString+= (transFuncTracebackIdentifier == STATE_NAME)  ? "TO_STATE:\t" + transFuncTracebackString + "\t":
                           (transFuncTracebackIdentifier == STATE_LABEL) ? "TO_LABEL:\t" + transFuncTracebackString + "\t": 
                           (transFuncTracebackIdentifier == STATE_GFF)   ? "TO_GFF:\t" + transFuncTracebackString   + "\t":
                           (transFuncTracebackIdentifier == START_INIT)  ? "TO_START:\t": "DIFF_STATE:\t" ;
            
            exFuncString+=  (transFuncCombineIdentifier == STATENAME) ? "COMBINE_STATE:\t" + transFuncCombineString + "\t":
                            (transFuncCombineIdentifier == STATELABEL) ? "COMBINE_LABEL:\t" + transFuncCombineString + "\t":
                            (transFuncCombineIdentifier == STATEGFF) ? "COMBINE_GFF:\t" + transFuncCombineString :"NO_COMBINE:\t";
        }
        
        if (transFuncScaling!=NULL){
            if (transFuncScaling->isAbsolute()){
                exFuncString+="SCALE:\t" + double_to_string(transFuncScaling->getAbsolute());
            }
            else{
                exFuncString+="SCALE:\t" + transFuncScaling->getName();
            }
        }
        
        exFuncString += " ]";
        
        return exFuncString;
    }
    
    
    //! Add function to externFunc
    //! \param funcName  std::string name given to function (as referenced in model)
    //! \param function Pointer to pt2StateFunc*, function to use
    //! \param scaling Pointer to weight, how to weight the functions results before applying to emission/transition
    void transitionFuncParam::setTransFunc(std::string& funcName, transitionFunc* function, weight* scaling){
        transFunc=function;
        transFuncName=funcName;
        transFuncScaling=scaling;
    }
    
    //!Set all the traceback and combine types for the traceback before calling the function
    //! \param tbIdent tracebackIdentifier 
    //! \param tbString Traceback string to use (State/GFF/Label)
    //! \param cbIdent  combineIdentifier
    //! \param cbString Combine string to use when editing traceback (State/GFF/Label)
    void transitionFuncParam::setTransTB(tracebackIdentifier tbIdent,std::string& tbString, combineIdentifier cbIdent, std::string& cbString){
        transFuncTracebackIdentifier=tbIdent;
        transFuncTracebackString=tbString;
        transFuncCombineIdentifier=cbIdent;
        transFuncCombineString=cbString;
        transFuncTraceback=true;
        return;
    }
    
    
    double transitionFuncParam::evaluate(const std::string* fullSequence, size_t pos, const std::string* partialSequence,size_t length){
        return (*transFunc)(fullSequence, pos, partialSequence,length);
    }
    
    
    
    emissionFuncParam::emissionFuncParam(std::string& functionName, StateFuncs* funcs,track* trk){
        emissionFunction=NULL;
        
        emissionFunction=funcs->getEmissionFunction(functionName);
        
        emissionFuncTrack=trk;
        trackName=trk->getName();
        trackNumber=trk->getIndex();
        emissionFuncName = functionName;
        return;
    }
    
    
    //!Create emissionFuncParam from a stringList parsed from model for given track and applying a weight as defined in the model
    //! \param lst  stringList from parsing the line in the function
    //! \param trcks tracks defined in the model
    //! \param wts pointer to weights as defined in the model
    //! \param funcs pointer to the state functions defined by user and referenced in the model
    
    emissionFuncParam::emissionFuncParam(){
        emissionFunction=NULL;
        //parse(lst,trcks, wts,funcs);
    }
    
    
    //!Parse the stringList from model import 
    bool emissionFuncParam::parse(stringList& lst, tracks& trcks, weights* wts, StateFuncs* funcs){
        
        size_t idx;
        
        //FUNCTION NAME (REQUIRED)
        if (lst.contains("FUNCTION")){
            idx=lst.indexOf("FUNCTION");
            idx++;
            emissionFuncName = lst[idx];
            if (funcs!=NULL){
                emissionFunction = funcs->getEmissionFunction(emissionFuncName);
            }
        }
        else{
            std::cerr << "Tag was parsed but contains no FUNCTION: . Please check the formatting of the tags\n" << std::endl;
            return false;
            
            //errorInfo(sCantParseLine, "Tag was parsed but contains no FUNCTION: . Please check the formatting of the tags\n");
        }
        
        //Implement Which track to pass to function
        
        if (lst.contains("TRACK")){
            idx=lst.indexOf("TRACK");
            idx++;
            trackName=lst[idx];
            emissionFuncTrack = trcks.getTrack(trackName);
            trackNumber = emissionFuncTrack->getIndex();
        }
        else{
            std::cerr << "Tag was parsed but contains no TRACK: . Please check the formatting of the tags\n" << std::endl;
            return false;
            
            //errorInfo(sCantParseLine, "Tag was parsed but contains no TRACK: . Please check the formatting of the tags\n");
        }
        
        
        
        if (lst.contains("SCALE")){
            idx=lst.indexOf("SCALE");
            idx++;
            if (isNumeric(lst[idx])){
                emissionFuncScaling = new(std::nothrow) weight;
                
                if (emissionFuncScaling==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                double tempValue;
                if (!stringToDouble(lst[idx], tempValue)){
                    std::cerr << "SCALE value could not be converted to numerical value: "<< lst[idx] << std::endl;
                    return false;
                }
                
                emissionFuncScaling->setAbsolute(tempValue);
            }
            else{
                std::string weightName=lst[idx];
                if (wts->count(weightName)){
                    emissionFuncScaling = (*wts)[weightName];
                }
            }
            
        }
        
        return true;
    }
    
    //!Print the externalFuncts to stdout
    void emissionFuncParam::print(){
        std::cout<<stringify()<<std::endl;
    }
    
    //!Get string representation of the externalFuncs
    //! \return std::string Representation of externalFuncs
    std::string emissionFuncParam::stringify(){
        std::string exFuncString;
        exFuncString+="[ FUNCTION:\t" + emissionFuncName + "\t";
        exFuncString+="TRACK:\t" + trackName + "\t";

        if (emissionFuncScaling!=NULL){
            if (emissionFuncScaling->isAbsolute()){
                exFuncString+="SCALE:\t" + double_to_string(emissionFuncScaling->getAbsolute());
            }
            else{
                exFuncString+="SCALE:\t" + emissionFuncScaling->getName();
            }
        }
        
        exFuncString += " ]";
        
        return exFuncString;
    }
    
    
    double emissionFuncParam::evaluate(const std::string* fullSequence, size_t pos){
        double val = (*emissionFunction)(fullSequence, pos);
        
        if (emissionFuncScaling!=NULL){
            val = emissionFuncScaling->getWeightedScore(val);
        }
        
        return val;
    }
    
    double emissionFuncParam::evaluate(sequences& seqs , size_t pos){
        double val = (*emissionFunction)(seqs.getUndigitized(trackNumber),pos);
        
        if (emissionFuncScaling!=NULL){
            val = emissionFuncScaling->getWeightedScore(val);
        }
        
        return val;
    }
	
	double emissionFuncParam::evaluate(sequence& seq , size_t pos){
        double val = (*emissionFunction)(seq.getUndigitized(),pos);
        
        if (emissionFuncScaling!=NULL){
            val = emissionFuncScaling->getWeightedScore(val);
        }
        
        return val;
    }
    
}