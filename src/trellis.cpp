//
//  newTrellis.cpp


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

#include "trellis.h"
namespace StochHMM{

    //  Input:  HMM* model, sequences* sequences, int type of trellis
    //  Types of Trellis:  Simple, Stochastic, Nth
    //   
    //  Stores pointer to trellis in trell variable
    // 
    //----------------------------------------------------------------------------//
    //! Create a trellis
    //! \param mdl Pointer to model
    //! \param seqs Pointer to sequences
    //! \param type  Type of trellis to setup

    trellis::trellis(model* mdl, sequences* seqs, trellisType type):path(mdl),stochPath(mdl){
        hmm=mdl;
        seq=seqs;
        
        stoch=false;
        simple=false;
        
        forwardCompleted=false;
        viterbiCompleted=false;
        backwardCompleted=false;
        posteriorCompleted=false;
        
        trell=NULL;
        
        switch (type){
            case 0:
                _initSimple();  //Init trell with pointer to simple trellis
                break;
            case 1:
                _initStochastic();  //Init trell to stochastic trellis
                break;
            default:
                break;
        }
        
        return;
    }
    
    

    //!Setup the trellis as a simple trellis and create the trellis
    bool trellis::_initSimple(){
        if (simple || stoch){
            simple=true;
            stoch=false;
            delete trell;
        }
        
        trell=new(std::nothrow) simpleTrellis(hmm,seq);
        if (trell==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        return true;
    }


    //!Setup the trellis as a stochastic trellis and create the trellis
    bool trellis::_initStochastic(){
        if (simple || stoch){
            simple=false;
            stoch=true;
            delete trell;
        }
        
        trell=new(std::nothrow) stochTrellis(hmm,seq);
        if (trell==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        return true;
    }


    //!Perform viterbi algorithm
    void trellis::viterbi(){
        
        
        (*this)->initViterbi();
        
        while ((*this)->moveToNext()){
            (*this)->calcViterbi();
    #ifdef DEBUG_VERBOSE
            std::cout <<"Current State: " <<  (*this)->currentState <<std::endl;
            std::cout <<"Sequence position: " <<  (*this)->sequencePosition <<std::endl;
            std::cout <<"Viterbi Score: " << (*this)->getViterbi()/log(2) <<std::endl << std::endl;
    #endif
        
        } 
        
        (*this)->calcEndViterbi();
        //(*this)->finalize();
        (*this)->resetCounters();
        
        viterbiCompleted=true;
    }
    
    
    
    //!Perform Nth viterbi algorithm
    void trellis::nthViterbi(size_t n){
        
        (*this)->initNthViterbi(n);
        
        while ((*this)->moveToNext()){
            (*this)->calcNthViterbi(n);
#ifdef DEBUG_VERBOSE
            std::cout <<"Current State: " <<  (*this)->currentState <<std::endl;
            std::cout <<"Sequence position: " <<  (*this)->sequencePosition <<std::endl;
            std::cout <<"Viterbi Score: " << (*this)->getViterbi()/log(2) <<std::endl << std::endl;
#endif
            
        }
        
        (*this)->calcNthEndViterbi(n);
        //(*this)->finalize();
        (*this)->resetCounters();
        
        viterbiCompleted=true;
        nthViterbiCompleted = true;
    }
    

    
    
    //!Perform forward algorithm
    void trellis::forward(){
        
        if (!hmm->isBasic()){
            this->forwardViterbi();
            return;
        }
        
        (*this)->initForward();
        
        while ((*this)->moveToNext()){
            (*this)->calcForward();


    #ifdef DEBUG_VERBOSE
            std::cout <<"Current State: " <<  (*this)->currentState <<std::endl;
            std::cout <<"Sequence position: " <<  (*this)->sequencePosition <<std::endl;
            std::cout <<"Forward Score: " << exp((*this)->getForward()) <<std::endl << std::endl;
    #endif

        
        } 

        (*this)->calcEndForward();

        


    #ifdef DEBUG_VERBOSE
        std::cout << "Sequence Probability: " << exp((*this)->probabilityOfSequence)  << std::endl;
    #endif
        (*this)->resetCounters();
        forwardCompleted=true;
        
        return;
    }

    //!Perform the backward algorithm
    void trellis::backward(){
        
        if (!hmm->isBasic()){
            this->viterbi();
        }
        
        (*this)->initBackward();
        
        
        while ((*this)->moveToPrevious()){
            
            (*this)->calcBackward();
            
            
#ifdef DEBUG_VERBOSE
            std::cout <<"Current State: " <<  (*this)->currentState <<std::endl;
            std::cout <<"Previous State: " << (*this)->previousState <<std::endl;
            std::cout <<"Sequence position: " <<  (*this)->sequencePosition - 1 <<std::endl;
            std::cout <<"Backward Score: " << exp((*this)->getBackward((*this)->sequencePosition-1,(*this)->currentState)) <<std::endl << std::endl;
#endif
            
            
        }
        // TODO:  Implement calcBeginBackward(); May not really need to do this because we already get P(sequence) from forward algorithm
        
        //(*this)->calcBeginBackward();   
        
        (*this)->resetCounters();
        backwardCompleted=true;
    }




    //!Perform the forward and Viterbi algorithm
    void trellis::forwardViterbi(){
        
        (*this)->initForwardViterbi();
        
        while ((*this)->moveToNext()){
            (*this)->calcForwardViterbi();
        }
        
        (*this)->calcEndForwardViterbi();
        
        (*this)->resetCounters();
        
        forwardCompleted=true;
        viterbiCompleted=true;
        
        return;
    }


    //!Perform both forward and Backward algorithms
    void trellis::posterior(){
        
        if (!hmm->isBasic()){
            this->forwardViterbi();
        }
        else{
            //Perform Forward
            forward();
        }
        
        backward();
        
        (*this)->calcPosterior();
        
        //(*this)->calcBeginBackward();
        posteriorCompleted=true;
        forwardCompleted=true;
        backwardCompleted=true;
        
    }

    
    //!Perform both forward, viterbi, and backward algorithm
    void trellis::decodeAll(){
        //Perform Forward & Viterbi
        forwardViterbi();
        backward();
        
        (*this)->calcPosterior();
        
        forwardCompleted=true;
        backwardCompleted=true;
        posteriorCompleted=true;
        viterbiCompleted=true;
        
    }

    //Perform a traceback on the trellis
    traceback_path& trellis::traceback(){
        if (path.size()!=0){
            path.clear();
        }
        
        (*this)->traceback(path);
        
        if (path.size() == 0){
            
            std::string info = "Path not valid"  + seq->getHeader() ;
            
            //errorInfo(sInvalidTraceback, info.c_str());
            std::cerr << info << std::endl;
            exit(1);
        }
        
        return path;
    }
    
    
    //Perform a nth viterbi traceback on the trellis
    traceback_path& trellis::traceback(size_t n){
        if (path.size()!=0){
            path.clear();
        }
        
        (*this)->traceback(path, n);
        
        if (path.size() == 0){
            
            std::string info = "Path not valid"  + seq->getHeader() ;
            
            //errorInfo(sInvalidTraceback, info.c_str());
            std::cerr << info << std::endl;
            exit(1);
        }
        
        return path;
    }



    //!Perform multiple stochastic Tracebacks on the trellis
    //!\param numberOfTracebacks  Number of tracebacks to perform
    //!\param type Type of decoding to use (enum decodingType)
    multiTraceback& trellis::stochasticTraceback(int numberOfTracebacks,decodingType type){
        if (paths.size()!=0){
            paths.clear();
        }
        
        int numberOfPaths(0);
        while(numberOfPaths<numberOfTracebacks){
            if (type==VITERBI){
                stochasticViterbiTraceback();
            }
            else{ //Forward Traceback
                stochasticForwardTraceback();
            }
            
            if (path.size()!=0){
                numberOfPaths++;
                paths.assign(stochPath);
                stochPath.clear();
            }
        }
        return paths;
    }


    //!Performs one stochastic traceback using viterbi data stores the traceback_path pointer to pathif path is already define then we delete it
    
    traceback_path& trellis::stochasticViterbiTraceback(){
        if (path.size()!=0){
            path.clear();
        }
        
        (*this)->traceStochViterbi(stochPath);
        
        return stochPath;
    }


    //!Performs one stochastic traceback using forward data stores the traceback_path pointer to pathif path is already define then we delete it
    traceback_path& trellis::stochasticForwardTraceback(){
        if (path.size()!=0){
            path.clear();
        }
        
        (*this)->traceStochForward(stochPath);
        
        return stochPath;
    }


    //!Perform a set number of tracebacks using the decodingType                                                           
    //!\param numberOfTracebacks  Number of tracebacks to perform
    //!\return multiTraceback&
    multiTraceback& trellis::stochasticViterbiTraceback(int numberOfTracebacks){
        if (paths.size()!=0){
            paths.clear();
        }
        
        int numberOfPaths(0);
        
        while(numberOfPaths<numberOfTracebacks){
            stochasticViterbiTraceback();
        
            if (path.size()!=0){
                numberOfPaths++;
                paths.assign(stochPath);
                stochPath.clear();
            }
        }
        return paths;
    }


   
    //!Perform a set number of tracebacks using the decodingType                                                           
    multiTraceback& trellis::stochasticForwardTraceback(int numberOfTracebacks){
        if (paths.size()!=0){
            paths.clear();
        }
        
        int numberOfPaths(0);
        
        while(numberOfPaths<numberOfTracebacks){
            stochasticForwardTraceback();
            
            if (path.size()!=0){
                numberOfPaths++;
                paths.assign(stochPath);
                stochPath.clear();
            }
        }
        return paths;
    }
    
    //!Perform a traceback using the posterior score
    traceback_path& trellis::posteriorTraceback(){
        if (path.size()!=0){
            path.clear();
        }
        
        if (backwardCompleted && forwardCompleted){
            (*this)->tracePosterior(path);
        }
        else{
            std::cerr << "Must complete forward and backward algorithm before tracing posterior";
        }
        
        
        return path;
    }
    
    traceback_path& trellis::stochasticPosterior(){
        if (path.size()!=0){
            path.clear();
        }
        
        if (backwardCompleted && forwardCompleted){
            (*this)->traceStochPosterior(path);
        }
        else{
            std::cerr << "Must complete forward and backward algorithm before tracing posterior";
        }
        return path;
    }

    //!Print the trellis to stdout
    void trellis::print(){
        (*this)->print();
        return;
    }





}