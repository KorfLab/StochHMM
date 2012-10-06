//
//  StochError.h

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
#ifndef StochHMM_StochErr_h
#define StochHMM_StochErr_h
#include <iostream>
#include <string>
#include <deque>
#include <exception>
#include <time.h>
#include <sstream>

namespace StochHMM{
    //! \file stochErr.h
    //! Define the exception handling and StochHMM errors

    //! \def #define error(x)  throw StochHMM::stochErr(x,__FUNCTION__,__FILE__,__LINE__);
    //! Error macro for throwing error
    //!\param x Error message to throw
    #define error(x)  throw StochHMM::stochErr(x,__FUNCTION__,__FILE__,__LINE__);
    
    
    
    //! \def #define errorInfo(x,y) throw StochHMM::stochErr(x,y,__FUNCTION__,__FILE__,__LINE__);
    //! Error macro for throwing error
    //!\param x Error message to throw
    //!\param y Additional error information to pass from function
    #define errorInfo(x,y) throw StochHMM::stochErr(x,y,__FUNCTION__,__FILE__,__LINE__);
    
    
    
    static const char* priority[] = 
    {
        "Debug",
        "Notice",
        "Warnings",
        "Critical",
        "Emergency"
    };
    
    //!Additional error information
    static const char* ErrInfo[] = 
    {
        "sNoError - No Error",
        "sUnknownError - Unknown Error",
        // General Logic Errors
        "sDomainError - Value passed to function was out of range",
        "sInvalidArgument - Invalid argument passed to function",
        "sLengthError - Container was too large",
        "sOutOfRange - Index used was out of range of container",
        // General Runtime Errors
        "sRangeError - Computation produced result that was out of range",
        "sOverflowError - Arithmetic overflow error",
        "sUnderflowError - Arithmetic underflow error",
        // General Exceptions
        "sOutOfMemory - Program ran out of memory",
        "sDynCasting - Dynamic casting error",
        //StochHMM General Errors
        "sFileNotFound - File not Found",
        "sCantOpenFile - Can't open the file",
        "sCantWriteFile - Can't write to file",
        "sUnexpectedEOF - Unexpected End Of File encountered",
        "sFileAlreadyExists - File already exists",
        "sCantParseLine - Can't parse line of information properly",
        "sTagIncorrect - Tag was formatted incorrectly or unable to parse",
        
        //StochHMM Model
        "sCantImportModel - Can't import the model",
        "sIncompletePathThroughModel - Model doesn't properly define state machine"
        
        //StochHMM Weights
        "sInvalidWeightFormatting - Invalid formatting of weights",
        
        //StochHMM Tracks
        "sMissingTrackDefinition - Definition of Track missing from the model file or Not able to parse the Track Definitions correctly",
        "sInvalidTrackFormatting - Invalid formatting of Track Definitions",
        "sCantHandleAmbiguousCharacter - Ambiguous character encountered while ambiguous characters are not allowed"
        "sAmbigousCharacterDoesntMatchTrack - Track name for Ambiguous character definition doesn't match any defined track"
        
        //StochHMM Templates
        "sMissingTemplate - No template section was defined or was unable to parse the template section",
        
        //StochHMM State
        "sMissingStateDefinition - State Definition is missing or unable to parse STATE DEFINITION section",
        "sCantImportState - Can't import the state",
        
        //StochHMM Initial State
        "sMissingInitialState - Initial  (INIT) state is missing or unable to parse the initial state",
        "sCantImportInitialState - Can't import the initial state",
        
        //StochHMM Transition
        
        //StochHMM Emission
        
        //StochHMM User Function Definitions
        
        //StochHMM Sequence
        "sCantDigitizeSequence - Can't digitize the sequence",
        "sInvalidFormatting - Invalid sequence formatting",
        
        //StochHMM Sequences
        "sDifferenSizeSequences - Sequence sizes are different sizes",
        
        //StochHMM External Definitions
        
        //StochHMM SeqTracks
        
        //StochHMM SeqJobs
        
        //StochHMM Trellis
        
        //StochHMM Traceback
        "sNoValidTraceback - No valid traceback through the trellis was possible",
        "sInvalidTraceback - Invalid traceback through the trellis"
        //StochHMM Output
    };
    
    //!\enum enum errStoch { sNoError=0,sUnknownError, sDomainError,sInvalidArgument, sLengthError, sOutOfRange, sRangeError,sOverflowError,sUnderflowError,sOutOfMemory,sDynCasting,sFileNotFound,sCantOpenFile,sCantWriteToFile,sUnexpectedEOF,sFileAlreadyExists,sCantParseLine,sTagIncorrect,sCantImportModel,sIncompletePathThroughModel,sInvalidWeightFormatting,sMissingTrackDefinition,sInvalidTrackFormatting,sCantHandleAmbiguousCharacter,sAmbigousCharacterDoesntMatchTrack,sMissingTemplate,sMissingStateDefinition,sCantImportState,sMissingInitialState,sCantImportInitialState,sCantDigitizeSequence,sInvalidFormatting,sDifferenSizeSequences, sNoValidTraceback,sInvalidTraceback};
    
    enum errStoch 
    { 
        sNoError=0,          // No Err
        sUnknownError,        // Unknown Error
        // General Logic Errors
        sDomainError,        // Value out of domain of function
        sInvalidArgument,    // Invalid argument passed to function
        sLengthError,        // Container size too large
        sOutOfRange,         // Index is out of range of container
        // General Runtime Errors
        sRangeError,         // Internal Computation out of range
        sOverflowError,      // Arithmetic overflow error
        sUnderflowError,     // Arithmetic underflow error
        //General Exceptions
        sOutOfMemory,        // bad_alloc (can't allocate memory)
        sDynCasting,         // Dynamic Casting Error
        
        //StochHMM Err General
        sFileNotFound,
        sCantOpenFile,
        sCantWriteToFile,
        sUnexpectedEOF,
        sFileAlreadyExists,
        sCantParseLine,
        sTagIncorrect,
        
        //StochHMM Model
        sCantImportModel,
        sIncompletePathThroughModel,
        
        //StochHMM Weights
        sInvalidWeightFormatting,
        
        //StochHMM Tracks
        sMissingTrackDefinition,
        sInvalidTrackFormatting,
        sCantHandleAmbiguousCharacter,
        sAmbigousCharacterDoesntMatchTrack,
        
        //StochHMM Templates
        sMissingTemplate,
        
        //StochHMM State
        sMissingStateDefinition,
        sCantImportState,
        
        //StochHMM Initial State
        sMissingInitialState,
        sCantImportInitialState,
        
        //StochHMM Transition
        
        //StochHMM Emission
        
        //StochHMM User Function Definitions
        
        //StochHMM Sequence
        sCantDigitizeSequence,
        sInvalidFormatting,
        
        //StochHMM Sequences
        sDifferenSizeSequences, 
        
        //StochHMM External Definitions
        
        //StochHMM SeqTracks
        
        //StochHMM SeqJobs
        
        //StochHMM Trellis
        
        
        //StochHMM Traceback
        sNoValidTraceback,
        sInvalidTraceback
        
        //StochHMM Output
    };
    
    
    //! \enum enum errPriority {ErrDebug, ErrNotice, ErrWarnings,ErrCritical, ErrEmergency};
    
    enum errPriority 
    {
        ErrDebug,
        ErrNotice,
        ErrWarnings,
        ErrCritical,
        ErrEmergency
    };
    
    
    
    //! \class stochErrInfo
    //! Error information about the thrown error
    class stochErrInfo{
    public:
        
        //!Create stochErrInfo
        inline stochErrInfo():sLine(0){};
        
        //!Create stochErrInfor with applicable information
        //!\param func  Function name that threw error
        //!\param file  File where error was thrown
        inline stochErrInfo(const char* func, const char* file, int ln=0): sFunction(func), sFile(file), sLine(ln){setTime();};
        
        
        
        //!Create stochErrInfor with applicable information
        //!\param info  Additional error information 
        //!\param func  Function name that threw error
        //!\param file  File where error was thrown
        inline stochErrInfo(const char* info, const char* func , const char* file,int ln=0): sFunction(func), sFile(file), sInfo(info), sLine(ln){setTime();};
        
        //Mutator
        void setTime();  //Get current time from system
        
        //!Set the time of error from character string
        //!\param tm  Time error occured
        inline void setTime(const char*  tm){sTime=tm;};
        
        //!Set the time of error from std::string
        //!\param tm  Time error occured
        inline void setTime(std::string& tm){sTime=tm;};
        
        //!Set the function name that threw the error with std::string
        //!\param fn Function that threw error
        inline void setFunction(std::string& fn){sFunction=fn;};
        
        //!Set the function name that threw the error with character string
        //!\param fn Function that threw error
        inline void setFunction(const char*  fn){sFunction=fn;};
        
        
        //!Set the file name that threw the error with std::string
        //!\param fn File where error was thrown
        inline void setFile(std::string& fl){sFile=fl;};
        
        //!Set the file name that threw the error with character string
        //!\param fn File where error was thrown
        inline void setFile(const char*  fl){sFile=fl;};
        
        
        //!Set the error information 
        //!\param info Additional Error Information
        inline void setInfo(std::string& info){sInfo=info;};
        
        //!Set the error information 
        //!\param info Additional Error Information
        inline void setInfo(const char*  info){sInfo=info;};
        
        std::string stringify();
        void print();
        
        //Accessor
        //!Get line where error occured
        inline int getLine(){return sLine;};
        
        
        //!Get std::string of time error occured
        //! \return std::string&
        inline const std::string& getTime(){return sTime;};
        
        //!Get Function that threw error
        //!\return std::string
        inline const std::string& getFunction(){return sFunction;};
        
        //!Get File that threw error
        //!\return std::string
        inline const std::string& getFile(){return sFile;};
        
        //!Get Additional information about the error
        //!\return std::string
        inline const std::string& getInfo(){return sInfo;};
        
    private:
        std::string sTime;
        std::string sFunction;
        std::string sFile;
        std::string sInfo;
        int         sLine;
        
    };
    
    
    //!\class stochErr : public std::exception
    //! Keeps error information when exception is thrown
    //! Allows for propagation of error, using C++ exception handling
    class stochErr: public std::exception {
    public:
        
        stochErr();
        stochErr(const errStoch, const char*, const char*, int);
        stochErr(const errStoch, const char*, const char*, const char*,int);
        //stochErr(stochErr&);
        
        ~stochErr() throw();
        
        const char* getPriority() const;
        const char* getError() const;
        const char* what() const throw();
        
        //!Get the Priority value of the error
        //! \return errPriority
        inline errPriority getPriorityValue() const{return sPriority;}
        
        //!Get the Error value of the Error
        //! \return errStoch
        inline errStoch  getErrorValue() const{return sErr;}
        
        stochErr operator= (const stochErr&);
        
        //!Get the error information for the ith thrown error
        //! \param it Index iterator for ith error
        inline stochErrInfo* operator[] (size_t it){return sInfo[it];};
        
        //!Get the error info for first error thrown
        inline stochErrInfo* front(){return sInfo.front();};
        
        //!Get the error info for last error thrown
        inline stochErrInfo* back(){return sInfo.back();};
        
        //!Add additional error information to back of the list
        inline void push_back(stochErrInfo* info){sInfo.push_back(info);};
        
        //!Add additional error information to front of the list
        inline void push_front(stochErrInfo* info){sInfo.push_front(info);};
        
        //! Remove the first error information
        inline void pop_front(){sInfo.pop_front();};
        
        //! Remove the last error information
        inline void pop_back() {sInfo.pop_back();};
        
        //! Get the amount of stochErrInfo
        inline size_t size(){return sInfo.size();};
        
    private:
        errStoch    sErr;  //Type of error
        errPriority sPriority; //The priority of the error
        std::deque<stochErrInfo*> sInfo;  //queue of error information
    };

}
#endif
