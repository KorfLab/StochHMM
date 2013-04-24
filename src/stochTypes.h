//
//  stochTypes.h
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

#ifndef StochHMM_stochTypes_h
#define StochHMM_stochTypes_h
namespace StochHMM{
    //! \file stochTypes.h
	//! Stores enumerated types used in StochHMM library

    //Enumerated Transition identifiers
    //!\enum enum transType {STANDARD , USER , INTERNAL , LEXICAL};
    //!Types of Transitions
    //! STANDARD = Transition has log probability value
    //! USER = Transition is duration dependent and calculated by user defined distribution
    //! INTERNAL = Transition is duration dependent and calculated by internal defined distibution
    //! LEXICAL =  Transition is dependent upon the preceeding sequence
    enum transType {STANDARD , DURATION , LEXICAL};
    
    
    //!\enum tracebackIdentifier { DIFF_STATE , STATE_NAME , STATE_LABEL , STATE_GFF , START_INIT };
    //! Traceback identifier describe how stochHMM will traceback when using duration dependent distributions
    //! Types of Tracebacks
    //! DIFF_STATE = Traceback current state until a different state is encountered
    //! STATE_NAME = Traceback until STATE_NAME changes from current state, essentially the same as DIFF_STATE
    //! STATE_LABEL = Traceback until a State with a different label is encountered
    //! STATE_GFF = Traceback until a State witha a different GFF tag is encountered
    //! START_INIT = Traceback until the start of the sequence
    enum tracebackIdentifier { DIFF_STATE , STATE_NAME , STATE_LABEL , STATE_GFF , START_INIT };
    
    
    //!enum combineIdentifier { FULL , STATENAME , STATELABEL , STATEGFF};
    //!Describes how a traceback will be processed
    //! Types of Combine Identifiers
    //! FULL = No traceback editing will occur
    //! STATENAME = Edit out all the states that aren't of given State Name
    //! STATELABEL= Edit out all the states that aren't of given State Label
    //! STATEGFF  = Edit out all the states that aren't of given GFF Tag
    enum combineIdentifier { FULL , STATENAME , STATELABEL , STATEGFF};

    //Enumerated Trellis Types
    //!\enum trellisType {SIMPLE, STOCH, NTH};
    //!Defines what type of trellis cells to use
    //!Types of trellis
    //! SIMPLE = Trellis will only contain single viterbi, forward, and backward values
    //! STOCH  = Trellis will contain multiple viterbi,forward,backward values and traceback probabilities for stochastic tracebacks
    //! NTH = Trellis will contain N viterbi values for top N tracebacks 
    enum trellisType {SIMPLE, STOCH, NTH};
    
    //!\enum decodingType {VITERBI, POSTERIOR};
    //!Type of decoding to perform.
    //! VITERBI = Traceback performed using viterbi value
    //! FORWARD = Traceback performed using stochastic forward value
    //! POSTERIOR = Traceback performed using posterior value
    enum decodingType {VITERBI, FORWARD, POSTERIOR};

    //Enumerated Emission Track types
    //!Track types
    //! UNDEFINED = NOT DEFINED BY USER
    //! ALPHA_NUM = Alphabet is letter or word based (discrete)
    //! REAL = Real number values (continuous)
    //! EXTERNAL = ???
    enum trackType {UNDEFINED, ALPHA_NUM , REAL , EXTERNAL};


    //Enumerated Unknown Emission Character Probability Types
    
    //!enum unknownCharScoringType { DEFINED_SCORE, AVERAGE_SCORE, LOWEST_SCORE, HIGHEST_SCORE, NO_SCORE};
    //! How to score unknown or ambiguous characters when encountered in emission or transition
    //! Scoring types:
    //! DEFINED_SCORE = user-defined score
    //! AVERAGE_SCORE = average of all possible scores
    //! LOWEST_SCORE = lowest of the possible scores
    //! HIGHEST_SCORE = highest of the possible scores
    //! NO_SCORE = No score will be defined (Produces error if unknown alphabet encountered)
    enum unknownCharScoringType { DEFINED_SCORE, AVERAGE_SCORE, LOWEST_SCORE, HIGHEST_SCORE, NO_SCORE};
    
    
    //enum valueType {PROBABILITY, LOG_ODDS, COUNTS, LOG_PROB, PERCENTAGE};
    //!Value types that are provided by user
    //!PROBABILITY = basic probability between [0,1]
    //!LOG_ODDS = log odds score
    //!COUNTS = word counts
    //!LOG_PROG = log2 value of probability
    //!PERCENTAGE = [0,100] or 100*Probability
    enum valueType {PROBABILITY, LOG_ODDS, COUNTS, LOG_PROB, PERCENTAGE}; 


}
#endif
