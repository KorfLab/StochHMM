//track.h
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

#ifndef TRACK_H
#define TRACK_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include "text.h"
#include "stochTypes.h"
#include "userFunctions.h"
#include "stochMath.h"
#include <limits>
#include <stdlib.h>

namespace StochHMM{


    class track;
    class tracks;
    
    
    /*! \class ambigCharacter
     \brief Define the ambiguous characters symbol and index number for digitizing ambiguous characters in the sequence
     For example in DNA N = [ACGT] = [0,1,2,3]

     */
    class ambigCharacter{
    public:
        ambigCharacter(track*, std::string&, std::vector<std::string>& ); //track, ambiguous character, unambiguous characters
        
        
        //! Get Ambiguous String/Alphabet character symbol
        //! \return std::string Symbol of character/word
        inline std::string getSymbol(){return symbol;};
        
        
        //FIXME: should be vector of shorts unless I expand from 256 to larger int
        //! Get the characters that the characters defines.
        //! For example in DNA N = [ACGT] = [0,1,2,3]
        //! \return std::vector<int> Digitized value of characters that are represented by the given symbol
        inline std::vector<int>& getDef(){return setDefinition;};
        
    private:
        std::string symbol;
        std::vector<int> setDefinition;
    };
    
       
    //! \class track 
    //! Defines types of data (real-value, text-sequence) used in the model
    //! and the alphabet that a text-sequence uses.  Tracks are used to digitize
    //! the sequence before decoding in HMM

    class track {
    public:
        track();
        
        //TODO: Complete definition of constructor
        track(TrackFuncs*); 
        
        friend class state;
        friend class model;
        friend class tracks;
        
        //MUTATOR
        bool parse(std::string&);
        bool parseAmbiguous(std::string&);
        
        //!Set the name of the track
        //! \param nm Name of the track
        inline void setName(std::string& nm){name=nm;};
        
        //!Set the Description of the track
        //! \param desc Description of track
        inline void setDescription(std::string& desc){description=desc;};
        
        //!Set the integer index value of the track in tracks
        //! \param indx  User defined index value;
        inline void setIndex(size_t indx){
            if (indx<std::numeric_limits<size_t>::max()){
                trackIndex=indx;
                return;
            }
            std::cerr << "Track index: " << indx << " is OUT_OF_RANGE\n";
            exit(1);
        }
            
        
        //!Set the alphabet type of the track (Real or alphanum)
        //! \param typ enum trackType
        inline void setAlphaType(trackType typ){alpha_type=typ;};
        
        bool addAlphabetChar(std::string&);
        bool addAlphabetChar(std::vector<std::string>& , std::vector<std::string>&);
        bool addAlphabetChar(const char *); 
        
        void addComplement(std::string&, std::string&);
        void addComplement(const char *, const char *);
        
        //! Set ambiguous character flag to true
        //! This will allow ambiguous characters to be processed in sequence
        //! Without this flag, only strict track characters or values are allowed
        inline void setAmbiguous(){ambiguous=true; return;};
        
        
        void addAmbiguous(std::string&,std::vector<std::string>&);
        
        
        //ACCESSOR
        
        //! Get the name of the track
        //! \return std::string Name of the track
        inline std::string getName(){return name;};
        
        //! Get the description of the track
        //! \return std::string Description of the track
        inline std::string getDescription(){return description;};
        
        //! Get the index of the track
        //! \return int Index of the track
        inline size_t getIndex(){
            if (trackIndex<std::numeric_limits<size_t>::max()){
                return trackIndex;
            }
            std::cerr << "Track Index is not set in track.  Set the track index with setIndex(size_t indx) before calling.\n";
            exit(1);
        };
        
        //! Get the alphabet type of the track
        //! \return trackType Alphabet type of the track
        //! \sa enum trackType
        inline trackType getAlphaType(){return alpha_type;}
        
        //! Get the number of characters defined in the track
        //! \return size_t Number of characters/words defined in the track
        inline size_t getAlphaSize(){return alphabet.size();};
        
        //! Get the size of the largest alphabet word
        //! \return size_t
        inline size_t getAlphaMax(){return maxSize;};
        
        
        std::string getAlpha(int);
        
        
        int symbolIndex(std::string&);
        
        int getComplementIndex(int val);
        int getComplementIndex(std::string&);
        
                
        std::string getComplementSymbol(std::string& character);
        std::string getComplementSymbol(int value);
        
        inline bool isComplementDefined(){return complementSet;}
        
        
        //! Check to see if ambiguous flag is set for the track
        //! \return true if ambiguous flag is set to handle ambig. characters
        //! \return false if not set
        inline bool isAmbiguousSet(){return ambiguous;};
        
        
        //! Check to see if the track is AlphaNum type and not a REAL Track
        //! \return true if it is AlphaNum type
        //! \return false if it is a REAL type
        inline bool isAlpha(){if (alpha_type == ALPHA_NUM){return true;} else {return false;}};
        
        //! Get the number of ambiguous characters that are defined
        //! \return size_t Number of ambiguous characters/words defined
        inline size_t getAmbiguousSize(){return ambiguousSymbols.size();};
        
        std::string getAmbiguousCharacter(int);
        
        
        //! Get the indices of characters that an ambiguous character represents
        //! \return std::vector<int>  
        inline std::vector<int> getAmbiguousSet(int val){return ambiguousSymbols[abs(val)-1].getDef();};
        
        void print();
        std::string stringify();
        std::string stringifyAmbig();
        std::string convertIndexToWord(size_t,size_t);

        
        //!Check if the TrackFunction is defined for this track
        //!\return true if the track has a trackFunc defined
        inline bool isTrackFuncDefined(){return trackFunctionDefined;};
        
        //! Get name of TrackFunc defined for track
        //! \return std::string Name of trackFunc defined
        inline std::string getTrackFunction(){return trackFunction;};
        
        //! Get name of Track to use for trackFunc
        inline std::string getTrackToUse(){return trackToUse;};
        
        
    private:
        std::string name;	/* Track Name */
        std::string description;	/* Track Desc */
        size_t trackIndex;     /*track number*/
        
        trackType alpha_type;	/* Track Type 1=string, 2=real_number  0=uninitialized*/
        
        //! Track Functions for defining Real Number Tracks
        bool trackFunctionDefined;
        bool complementSet;
        
        std::string trackToUse;
        std::string trackFunction;
        
        std::vector<std::string> alphabet;  //Contains the corresponding symbol,letter, word that is referenced in the seq by index
        std::map<int,int> complementAlphabet;
        
        size_t maxSize;  //Maximum size of the alphabet words used.
        bool ambiguous;
        int defaultAmbiguous;
        
        //Contains the ambiguous characters defined by user corresponding to position in the
        //array. Where index 0=-1, 1=-2... so on.
        std::vector<ambigCharacter> ambiguousSymbols;
        
        std::map<std::string,int> symbolIndices;
        
        void _splitAmbiguousList(std::vector<std::pair<std::string ,std::vector<std::string> > >&, const std::string&);
    };

    
    class tracks{
    public:
        
        //MUTATOR
        void push_back(track*);
        
        //ACCESSOR
        size_t indexOf(const std::string&);
        size_t size(){return trks.size();};
        track* getTrack(const std::string&);
        bool isTrackDefined(const std::string&);
        track* operator[](size_t i){return trks[i];};
        
        void print();
        std::string stringify();
        
    private:
        std::vector<track*> trks;
        std::map<std::string,size_t> index;
    };




}
#endif /*TRACK_H*/