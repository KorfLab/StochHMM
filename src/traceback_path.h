///traceback_path.h

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

#ifndef TRACEBACK_PATH_H
#define TRACEBACK_PATH_H

#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdint.h>
#include <stdlib.h>
#include <cassert>
#include "text.h"
#include "options.h"
#include "hmm.h"
#include "stochMath.h"
namespace StochHMM{

    //! \struct gff_feature
    //!Each gff_feature represents a single GFF line
    struct gff_feature{
        std::string seqname;    //Column 1
        std::string source;     //Column 2
        std::string feature;    //Column 3
        size_t start;              //Column 4
        size_t end;                //Column 5
        char score;             //Column 6
        char strand;            //Column 7
        char frame;             //Column 8
        std::string attribute;  //Column 9
    };


    //! Perform traceback of traceback table
    //! Stores one traceback path for a sequence
    class traceback_path{
    public:
        traceback_path(model*);
        
		//!Add state to traceback
        void push_back(int);
		
		//!Erase the traceback
        void clear();
		
		//!Get the size of traceback
		//!\return size_t
        size_t size() const;
		
		//!Returns the state index at a given position (it) within the traceback sequence
        inline int val(size_t it){
			if (it>=trace_path.size()){
				std::cerr << "Out of Range index\n ";
				exit(2);
			}
			return trace_path[it];
		}
		
		//!Get the model used for the decoding
        //! \return model
		inline model* getModel() const {return hmm;};
		inline void setModel(model* mdl){hmm=mdl;};
        
		
		//!Print the path to file stream
        void fprint_path(std::ofstream&);
        
		//! Get traceback as vector of state labels
		//! \param [out] std::vector<std::string> Vector of Labels
        void label(std::vector<std::string>&);
		
		//! Get traceback as string of state labels
        void label(std::string&);
        
		
		//! Get traceback as vector of gff_features
		//! \param[out] pth reference to vector of gff_feature
		//! \param[in] sequenceName Name of Sequence to be used in GFF
        void gff(std::vector<gff_feature>&,std::string&);
		
		//!Get names of traceback path
		//!\param [out] pth vector<string>
        void name(std::vector<std::string>&);
		
		//! Get the path to std::vector<int>
		//! \param [out] pth std::vector<int> that represents path
        void path(std::vector<int>&);
        
		//! Print the traceback path as path to stdout using cout
        //! Path numbers correspond to state index in model
		void print_path() const ;
		
		//! Print the traceback path as state labels
		//! State labels
        void print_label() const ;
		
		//!Outputs the gff formatted output for the traceback to stdout
		//!Allows user to provide additional information, that may be
		//!pertinent to stochastic tracebacks
		//!\param[in] sequence_name  Name of sequence used
		//!\param[in] score score to use in the GFF output
		//!\param[in] ranking  Rank of traceback
		//!\param[in] times Number of times that traceback occurred
		//!\param[in] posterior Posterior probability score
        void print_gff(std::string,double,int,int,double) const ;
        
		//!Outputs the gff formatted output for the traceback
		void print_gff(std::string) const ;
        
		//!Get the score that is associated with the traceback
        inline double getScore(){
            return score;
        }
        
		//!Set the score for the traceback;
        inline void setScore(double scr){
            score = scr;
            return;
        };
        
        //double path_prob (const HMM&, sequence&);  //Need to rewrite function
        inline int operator[](size_t val) const {return trace_path[val];};
        bool operator== (const traceback_path&) const;
        bool operator<  (const traceback_path&) const;
        bool operator>  (const traceback_path&) const;
        bool operator<= (const traceback_path&) const;
        bool operator>= (const traceback_path&) const;
    private:
        model* hmm;
        std::vector<int> trace_path;
        double score;
    };

    //Rows are the Sequence Position
    //Columns are the States
    //! \def std::vector<std::vector<int> > heatTable;
    //! Table for heat map data generated from multiple tracebacks
    typedef std::vector<std::vector<int> > heatTable;

    //! \class multiTraceback
    //! Contains multiple tracebacks.  Will store them in sorted unique list (sorted in order of number of occurances);
    class multiTraceback{
    public:
        multiTraceback();
        ~multiTraceback();
        
        //Iteratorate through map;
        void begin();
        void end();
        void operator++();
        void operator--();
        void operator=(size_t);
        
        
        void print_path();
        void print_label();
        void print_gff(std::string&);
        void print_hits();
        
        
        //Access the data at a point
        traceback_path path();
        int counts();
        traceback_path operator[](size_t);
        
        
        //Assign
        void assign(traceback_path&);
        
        //Finalize multiTraceback (Sort and setup Iterators);
        void finalize();
            
        inline void clear(){paths.clear();return;};
        inline size_t size(){return paths.size();};
        
        heatTable* get_hit_table();

    private:
        size_t vectorIterator;
        size_t maxSize;
        bool isFinalized;
        std::vector<std::map<traceback_path,int>::iterator> pathAccess;
        std::map<traceback_path,int> paths;
        heatTable* table;
    };
     
    bool sortTBVec(std::map<traceback_path,int>::iterator, std::map<traceback_path,int>::iterator);
}
#endif /*TRACEBACK_PATH_H*/



