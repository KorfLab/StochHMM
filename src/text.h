//text.h

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
#ifndef TEXT_H
#define TEXT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

namespace StochHMM{
	
	/*! \defgroup Text  Text Handling
	 *  
	 *
	 */
	
	
    
    /*! \class stringList
     \brief Stringlist is an list of strings that contains parsed comments and 
     can be used to split string or remove whitespace from a string
     
     */

    class stringList{
    public:
        stringList();
        stringList(std::string&,std::string&,std::string&,bool);
        
        //MUTATORS
        //! Set remove whitespace flag. When splitting string it will remove Whitespace first and then split the string.
        inline void setRemoveWS(){removeWS=true;}; //!<
        
        
        //! Unset remove whitespace flag. When splitting string it will remove Whitespace first and then split the string.
        inline void unsetRemoveWS(){removeWS=false;};
        
        //! Set whitespace character. When splitting string it will remove defined whitespace characters before splitting the string.
        //! \param txt String of whitespace characters to remove from sequence
        inline void setWhitespace(const std::string& txt){whitespace=txt;};
        
        //! Set delimiter characters.  String will be split on these characters
        //! \param txt String of text-delimiters to use in splitting string
        inline void setDelimiters(const std::string& txt){delimiters=txt;};
        
        void getLine(std::istream&);
        void processString(std::string&); //Remove comment, whitespace, and split string
        void splitString(const std::string&); //Split string with delimiter defined in stringList
        void splitString(const std::string&, const std::string&); //Split string with given delimiter
        void splitString(const std::string&, size_t);
        void splitND(const std::string&, const std::string&);
        void removeLWS(const std::string&);
        void removeComments();
        
        //! Clears all values from the stringList, including comments, whitespace, and delimiters
        inline void clear(){line.clear();comment.clear();whitespace.clear();delimiters.clear();};
        
        bool fromTable(std::string&);
        bool fromTable(std::istream&);
        bool fromTrack(std::string&);
        bool fromTrack(std::istream&);
        bool fromAlpha(const std::string&,size_t);
        bool fromAlpha(std::istream&,size_t);
        bool fromTxt(std::string&);
        bool fromTxt(std::istream&);
        bool fromNext(std::string&);
        bool fromNext(std::istream&);
        bool fromDef(std::string&,std::string&, std::string&);
        bool fromDef(std::istream&,std::string&, std::string&);
        
        //! Add string to StringList
        //! \param txt String to add to the stringList
        inline void operator<< (std::string& txt){line.push_back(txt);};
        
        //! Add string to StringList
        //! \param txt String to add to the stringList
        inline void push_back  (std::string& txt){line.push_back(txt);};
        
        std::string pop_ith(size_t);
        
        //! Copy StringList
        //! \param lst stringList to copy
        inline void operator= (stringList& lst){line=lst.line; comment=lst.comment; whitespace=lst.whitespace; delimiters=lst.delimiters;};
        
        //ACCESSORS
        //! Access string at iterator 
        //! \param iter Position to return;
        inline std::string& operator[](size_t iter){return line[iter];};
        
        //! Return any comments parsed from the line
        inline std::string& getComment(){return comment;};
        
        //! Returns the list as a Vector of Doubles
        std::vector<double> toVecDouble();
        
        void toVecInt(std::vector<int>&);

        
        size_t indexOf(const std::string&);
        size_t indexOf(const std::string&,size_t);
        bool contains(const std::string&);
		bool containsExact(const std::string& txt);
        void print();
        std::string stringify();
        
        //! Return the amount of strings in the stringList 
        inline size_t size(){return line.size();};
    
        
    private:
        void _splitIndividual(std::string&, size_t);
        bool foundAlphaDelimiter(const std::string&);
        
        std::vector<std::string> line;
        std::string comment;
        
        bool removeWS;
        std::string whitespace;
        std::string delimiters;
    };


   
    
    stringList& splitString(const std::string&,const std::string&);
    
    void clear_whitespace(std::string&, std::string);
    void replaceChar(std::string&, char ,char);
    void removeLeadingWS(std::string&,const std::string&);
    std::string parseComment(std::string&,char);
    void getKeyValue(std::string&,std::string&,std::string&);

    std::string join(std::vector<int>&, char);
	std::string join(std::vector<size_t>&, char);
    std::string join(std::vector<short>&,char);
    std::string join(std::vector<double>&, char);
    std::string join(std::vector<std::string>&, char);

    std::string int_to_string(int);
    std::string int_to_string(size_t);
    std::string double_to_string(double);
    std::string double_to_string(float);
    
    bool stringToInt(std::string&, int&);
	bool stringToInt(std::string&, size_t&);
    bool stringToDouble(std::string&, double&);

    bool isNumeric(const std::string&);


    void split_line(std::vector<std::string>&, std::string& );

    stringList extractTag(std::string&);
        std::pair<size_t,size_t> balanced_brackets(const std::string&, const std::string& );
    std::pair<size_t,size_t> balanced_brackets(const std::string&, const std::string& , size_t );
    
    std::string slurpFile(std::string&);

    /*
     inline double convertToDouble (const std::string& s){
     istd::stringstream i(s);
     double x;
     if (!(i>>x)){
     throw BadConversion("convertToDouble(\"" + s "\")"));
     }
     return x;
     }
     */
	
	/**@}*/ //End of Text Handling Group
    
}
#endif /*TEXT_H*/
