//text.cpp
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

#include "text.h"
namespace StochHMM{
    
    //! Create a new stringList
    //! Remove whitespace set to true by default.
    stringList::stringList(){
        removeWS=true;
        return;
    }
    
    //! Create a new StringList
    //! \param txt String to be split
    //! \param ws Whitespace characters to remove from string before splitting
    //! \param del Delimiter characters to use to split string
    //! \param remove True will remove whitespace, False will leave whitespace
    
    stringList::stringList(std::string& txt, std::string& ws, std::string& del,bool remove):whitespace(ws),delimiters(del),removeWS(remove)
    {
        splitString(txt);
        return;
    }
    
    //! Searches for the string and returns bool if found or not
    //! \param txt String to search stringList
    bool stringList::contains(const std::string& txt){
        for (size_t i=0;i<line.size();i++){
            if (line[i].find(txt) != std::string::npos){
                return true;
            }
        }
        return false;
    }
    
    //FIXME: Fix return value to size_t max not -1
    //! Searches the stringList for matching string and returns the index position of first match 
    //! \param txt String to search stringList for
    size_t stringList::indexOf(const std::string&txt) {
        for (size_t i=0;i<line.size();i++){
            if (line[i].find(txt) != std::string::npos){
                return i;
            }
        }
        return -1;
    }
    
    //! Searches the stringList for matching string and returns the index position of a given string from the starting position
    //! \param txt String to search stringList for
    //! \param pos Position to start search from
    size_t stringList::indexOf(const std::string&txt,size_t pos){
        for (size_t i=pos;i<line.size();i++){
            if (line[i].find(txt) != std::string::npos){
                return i;
            }
        }
        return -1;
    }
    
    //! Removes Comments, Removes Whitespace, and splits the string
    //! Whitespace and Split delimiters are required to be previously set in stringList
    void stringList::processString(std::string& txt){
        line.clear();
        comment.clear();
        
        //Extract Comments
		comment=parseComment(txt, '#');
        
        //Remove Whitespace
        if (removeWS){
            clear_whitespace(txt, whitespace);
        }
		
        splitString(txt);
        
        return;
    }
    
    //! Split string using delimiter
    //! Splits the string up according to character delimiters defined in stringList.
    //! The char delimiters are deleted from the returned value.  Delimiters allows you supply multiple
    //! delimiters in a single string.
    //! \param txt String to be split
    void stringList::splitString(const std::string& txt){
        
        line.clear();
        comment.clear();
        
		size_t found=-1;
		size_t initial=-1;
        do {
            found++;
            initial=found;
            
            found=txt.find_first_of(delimiters.c_str(),found);
            std::string st=txt.substr(initial,found-initial);
            
            if (st.size()>0){
				line.push_back(st);
			}
            
        } while (found!=std::string::npos);
        
        return;
    }
    
    //! Split string using delimiter
    //! Splits the string up according to character delimiters supplied in the delimiter string.
    //! The char delimiters are deleted from the returned value.  Delimiters allows you supply multiple
    //! delimiters in a single string.
    //! \param txt String to be split
    //! \param del Delimiters to use to split evaluated as single characters, not whole string
    void stringList::splitString(const std::string& txt, const std::string& del){
        
        delimiters=del;
        splitString(txt);
        return;
    }
    
    void stringList::splitString(const std::string& txt,size_t charSize){
        for(size_t i=0;i<txt.size();i+=charSize){
            line.push_back(txt.substr(i,charSize));
        }
        return;
    }
    
    //! Non-destructive string split
    //! Splits the string up using an additional string as the delimiter. Unlike
    //! stringSplit the delimiter isn't deleted from the returned value;
    //! \param txt  String to be split
    //! \param del  Delimiter string to use to split
    void stringList::splitND(const std::string& txt,const std::string& del){
		line.clear();
        comment.clear();
        size_t initial=txt.find(del);
        
        do {
            size_t found=txt.find(del,initial+1);
            std::string st=txt.substr(initial,found-initial);
            initial=found;
            if (st.size()>0){
				line.push_back(st);
			}
            
        } while (initial!=std::string::npos);
        
        return;
    }
    
    
    //!Remove leading whitespace
    //!Removes the leading whitespace from the sequence
    //! \param ws String of whitespace characters to remove
    void stringList::removeLWS(const std::string& ws){
        for(size_t i=0;i<this->size();i++){
            removeLeadingWS((*this)[i],ws);
        }
        return;
    }
    
    
    //! Remove whitespace characters '\\t' and then splits string by delimiters ':' or '\\n'
    //! \param txt String to be split 
    bool stringList::fromTxt(std::string& txt){		
        setWhitespace("\t");
        setDelimiters(":\n");
        processString(txt);
        
        if (size()>0){
            return true;
        }
        else{
            return false;
        }
    }
    
    //! Remove whitespace characters '\\t' and then splits string by delimiters ':' or '\\n'
    //! Uses std::istream::getline(...,\\n) and splits string
    //! \param in std::istream reference 
    bool stringList::fromTxt(std::istream& in){
        std::string temp;
        getline(in, temp, '\n');
        return fromTxt(temp);
    }
    
    
    //! Remove whitespace characters '\t' and then splits string by delimiters ':', '\n', comma or space
    //! \param txt String to be split 
    bool stringList::fromTrack(std::string& txt){
		setWhitespace("\t");
        setDelimiters(":, \n");
        processString(txt);
        
        if (size()>0){
            return true;
        }
        
        return false;
    }
    //! Remove whitespace characters '\t' and then splits string by delimiters ':', '\n', comma or space
    //! \param in std::istream reference
    bool stringList::fromTrack(std::istream& in){
        std::string temp;
        getline(in, temp, '\n');
        return fromTrack(temp);
    }
    
    //! Remove whitespace characters '\t','\n' or space and then splits string by delimiters ':' or comma
    //! \param txt String to be split 
    bool stringList::fromNext(std::string& txt){
		setWhitespace("\t\n ");
        setDelimiters(":,");
        processString(txt);
        
        if (size()>0){
            return true;
        }
        
        return false;
    }
    
    //! Remove whitespace characters '\t','\n' or space and then splits string by delimiters ':' or comma
    //! \param txt std::istream reference
    bool stringList::fromNext(std::istream& in){
        std::string temp;
        getline(in, temp, '\n');
        return fromNext(temp);
    }
    
    
    //! Remove whitespace characters '\t' and then splits string by delimiters ':' or '\n'
    //! \param txt String to be split
    bool stringList::fromAlpha(const std::string& txt,size_t alpha){
        line.clear();
		comment.clear();
		
		if (foundAlphaDelimiter(txt)){
            splitString(txt,";, \t");
        }
        else{
            splitString(txt,alpha);
        }
        
        if (this->size()>0){
            return true;
        }
        
        return false;
    }
    
    //! Returns true if text delimiters ":,[space]\\t" are found in the string
    //! \param txt String used to find delimiter
    bool stringList::foundAlphaDelimiter(const std::string& txt){
        for(int i=0;i<txt.size();i++){
            switch(txt[i]){
                    case ':':
                    return true;
                    case ',':
                    return true;
                    case ' ':
                    return true;
                    case '\t':
                    return true;
                default:
                    break;
            }
            
        }
        return false;
    }
    
    

    //! Remove whitespace characters '\t' and then splits string by delimiters ':' or '\n'
    //! \param in std::istream stream
    //! \param alpha Size of alphabet
    bool stringList::fromAlpha(std::istream& in, size_t alpha){
        std::string temp;
        getline(in, temp, '\n');
        return fromAlpha(temp,alpha);
    }
    
    
    //! Splits a definition from model file
    //! \param txt string to split
    //! \param ws Whitespace characters to use
    //! \param del Delimiters characters to use
    bool stringList::fromDef(std::string& txt, std::string& ws, std::string& del){
        whitespace=ws;
        delimiters = del;
        processString(txt);
        
        if (size()>0){
            return true;
        }
        else{
            return false;
        }
    }
    
    
    //! Splits a definition from model file
    //! \param in std::istream to get line from
    //! \param ws Whitespace characters to use
    //! \param del Delimiters characters to use
    bool stringList::fromDef(std::istream& in, std::string& ws, std::string& del){
        whitespace=ws;
        delimiters = del;
        std::string temp;
        getline(in, temp, '\n');
        return fromDef(temp, ws, del);
    }
    
    
    std::string stringList::pop_ith(size_t pos){
        if (pos>=line.size()){
            return "";
        }
        
        std::string temp=line[pos];
        line.erase(line.begin()+pos);
        return temp;
    }
    
    
    //! Returns the values in the stringList as std::vector of doubles
    std::vector<double> stringList::toVecDouble(){
        std::vector<double> temp;
        for(size_t iter=0;iter<line.size();iter++){
            
            double val;
            
            stringToDouble(line[iter], val);
            temp.push_back(val);
        }
        return temp;
    }
    
    //! Returns the values in the stringList as std::vector of integers
    void stringList::toVecInt(std::vector<int>& ret_val){
        
        for(size_t iter=0;iter<line.size();iter++){
            
            int val(0);
            
            if (!stringToInt(line[iter], val)){
                std::cerr << "Couldn't convert " << line[iter] << " to an integer\n";
            }
            else{
                ret_val.push_back(val);
            }
            
        }
        
        return;
    }

    
    
    //! Print each line to stdout
    void stringList::print(){
        for(size_t i=0;i<line.size();i++){
            std::cout << line[i] << std::endl;
        }
        std::cout << "#" << comment << std::endl;
    }
    
    
    //! Joins the stringList into string using "\t" and return string
    //! \return string of stringList joined by "\\t"
    std::string stringList::stringify(){
        std::string output=join(line,'\t');
        if (comment.size()>0){
            output+="#"+ comment;
        }
        
        return output;
    }
    
    
    //!Parses out the comments and stores the comment delimited by "#" 
    //!Comment can be accessed by command getComment();
    void stringList::removeComments(){
        for(size_t iter=0;iter<line.size();iter++){
            comment+=parseComment(line[iter], '#');
        }
    }
    
    
    //! Find first comment character and then return everything following the character
    std::string parseComment(std::string& txt, char commentChar){
        std::string comment;
        size_t commentPos=txt.find_first_of(commentChar);
        if (commentPos!=std::string::npos){
            comment=txt.substr(commentPos);
            txt=txt.substr(0,commentPos);
        }
        return comment;
    }

    //! Given a string, and a white space character, it will remove all the whitespace characters from the string
    //! \param input String to remove whitespace from
    //! \param white String of whitespace characters to remove
    void clear_whitespace(std::string &input,std::string white){
        size_t found;
        found=input.find_first_of(white);
        //found=input.find_first_of("\t\n ");
        while(found!=std::string::npos){
            input.erase(found,1);
            //found=input.find_first_of("\t\n ");
            found=input.find_first_of(white);
        }
        return;
    }
    
    //! Removes leading whitespace characters from a string
    //! \param txt String user wants to remove whitespace from
    //! \param ws  String containing whitespace characters to remove 
    void removeLeadingWS(std::string& txt,const std::string& ws){
        size_t start = txt.find_first_not_of(ws);
        if (start==std::string::npos){
            txt="";
        }
        else{
            txt=txt.substr(start);
        }
        return;
    }
    
    
    //! Parses key and value from a line
    //! Where key is delimited by <<KEY>> = Value
    //! \param txt String to extract key value from
    //! \param key String to assign key to
    //! \param value String to assign value to
    void getKeyValue(std::string& txt,std::string& key,std::string& value){
        size_t found=txt.find("<<");
        
        if (found==std::string::npos){
            removeLeadingWS(txt,"\t \n");
            value=txt;
            return;
        }
        else{
            size_t ending=txt.find(">>");
            if (ending!=std::string::npos){
                key=txt.substr(found+2,ending-(found+2));
                stringList lst;
                lst.splitString(txt,"=");
                removeLeadingWS(lst[1],"\t \n");
                value= lst[1];
                return;
            }
            else{
                std::cerr << "Missing closing brackets on Key\n";
                exit(1);
            }
            
        }
    }
    
    
    //! Splits string using delimiters and return stringList
    stringList& splitString(const std::string& txt, const std::string& delimiters){
        static stringList lst;
        lst.clear();
        lst.splitString(txt,delimiters);
        return lst;
    }
    
    
    //! Replace a given character with another character in a string
    //! \param txt String to use have characters replaced
    //! \param ch Character to search string for
    //! \param replaceCh Character to replace found ch with
    void replaceChar(std::string& txt, char ch, char replaceCh){
        size_t found = txt.find(ch);
        while(found!=std::string::npos){
            txt[found]=replaceCh;
            found=txt.find(ch);
        }
        return;
    }
    
    //! Converts a vector of ints into a string delimited by a character c
    //! \param input Vector of integers to be converted
    //! \param c Character to use as a delimiter
    std::string join(std::vector<int> &input, char c){
        std::string out;
        if (input.size()==0){
            out="";
            return out;
        }
        else if (input.size()==1){
            out=int_to_string(input[0]);
            return out;
        }
        else{
            out=int_to_string(input[0]);
            for(int i=1;i<input.size();i++){
                out+=c;
                out+=int_to_string(input[i]);
            }
            return out;
        }
    }

    //! Convert an integer to a string
    //! \param input Integer you want to convert to string;
    std::string int_to_string(int input){
        std::stringstream ss;
        ss << input;
        std::string s=ss.str();
        return s;
    }
    
    
    //! Convert an size_t to a string
    //! \param input Integer you want to convert to string;
    std::string int_to_string(size_t input){
        std::stringstream ss;
        ss << input;
        std::string s=ss.str();
        return s;
    }
    
    
        
    //! Convert a double to a string
    //! \param input Double you want to convert to a string
    std::string double_to_string(double input){
        std::stringstream ss;
        ss << input;
        std::string s=ss.str();
        return s;
    }
    
    //! Convert a double to a string
    //! \param input Double you want to convert to a string
    std::string double_to_string(float input){
        std::stringstream ss;
        ss << input;
        std::string s=ss.str();
        return s;
    }
    
    
    //!Convert string to integer
    //!\param txt Text representation of integer
    //!\param val Integer to be assigned
    //!\return true if conversion is valid
    //!\return false if conversion can't be performed
    bool stringToInt(std::string& txt, int& val){
        std::istringstream input(txt);
        if (!(input >> val)){
            return false;
        }
        
        return true;
    }
    
    
    //!Convert string to double
    //!\param txt Text representation of double
    //!\param val Integer to be assigned
    //!\return true if conversion is valid
    //!\return false if conversion can't be performed
    bool stringToDouble(std::string& txt, double& val){
        std::istringstream input(txt);
        if(!(input >> val)){
            return false;
        }
        return true;
    }
    
    

        
    //! Converts a vector of shorts into a string delimited by a character c
    //! \param input Vector of shorts to be converted
    //! \param c Character to use as a delimiter
    std::string join(std::vector<short> &input, char c){
        std::string out;
        if (input.size()==0){
            out="";
            return out;
        }
        else if (input.size()==1){
            out=int_to_string(input[0]);
            return out;
        }
        else{
            out=int_to_string(input[0]);
            for(int i=1;i<input.size();i++){
                out+=c;
                out+=int_to_string(input[i]);
            }
            return out;
        }
    }

    //! Converts a vector of doubles into a string delimited by a character c
    //! \param input Vector of doubles to be converted
    //! \param c Character to use as a delimiter
    std::string join(std::vector<double> &input, char c){
        std::string out;
        if (input.size()==0){
            out="";
            return out;
        }
        else if (input.size()==1){
            out=double_to_string(input[0]);
            return out;
        }
        else{
            out=double_to_string(input[0]);
            for(int i=1;i<input.size();i++){
                out+=c;
                out+=double_to_string(input[i]);
            }
            return out;
        }
    }
    
    //! Converts a vector of strings into a string delimited by a character c
    //! \param input Vector of strings to be converted
    //! \param c Character to use as a delimiter
    std::string join(std::vector<std::string> &input, char c){
        std::string out;
        size_t sz=input.size();
        if (sz==1){
            out = input[0];
        }
        else if(sz>1){
            out+=input[0];
            for(size_t i=1;i<sz;i++){
                out+= c + input[i];
            }
        }
        
        return out;
    }
    
    
    //! Splits a line into a vector of string using delimiters ' \",[space]\\n\\t'
    //! \param line vector of strings to split input into
    //! \param input String to be split using delimiters ' \",[space]\\n\\t'
    void split_line(std::vector<std::string> &line,std::string &input){
        
        //split line accoring to delimiters;
        size_t found=input.find_first_of("\", \n\t");
        while(found!=std::string::npos){
            if (found>0){
                line.push_back(input.substr(0,found));
                //cout << line.back() << endl;
                input=input.substr(found+1);
                //cout << input << endl;
            }
            else{
                input.erase(found,1);
                //cout << input <<endl;
            }
            found=input.find_first_of("\", \n\t");
            //cout <<input <<endl;
        }
        if (input.size()>0){
            line.push_back(input);
        }
        return;
    }
    
    //! Parse a line and extract a bracketed tag from the model file
    //! Returns a stringList which contains the tag split using ":\\t[space]"
    //! \param txt String to be have tag extracted from
    stringList extractTag(std::string& txt){
        stringList lst;
        std::pair<size_t,size_t> tagCoord = balanced_brackets(txt,"[]");
        if (tagCoord.first!=tagCoord.second){
            std::string tag = txt.substr(tagCoord.first+1,tagCoord.second-tagCoord.first-1);
            txt.erase(tagCoord.first,tagCoord.second-tagCoord.first+1);
            lst.splitString(tag,":\t ");
        }
        return lst;
    }

    //! Returns a pair of size_t values that describe the coordinates of between brackets
    //! \param text  String to use to search for brackets
    //! \param brackets String of length two containing opening and closing bracket
    //! \param offset Offset of where to start searching for balanced brackets
    std::pair<size_t,size_t> balanced_brackets(const std::string& text, const std::string& brackets, size_t offset){
        char opening = brackets[0];
        char closing = brackets[1];
        
        
        int currentTotal(0);
        
        size_t start;
        size_t found;
        found=text.find_first_of(opening,offset);
        
        if (found!=std::string::npos){
            start=found;
            currentTotal++;
        }
        else{
            return std::make_pair(0,0);
        }
        
        while(currentTotal!=0){
            found++;
            found=text.find_first_of(brackets,found);
            if (found!=std::string::npos){
                if (text[found]==opening){
                    currentTotal++;
                }
                else if (text[found]==closing){
                    currentTotal--;
                }
            }
            else{
                return std::make_pair(0,0);
            }
        }
        
        return std::make_pair(start,found);
    }
    
    
    //! Returns a pair of size_t values that describe the coordinates of between brackets from start of the string
    //! \param text  String to use to search for brackets
    //! \param brackets String of length two containing opening and closing bracket
    std::pair<size_t,size_t> balanced_brackets(const std::string& text, const std::string& brackets){
        return balanced_brackets(text,brackets,0);
    }



    //! Is the value of string numeric
    //! \param str  String to determine if it is numeric
    bool isNumeric(const std::string& str){
        size_t found;
        found = str.find_first_not_of("0123456789.-eE");
        
        if (found!=std::string::npos){
            return false;
        }
        else{
            return true;
        }
    }
    
    //!Slurp a file into a string 
    //! \param file Filename
    //! \return string that contains complete file
    std::string slurpFile(std::string& file){
        std::ifstream in(file.c_str(),std::ifstream::in);
        if (!in.good()){
            std::cerr << "File doesn't exist:" << file;
            exit(1);
        }
            
        std::stringstream sstr;
        sstr << in.rdbuf();
        return sstr.str();
    }

    
}
