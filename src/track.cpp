//track.cpp
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

#include "track.h"
namespace StochHMM{
    
    // FIXME:  Check if ambiguous before allowing to get character.
    //!Create an ambiguous character
    //! For example in DNA N = [ACGT] = [0,1,2,3]
    //! \param tr Track to use for ambiguity code
    //! \param ambChar String representation of the character/symbol/word for ambiguous character
    //! \param defs vector of strings where each string is a non-ambiguous character defined in the track
    ambigCharacter::ambigCharacter(track* tr, std::string& ambChar, std::vector<std::string>& defs){
        symbol=ambChar;
        for(size_t i=0;i<defs.size();i++){
            setDefinition.push_back(tr->symbolIndex(defs[i]));
        }
        sort(setDefinition.begin(),setDefinition.end());
        return;
    }
	
    
    //!Create a track
    track::track(){
        alpha_type=UNDEFINED;
        trackIndex=std::numeric_limits<size_t>::max();
        ambiguous=false;
        defaultAmbiguous=-1;
        trackFunctionDefined= false;
        maxSize=0;
		max_ambiguous =0;
		max_unambiguous =0;
        complementSet=false;
    }
	
    //!Get the letter/word that has a given digitized value
    //! \param iter integer value of digital symbol
    //! \return std::string The string value in the undigitized sequence that is associated with the integer digitized value;
    std::string track::getAlpha(int iter){
		if (iter<=max_unambiguous){
			return alphabet[iter];
		}
		else if (iter <= max_ambiguous){
			return getAmbiguousCharacter(iter);
		}
		else{
			return "";
		}
    }
    
    //FIXME: Have it check before adding value
    //! Add a letter/word symbol to the track
    //! \param character Word or symbol used in undigitized sequence
    bool track::addAlphabetChar(std::string& character){
        
        if (alphabet.size() >= 255){
            std::cerr << "Alphabet limit reached.   Unable to add additional characters to the track:\t" << character << std::endl;
            return false;
        }
        
        
        if (character.size()>maxSize){maxSize=character.size();};
        
		alphabet.push_back(character);
        
		size_t index = alphabet.size()-1;
        
		symbolIndices[character]=index;
		
		max_unambiguous = index;
		unambiguous.push_back(index);
        
		return true;
    }
    
    bool track::addAlphabetChar(const char *character){
        std::string string_character(character);
        return addAlphabetChar(string_character);
    }
    
	
	//FIXME:  Need to fix the code below and test
	//Complements added by Ken
    bool track::addAlphabetChar(std::vector<std::string>& characters, std::vector<std::string>& complement_char){
        
        if (characters.size() != complement_char.size()){
            //Error Message
        }
        
        
        for(size_t i=0;i<characters.size();i++){
            addAlphabetChar(characters[i]);
        }
        
        for (size_t j=0;j<characters.size();j++){
            addComplement(characters[j], complement_char[j]);
        }
        
        complementSet = true;
        
        return true;
    }
    
    void track::addComplement(std::string& character, std::string& complement){
        int index = symbolIndex(character);
        int comp  = symbolIndex(complement);
        complementAlphabet[index]=comp;
        complementSet = true;
        
        return;
    }
    
    void track::addComplement(const char *character, const char *complement) {
        std::string string_character(character);
        std::string string_complement(complement);
        addComplement(string_character, string_complement);
        
        return;
    }
    
    //!Add an ambiguous character/word definition to the track
    //! \param ambChar  word/symbol fore the ambiguous character
    void track::addAmbiguous(std::string& ambChar, std::vector<std::string>& defs){
		if (defaultAmbiguous==-1){
			defaultAmbiguous = max_unambiguous+1;
		}
        ambigCharacter amb(this,ambChar,defs);
        ambiguousSymbols.push_back(amb);
        
		int index = (int) symbolIndices.size();  //Get the index position and new digital reference value
		
		if (index >= 255){
			std::cerr << "Maximum number of discrete symbols reached at 255\n";
			exit(2);
		}
		
        symbolIndices[ambChar]=index;
		max_ambiguous = index;
        return;
    }
	
	
    //! Get symbol assigned integer value
    //! \param symbol word/letter/symbol that we want to get it's assigned integer value
    uint8_t track::symbolIndex(std::string& symbol){
        if (symbolIndices.count(symbol)==0){  //If isn't found in the hash
            if (ambiguous){ //Return default character if ambiguous is set
                return defaultAmbiguous;
            }
            else{
                std::cerr << "Encountered an ambiguous character in the sequence.  No ambiguous characters are allowed because they weren't set in the model.   To allow ambiguous characters, please add an \" Ambiguous Character Definition\" to the model" << std::endl;
                exit(1);
                //return -1;
                //errorInfo(sCantHandleAmbiguousCharacter, "Ambiguous character handling is off but ambiguous character was encountered");
            }
        }
        else{
            return symbolIndices[symbol];
        }
    }
    
    
    //FIXME: Change return value so only returns true if parse is OK
    //! Parse a string representation of track to define a tracks parameters
    //! \param txt Line from model that describes a track
    //! \return true if the track was parsed properly
    bool track::parse(std::string& txt){
        stringList lst;
        stringList tag = extractTag(txt);
        size_t index;
        
        lst.fromNext(txt);
        setName(lst[0]);
        setDescription(lst.getComment());
        
        if (lst[1].compare("REAL_NUMBER")==0){
            setAlphaType(REAL);
            if (tag.size()>0){
                if (tag.contains("FUNCTION")){
                    index=tag.indexOf("FUNCTION");
                    trackFunction=tag[index+1];
                }
                else{
                    std::cerr << "Real number track function tag must contain FUNCTION: and USE: . Please check the formatting of your tag.  Here is the tag as parsed: " <<  tag.stringify() << std::endl;
                    return false;
                    
                }
                
                if (tag.contains("USE")){
                    index=tag.indexOf("USE");
                    trackToUse=tag[index+1];
                }
                else{
                    std::cerr << "Real number track tag must contain FUNCTION: and USE: . Please check the formatting of your tag.  Here is the tag as parsed: " <<  tag.stringify() << std::endl;
                    return false;
                }
                trackFunctionDefined=true;
            }
        }
        else{
            setAlphaType(ALPHA_NUM);
			
            for(int i=1;i<lst.size();i++){
                if (!addAlphabetChar(lst[i])){
                    std::cerr << "Track import failed, because number of symbols exceeded 255. Alternatively, you can create a real number track for different emissions" << std::endl;
                    return false;
                }
                
            }
        }
        return true;
    }
	
    //! Print the string representation of the track to stdout
    void track::print(){
        std::cout << stringify() << std::endl;
    }
    
    //! Get the string representation of the track
    //! \return std::string Definition of the track as in model
    std::string track::stringify(){
        std::string output;
        output+=name + ":\t";
        
        if (alpha_type == ALPHA_NUM){
            output+=join(alphabet,',');
        }
        else{
            output+="REAL_NUMBER";
            if (trackFunctionDefined){
                output+="\t[FUNCTION:\t" + trackFunction;
                output+="\tUSE:\t" + trackToUse + "]";
            }
        }
        output+="\n";
        
        return output;
    }
    
    //!Get the string representation of the ambiguous character definitions as in model file
    //! \return std::string
    std::string track::stringifyAmbig(){
        std::string output;
        output+=name + ":\t";
        for (size_t i = max_unambiguous+1; i <= max_ambiguous; i++){
            if (i>max_unambiguous+1){ output+= ",";}
            
            output+=getAmbiguousCharacter(i);
            output+="{";
            
            std::vector<size_t>& regChar = getAmbiguousSet(i);
            for(size_t k = 0; k<regChar.size();k++){
                if (k>0){output+=",";}
                output+=alphabet[regChar[k]];
            }
            output+="}";
        }
        return output;
    }
    
    
    
    std::string track::convertIndexToWord(size_t wordIndex, size_t order){
        std::string output="";
        
        if (order == 0){
            return "";
        }
        
        size_t currentOrder = order;
        
        while (currentOrder>1){
            double dreg=POWER[currentOrder-1][alphabet.size()-1];
            size_t temp = floor ((double) wordIndex / dreg);
            output+=alphabet[temp];
            if (maxSize!=1){
				output += ",";
            }
			
            wordIndex-=temp*dreg;
            currentOrder--;
        }
        
        output+=alphabet[wordIndex];
        
        return output;
    }
	
	
	void track::convertIndexToDigital(size_t wordIndex, size_t order, uint8_t word[]){
        if (order == 0){
			word[0]=wordIndex;
			return;
        }
		
		std::cout << alphabet.size() << std::endl;
        
        size_t currentOrder = order;
        
        while (currentOrder>1){
            double dreg=POWER[currentOrder-1][alphabet.size()-1];
            size_t temp = floor ((double) wordIndex / dreg);
			word[currentOrder-1] = temp;
            wordIndex-=temp*dreg;
            currentOrder--;
        }
        
        word[0] = wordIndex;
        return;
    }
    
    
    //FIXME: Change return value so only returns true if parse is OK
    //! Parse the ambiguous character definitions from model file
    //! \param txt String representation of ambiguous character definition as in model file
    //! \return true if the ambiguous characters were properly parsed
    bool track::parseAmbiguous(std::string& txt){
        std::vector<std::pair<std::string,std::vector<std::string> > > temp;
        
        _splitAmbiguousList(temp, txt);
        setAmbiguous();
        for (size_t i=0;i<temp.size();i++){
            addAmbiguous(temp[i].first, temp[i].second);
        }
        
        return true;
    }
    
    void track::_splitAmbiguousList(std::vector<std::pair<std::string,std::vector<std::string> > >& results, const std::string& text){
        
        size_t opening;
        size_t closing;
        size_t start=0;
        
        opening=text.find_first_of('{');
        while(opening!=std::string::npos){
            std::pair<std::string,std::vector<std::string> > amb;
            amb.first=text.substr(start,opening-start);
            
            closing=text.find_first_of('}',opening);
            if (closing!=std::string::npos){
                std::string tempString=text.substr(opening+1,closing-opening-1);
                split_line(amb.second, tempString);
            }
            start=text.find_first_not_of(',',closing+1);
            results.push_back(amb);
            opening=text.find_first_of('{',closing);
        }
        
        return;
    }
    
    
    //! Get the string representation of the ambigous character defined by integer value
    //! If ambiguous character isn't defined, return value is "*"
    //! \param val Integer value representing the ambiguous character
    std::string track::getAmbiguousCharacter(int val){
        if (getAmbiguousSize()==0){
            return "*";
        }
        
        return ambiguousSymbols[(val-max_unambiguous)-1].getSymbol();
    }
    
    
    //! Add track to tracks container
    //! \param tk Pointer to track to be added
    void tracks::push_back(track* tk){
        std::string& name= tk->name;
        
        if (!index.count(name)){
            index[tk->getName()]=trks.size();
            trks.push_back(tk);
        }
        else{
            std::cerr << "Track with name: " << name << " already exists.  Cannot add tracks with the same name\n";
            exit(1);
        }
        
    }
    
    //!Get iterator index of track by tracks name
    //! \param name Name of the track
    //! \return size_t Iterator to track within the tracks
    //! \return -1 if track doesn't exist in tracks
    size_t tracks::indexOf(const std::string& name){
        if (index.count(name)){
            return index[name];
        }
        else{
            return SIZE_MAX;
        }
    }
    
    
    //!Get pointer to track from the track name
    //! \param name Name of the track
    //! \return pointer to track if it exists, NULL otherwise
    track* tracks::getTrack(const std::string& name){
        if (index.count(name)){
            return trks[index[name]];
        }
        else{
            return NULL;
        }
    }
    
    bool tracks::isTrackDefined(const std::string& name){
        if (index.count(name)){
            return true;
        }
        
        return false;
    }
    
    //!Print the each track in tracks to stdout
    void tracks::print(){
        std::cout << stringify() << std::endl;
    }
    
    //! Get string representation of each track in tracks
    //! \return std::string String representation of tracks as in model file
    std::string tracks::stringify(){
        std::string trackString;
        std::string ambigString;
        std::string lnSep(50,'=');
        trackString+="TRACK SYMBOL DEFINITIONS\n" + lnSep + "\n";
        //ambigString+="AMBIGUOUS SYMBOL DEFINITIONS\n" + lnSep + "\n";
        
        for(size_t i=0;i<trks.size();i++){
            trackString+=trks[i]->stringify();
            if (trks[i]->isAmbiguousSet()){
                ambigString+=trks[i]->stringifyAmbig();
            }
        }
        if (!ambigString.empty()){
            ambigString="AMBIGUOUS SYMBOL DEFINITIONS\n" + lnSep + "\n"+ ambigString;
            trackString+="\n" + ambigString + "\n";
        }
        else{
            trackString+="\n";
        }
        
        return trackString;
    }
    
    
    //! Get the complement alphabet character digitized value given a value
    //! \param val Value of character to get complement of
    //! \return int value of complement
    uint8_t track::getComplementIndex(uint8_t val){
        if (complementSet){
            if (complementAlphabet.count(val)){
                return complementAlphabet[val];
            }
            else{
                std::cerr<< "Complement of " << val << " is not set in track\n";
                return -1;
            }
        }
        else{
            std::cerr << "No complements are set in the track\n";
            exit(1);
        }
    }
    
    
    
    //! Get the complement alphabet character digitized value given the string
    //! \param character String of alphanumerical symbol
    //! \return int Defined complement string symbol of symbol
    uint8_t track::getComplementIndex(std::string& character){
        if (complementSet){
            if (symbolIndices.count(character)){
                int characterIndex = symbolIndex(character);
                return complementAlphabet[characterIndex];
            }
            else{
                std::cerr<< "Complement of " << character << " is not set in track\n";
                return -1;
            }
        }
        else{
            std::cerr << "No complements are set in the track\n";
            exit(1);
        }
    }
    
    
    //! Get the complement alphanumerical string of a given integer value;
    //! \param value Integer value of a symbol
    //! \return std::string Defined complement string symbol of symbol with digitized value
    std::string track::getComplementSymbol(uint8_t value){
        if (complementSet){
            if (complementAlphabet.count(value)){
                int complement_value = complementAlphabet[value];
                return getAlpha(complement_value);
            }
            else{
                std::cerr<< "Complement of " << value << " is not set in track\n";
                return " ";
            }
        }
        else{
            std::cerr << "No complements are set in the track\n";
            exit(1);
        }
		
    }
    
    
    //! Get the compelment alphabet character digitized value given the string
    //! \param character String of alphanumerical symbol
    //! \return std::string Defined complement string symbol
    std::string track::getComplementSymbol(std::string& character){
        if (complementSet){
            if (symbolIndices.count(character)){
                int characterIndex = symbolIndex(character);
                int complement_value = complementAlphabet[characterIndex];
                return getAlpha(complement_value);
            }
            else{
                std::cerr<< "Complement of " << character << " is not set in track\n";
                return " ";
            }
        }
        else{
            std::cerr << "No complements are set in the track\n";
            exit(1);
        }
    }
    
	
}
