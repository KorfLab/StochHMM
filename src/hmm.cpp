//hmm.cpp
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
#include "hmm.h"
namespace StochHMM{
	
    
    //! Import multiple models
    //! \param modelFile Path to multiple model file
    //! \param funcs Pointer to state Functions
    void models::importModels(std::string& modelFile, StateFuncs* funcs){
        std::ifstream MOD;
        MOD.open(modelFile.c_str());
        if (MOD.fail()){
            std::string info = "Model file not found: " + modelFile;
            std::cerr << info << std::endl;
            exit(1);
        }
        
        //Input is assigned to "input" & get first line
        std::string input;
        
        //Check file header for correct type
        getline(MOD, input, '\n');
#ifdef DEBUG_MODEL
        std::cout << input <<std::endl;
#endif
        if (input.compare("#STOCHHMM MODEL FILE")==0){
            MOD.close();
            model* temp= new(std::nothrow) model;
            
            if (temp==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            hmms.push_back(temp);
            hmms[0]->import(modelFile,funcs);
        }
        else if (input.compare("#STOCHHMM MODELS")==0){
            std::vector<std::string> filenames;
            
            while (getline(MOD,input,'\n')){
#ifdef DEBUG_MODEL
                sd::cout << input <<endl;
#endif
                //Need to check input to make sure its valid;
                //cout << input <<endl;
                if (input.compare("")==0){
                    break;
                }
                else{
                    filenames.push_back(input);
                }
            }
            MOD.close();
            
            for(int i=0;i<filenames.size();i++){
                model* temp=new(std::nothrow) model;
                
                if (temp==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                hmms.push_back(temp);
            }
            
            for(int i=0;i<filenames.size();i++){
                hmms[i]->import(filenames[i],funcs);
                hmms[i]->print();
            }
        }
        else{
            std::cerr << "Header for " << modelFile <<  "does not indicate that it is a StochHMM model file.\n";
            std::cerr << "Header should be:  #STOCHHMM MODEL FILE \n" << "Header given: " << input <<std::endl;
        }
        return;
    }
	
	
	// TODO:  Change to getModelByAttrib();
	//
	////----------------------------------------------------------------------------//
	//// Description: getGCModel(float gc)
	//// Returns a pointer to the model that has the closest GC range(calculated from
	//// midpoint to the sequence's gc percentage
	////
	////----------------------------------------------------------------------------//
	//model* models::getGCModel(float gc){
	//
	//    model *model=NULL;
	//    gc*=100.0f;  //Get percentage
	//    float min=100.0f;
	//    size_t minIterator=0;
	//
	//    for(size_t iter=0; iter<models.size(); iter++){
	//        //get distance from midpoint to sequences GC%
	//        float start=models[iter]->getLowerRange();
	//        float stop=models[iter]->getUpperRange();
	//        float midpoint=(stop-start)/2.0f;
	//        midpoint+=start;
	//
	//        float distance =fabs(midpoint - gc);
	//
	//        if (distance<min){  //if it is a minimum keep track of it
	//            minIterator=iter;
	//            min=distance;
	//        }
	//    }
	//
	//    if (min<100.0f){   //If we've seen a min then return it else return NULL
	//        model=models[minIterator];
	//    }
	//
	//    return model;   //return pointer to HMM
	//}
	
    //!Get pointer to model at index position
    //! \param iter Index iterator for position
    //! \return pointer to model
    model* models::getModel(size_t iter) {
        if (iter<this->size()) {
            return hmms[iter];
        }
        else{
            return NULL;
        }
    }
	
    
    
    //!Create a model
    model::model(){
        ending=new(std::nothrow) state();
        
        if (ending==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        finalized=false;
        basicModel=true;
        initial=NULL;
        scaling=NULL;
        templatedStates=NULL;
        attribTwo=false;
        range[0]=-INFINITY;
        range[1]=-INFINITY;
        
        return;
    }
    
	//    //TODO: more model constructors based on additional templates and weights to pass to import/parse
	//    //!Create a model file
	//    //! \param modelFile Filename of model file
	//    //! \param funcs  Pointer to State functions defined by programmer
	//    model::model(std::string& modelFile, StateFuncs* funcs){
	//        ending=new state();
	//        finalized=false;
	//        basicModel=true;
	//        initial=NULL;
	//        scaling=NULL;
	//        templatedStates=NULL;
	//        attribTwo=false;
	//        import(modelFile,funcs);
	//        return;
	//    }
    
    //TODO: multiple import function for templates and weights...
    //!Import the model file and parse it
    //! \param modelFile Model filename
    //! \param funcs  Pointer to State functions defined by programmer
    bool model::import(std::string& modelFile, StateFuncs* funcs){
        std::string modelString=slurpFile(modelFile);
        return parse(modelString,funcs,NULL,NULL);
    }
    
    
    bool model::import(std::string& modelFile, StateFuncs* funcs, templates* tmpls, weights* scl){
        std::string modelString=slurpFile(modelFile);
        return parse(modelString, funcs, tmpls, scl);
    }
    
    bool model::import(std::string& modelFile){
        std::string modelString=slurpFile(modelFile);
        return parse(modelString,NULL,NULL,NULL);
    }
    
	
    //!Parses text model file
    //!Splits the model into sections that are then parsed by the individiual classes
    //!parse() functions.
    bool model::parse(const std::string& model, StateFuncs* funcs, templates* tmpls, weights* scl ){
        
        templatedStates=tmpls;
        scaling=scl;
		
        //std::cout << model <<std::endl;
        
        size_t header = model.find("MODEL INFORMATION");
        size_t track  = model.find("TRACK SYMBOL DEFINITIONS");
        size_t ambig  = model.find("AMBIGUOUS SYMBOL DEFINITIONS");
        size_t templ  = model.find("TEMPLATED STATES");
        size_t scale  = model.find("SCALING VALUES");
        size_t st = model.find("STATE DEFINITIONS");
        size_t blank;
        size_t nlChar;
        
        //Parse Model Informaton (Optional)
        if (header!=std::string::npos){
            
            blank=model.find("\n\n",header); //Get coordinates of splitting Model Information
            
            size_t nlCharEq = model.rfind("####\n",blank);
            size_t nlCharNum= model.rfind("====\n",blank);
            //Check for optional dividing line
            if (nlCharEq!=std::string::npos){
                nlChar=nlCharEq+5;
            }
            else if (nlCharNum!=std::string::npos){
                nlChar=nlCharNum+5;
            }
            else{  //No divider line
                nlChar=model.find("\n",header);
                nlChar++;
            }
            
            std::string head = model.substr(nlChar,blank-nlChar);
            if (!_parseHeader(head)){
                return false;
            }
        }
        
        
        //Parse Tracks (Required)
        if (track!=std::string::npos){
            
            blank=model.find("\n\n",track);
            
            
            size_t nlCharEq = model.rfind("####\n",blank);
            size_t nlCharNum= model.rfind("====\n",blank);
            //Check for optional dividing line
            if (nlCharEq!=std::string::npos){
                nlChar=nlCharEq+5;
            }
            else if (nlCharNum!=std::string::npos){
                nlChar=nlCharNum+5;
            }
            else{  //No divider line
                nlChar=model.find("\n",track);
                nlChar++;
            }
            
            
            std::string trck (model.substr(nlChar,blank-nlChar));
            
            if (!_parseTracks(trck)){
                return false;
            }
        }
        else{
            std::cerr << "Required section: TRACK SYMBOL DEFINITIONS missing from the model" << std::endl;
            return false;
        }
		
		//        else{
		//            errorInfo(sMissingTrackDefinition, "Required section: TRACK SYMBOL DEFINITIONS missing from the model")
		//        }
        
        //Parse Ambiguous Characters (Optional)
        if (ambig!=std::string::npos){
            blank=model.find("\n\n",ambig);
            
            size_t nlCharEq = model.rfind("####\n",blank);
            size_t nlCharNum= model.rfind("====\n",blank);
            //Check for optional dividing line
            if (nlCharEq!=std::string::npos){
                nlChar=nlCharEq+5;
            }
            else if (nlCharNum!=std::string::npos){
                nlChar=nlCharNum+5;
            }
            else{  //No divider line
                nlChar=model.find("\n",ambig);
                nlChar++;
            }
            
            std::string amb(model.substr(nlChar,blank-nlChar));
            
            if (!_parseAmbiguous(amb)){
                return false;
            }
            //std::cout << ambig << "\t" << amb << std::endl;
        }
        
        //Parse Scaling Values (Optional)
        if (scale!=std::string::npos){
            blank=model.find("\n\n",scale);
            
            
            size_t nlCharEq = model.rfind("####\n",blank);
            size_t nlCharNum= model.rfind("====\n",blank);
            //Check for optional dividing line
            if (nlCharEq!=std::string::npos){
                nlChar=nlCharEq+5;
            }
            else if (nlCharNum!=std::string::npos){
                nlChar=nlCharNum+5;
            }
            else{  //No divider line
                nlChar=model.find("\n",scale);
                nlChar++;
            }
            
            
            
            std::string scaleTxt = model.substr(nlChar,blank-nlChar);
            if(!_parseScaling(scaleTxt)){
                return false;
            }
            //std::cout << scale << "\t" << scaleTxt << std::endl;
        }
		
		
        
        //Parse Templated Tracks (Optional)
        if (templ!=std::string::npos){
            blank=model.find("\n\n",templ);
            
            
            size_t nlCharEq = model.rfind("####\n",blank);
            size_t nlCharNum= model.rfind("====\n",blank);
            //Check for optional dividing line
            if (nlCharEq!=std::string::npos){
                nlChar=nlCharEq+5;
            }
            else if (nlCharNum!=std::string::npos){
                nlChar=nlCharNum+5;
            }
            else{  //No divider line
                nlChar=model.find("\n",templ);
                nlChar++;
            }
            
            std::string tempTxt = model.substr(nlChar,blank-nlChar);
            if (!_parseTemplates(tempTxt)){
                return false;
            }
            //std::cout << templ << "\t" << tempTxt << std::endl;
        }
        
        
        //Parse State Definitions (Required)
        if (st!=std::string::npos){
            size_t blankNum = model.find("####\n",st);
            size_t blankEq  = model.find("====\n",st);
            
            blank=model.find("####\n",st);
            
            
            if (blankEq!=std::string::npos){
                blank=blankEq+5;
            }
            else if (blankNum!=std::string::npos){
                blank=blankNum+5;
            }
            else{
                blank=model.find("\n",st);
                blank++;
            }
            
            
            nlChar=model.find("\n//END");
            if (nlChar==std::string::npos){
                nlChar=model.size()-1;
            }
            
            std::string stateTxt= model.substr(blank,nlChar-blank);
            if (!_parseStates(stateTxt,funcs)){
                return false;
            }
            
        }
        else{
            std::cerr << "Required sections <STATE DEFINITIONS> missing from the model" << std::endl;
            return false;
            //errorInfo(sMissingStateDefinition, "Required sections <STATE DEFINITIONS> missing from the model")
        }
        
        return true;
    }
	
	
    //Parse model header information
    bool model::_parseHeader(std::string& txt){
        stringList lst;
        size_t index;
        bool first(false);
        bool second(false);
        
        //std::string headers[] = {"MODEL_NAME", "MODEL_DESCRIPTION", "MODEL_CREATION_DATE","MODEL_CREATION_COMMAND", "MODEL_AUTHOR","MODEL_NUM_ATTRIB","MODEL_NUM_ATTRIB_UPPER","MODEL_NUM_ATTRIB_LOWER"};
        std::string headers[] = {"NAME", "DESCRIPTION", "CREATION_DATE","CREATION_COMMAND", "AUTHOR","NUM_ATTRIB","UPPER","LOWER"};
        
        std::string* head[]={&name, &desc, &date, &command,&author};
        
        lst.fromTxt(txt);
        
        
        for(int i=0;i<5;i++){
            if (lst.contains(headers[i])){
                index = lst.indexOf(headers[i]);
                if (index+1 < lst.size()){
                    index++;
                    (*head[i])=lst[index];
                    
                }
                else{
                    std::cerr << "Couldn't parse " << headers[i] << " from \"MODEL INFORMATION\" section." << std::endl;
                    return false;
                }
            }
        }
        
        //Set Numerical Attributes of Model
        if (lst.contains(headers[5])){
            index = lst.indexOf(headers[5]);
            
            if (index+1 < lst.size()){
                index++;
                
                
                double tempValue;
                if (!stringToDouble(lst[index], tempValue)){
                    std::cerr << "Numerical attribute couldn't be converted to numerical value: " << lst[index] << std::endl;
                    return false;
                }
                range[0]=tempValue;
                
            }
            else{
                std::cerr << "Couldn't parse " << headers[5] << " value from \"MODEL INFORMATION\" section." << std::endl;
                return false;
            }
        }
        else if (lst.contains(headers[6]) && lst.contains(headers[7])){
            index = lst.indexOf(headers[6]);
            
            if (index+1<lst.size()){
                index++;
                
                double tempValue;
                if (!stringToDouble(lst[index], tempValue)){
                    std::cerr << "Numerical attribute couldn't be converted to numerical value: " << lst[index] << std::endl;
                    return false;
                }
                range[1]=tempValue;
                second=true;
            }
            else{
                std::cerr << "Couldn't parse " << headers[6] << " value from \"MODEL INFORMATION\" section." << std::endl;
                return false;
            }
            
            
            index = lst.indexOf(headers[7]);
            
            if (index+1<lst.size()){
                index++;
                
                
                double tempValue;
                if (!stringToDouble(lst[index], tempValue)){
                    std::cerr << "Numerical attribute couldn't be converted to numerical value: " << lst[index] << std::endl;
                    return false;
                }
                range[0]=tempValue;
                first=true;
            }
            else{
                std::cerr << "Couldn't parse " << headers[6] << " value from \"MODEL INFORMATION\" section." << std::endl;
                return false;
            }
            
            if (first && second){
                attribTwo=true;
            }
            else{
                std::cerr << "Unable to parse both UPPER and LOWER" << std::endl;
                return false;
            }
        }
        else if (lst.contains(headers[6])){
            index = lst.indexOf(headers[6]);
            
            if (index+1<lst.size()){
                index++;
                
                double tempValue;
                if (!stringToDouble(lst[index], tempValue)){
                    std::cerr << "Numerical attribute couldn't be converted to numerical value: " << lst[index] << std::endl;
                    return false;
                }
                range[0]=tempValue;
                
                
            }
            else{
                std::cerr << "Couldn't parse " << headers[6] << " value from \"MODEL INFORMATION\" section." << std::endl;
                return false;
            }
        }
        else if (lst.contains(headers[7])){
            index = lst.indexOf(headers[6]);
            
            if (index+1<lst.size()){
                index++;
                
                double tempValue;
                if (!stringToDouble(lst[index], tempValue)){
                    std::cerr << "Numerical attribute couldn't be converted to numerical value: " << lst[index] << std::endl;
                    return false;
                }
                range[0]=tempValue;
            }
            else{
                std::cerr << "Couldn't parse " << headers[7] << " value from \"MODEL INFORMATION\" section." << std::endl;
                return false;
            }
        }
        
        return true;
    }
    
	
    bool model::_parseTracks(std::string& txt){
        stringList lst;
        lst.splitString(txt, "\n");
        for(size_t i=0;i<lst.size();i++){
            track* trk=new(std::nothrow) track();
            
            if (trk==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            if (trk->parse(lst[i])){
                trk->setIndex(trcks.size());
                trcks.push_back(trk);
            }
            else {
                std::string info = "Couldn't parse new track line.  Please check formatting of : " + lst[i];
                std::cerr << info << std::endl;
                exit(1);
                
                //errorInfo(sCantParseLine, info.c_str());
                
                //std::cerr << "Couldn't parse new track line.  Please check formatting of : " << lst[i] << std::endl;
                //exit(1);
            }
            
        }
        return true;
    }
	
    bool model::_parseAmbiguous(std::string& txt){
        stringList lst;
        lst.splitString(txt, "\n");
        for(size_t i=0;i<lst.size();i++){
            stringList ln;
            ln.splitString(lst[i],":\t ");
            track* trk = getTrack(ln[0]);
            if (trk!=NULL){
                if (!trk->parseAmbiguous(ln[1])){
                    std::cerr << "Couldn't parse the Ambiguous section for " << ln[0] << std::endl;
                    return false;
                }
            }
            else{
                std::string info = "Ambiguous Characters Section\nSupplied track name doesn't correspond to any previously parsed tracks.\nPlease check the formatting and names.\n Unfound Name: " + ln[0];
                //errorInfo(sAmbigousCharacterDoesntMatchTrack, info.c_str());
                
                std::cerr << info << std::endl;
                return false;
            }
            
        }
        //lst.print();
        return true;
    }
    
    
    bool model::_parseScaling(std::string& txt){
        
        if (scaling==NULL){
            scaling = new(std::nothrow) weights;
            
            if (scaling==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        
        stringList lst;
        lst.splitND(txt, "SCALE:");
        for(size_t i=0;i<lst.size();i++){
            weight* wt = new(std::nothrow) weight;
            
            if (wt==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            if (!wt->parse(lst[i])){
                return false;
            }
            
            
            scaling->addWeight(wt);
        }
        //lst.print();
        
        return true;
    }
    
    
    bool model::_parseTemplates(std::string& txt){
        if (!txt.empty()){
            templatedStates = new(std::nothrow) templates();
            
            if (templatedStates==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            if (!templatedStates->parse(txt)){
                std::cerr << "Unable to parse Templates for the model" << std::endl;
                return false;
            }
        }
        
        return true;
    }
    
    bool model::_parseStates(std::string& txt, StateFuncs* funcs){
        //1. split sections and identify any template sections
        //2. get state names list
        //3. create and parse states
        
        stringList stats;
        _splitStates(txt,stats);
        
        
        stringList NameList;
        _getOrderedStateNames(stats,NameList);
        
        for(size_t iter=0;iter<stats.size();iter++){
            state* st = new(std::nothrow) state;
            
            if (st==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            if (!st->parse(stats[iter],NameList,trcks,scaling, funcs)){
                delete st;
                return false;
            }
            
            
            if (st->getName() == "INIT"){
                initial=st;
                statesByName[st->getName()]=st;
            }
            else{
                states.push_back(st);
                statesByName[st->getName()]=st;
            }
        }
        
        
        //Post process states to create final state with only transitions from filled out.
        finalize();
        
        
        return true;
    }
    
    
    //!Add state to model
    //! \param st Pointer to state
    void model::addState(state* st){
        states.push_back(st);
        statesByName[st->getName()]=st;
        return;
    };
	
    
    //Split states into individual state strings
    bool model::_splitStates(std::string& txt ,stringList& sts){
        
        size_t start=0;
        size_t end=0;
        
        while(start!=std::string::npos){
            end=txt.find("STATE:",start+1);
            
            std::string st = txt.substr(start,end-start);
            
            clear_whitespace(st, "#");
            
            if (st.find("TEMPLATE:")!=std::string::npos){
                
                stringList tmpd;
                
                if (!_processTemplateState(st,tmpd)){ //! Get templated states using the defined parameters
                    std::cerr << "Unable to process template information" << std::endl;
                    return false;
                }
                
                for(size_t i=0;i<tmpd.size();i++){
                    sts.push_back(tmpd[i]);
                }
            }
            else{
                sts.push_back(st);
            }
			
            start=txt.find("STATE:",end);
            
        }
        return true;
    }
    
    bool model::_processTemplateState(std::string& txt, stringList& lst){
        if (templatedStates==NULL){
            std::cerr << "No template provided for State Definitions.  Please provide the template in the model or when calling HMM class" << std::endl;
            exit(1);
            
			//errorInfo(sMissingTemplate, "No template provided for State Definitions.  Please provide the template in the model or when calling HMM class");
        }
        std::string templateName;
        std::string templateIdentifier;
        std::map<std::string,std::string> parameters;
        
        lst.splitString(txt,"\n");
        
        //Process Template Name and Identifier
        clear_whitespace(lst[0], " \n");
        stringList nmid;
        nmid.splitString(lst[0],":\t");
        nmid.print();
        
        //Extracting TEMPLATE NAME
        if (nmid.contains("TEMPLATE")){
            size_t index=nmid.indexOf("TEMPLATE");
            index++;
            templateName=nmid[index];
        }
        else{
            std::cerr << "State Template definition doesn't contain \"TEMPLATE:\" keyword in first line of State definition of the template.   Please check formatting\n";
            return false;
        }
        
        //Extracting TEMPLATE IDENTIFIER
        if (nmid.contains("IDENTIFIER")){
            size_t index=nmid.indexOf("IDENTIFIER");
            index++;
            templateIdentifier=nmid[index];
        }
        else{
            std::cerr << "State Template definition doesn't contain \"IDENTIFIER:\" keyword in first line of State definition of the template.   Please check formatting\n";
            return false;
        }
		
        //Process Template Variables
        std::string lastTag;
        for(size_t i=1;i<lst.size();i++){
            std::string key;
            std::string value;
            getKeyValue(lst[i],key,value);
            if (key.empty()){
                parameters[lastTag]+="\n" + value;
            }
            else{
                lastTag=key;
                parameters[lastTag]=value;
            }
        }
        
        
        //Get filled out model template
        std::string filledTemplate = templatedStates->getTemplate(templateName, templateIdentifier, parameters);
        
        //Split filled out template into individual states
        stringList sts;
        
        if (!_splitStates(filledTemplate, sts)){ //Split the templated states
            std::cerr << "Unable to split templated states\n" << sts.stringify() << std::endl;
            return false;
        }
        
        return true;
    }
    
    
    bool model::_getOrderedStateNames(stringList& states, stringList& names){
        for(size_t i=0;i<states.size();i++){
            size_t nameHeader=states[i].find("NAME:");
            size_t nameLineEnding=states[i].find_first_of("\n",nameHeader);
            std::string name = states[i].substr(nameHeader+5,nameLineEnding-(nameHeader+5));
            clear_whitespace(name, " \t\n");
            if (names.contains(name)){
                std::cerr << "State with name of: " << name << " is defined twice in the model\n";
                return false;
            }
            else{
                names.push_back(name);
            }
        }
        return true;
    }
    
    
    //!Get pointer to state by state name
    //!\param txt String name of state
    //!\return pointer to state if it exists;
    //!\return NULL if state doesn't exist in model
    state* model::getState(const std::string& txt){
        if (statesByName.count(txt)){
            return statesByName[txt];
        }
        else {
            return NULL;
        }
    }
    
    //!Finalize the model before performing decoding
    //!Sets transitions, checks labels, Determines if model is basic or requires intermittent tracebacks
    void model::finalize(){
        if (!finalized){
            //Add States To Transitions
            std::set<std::string> labels;
            std::set<std::string> gff;
            std::set<std::string> name;
            
            //Create temporary hash of states for layout
            for(size_t i=0;i<states.size();i++){
                labels.insert(states[i]->getLabel());
                gff.insert(states[i]->getGFF());
                gff.insert(states[i]->getName());
            }
			
			for (size_t i=0; i < states.size() ; ++i){
				states[i]->setIter(i);
			}
			
            
            
            //Add states To and From transition
            
            for(size_t i=0;i<states.size();i++){
				states[i]->checkLabels(labels,gff,name);
                _addStateToFromTransition(states[i]);
                
            }
			
			_addStateToFromTransition(initial);
			
			
			//Now that we've seen all the states in the model
            //We need to fix the States transitions vector transi, so that the state
            //iterator correlates to the position within the vector
            for(size_t i=0;i<states.size();i++){
                states[i]->_finalizeTransitions(statesByName);
            }
            
            //Check to see if model is basic model
			//Meaning that the model doesn't call outside functions or perform
			//Tracebacks for explicit duration.
            for(size_t i=0;i<states.size();i++){
                std::vector<transition*>* transitions = states[i]->getTransitions();
                for(size_t trans=0;trans<transitions->size();trans++){
                    if ((*transitions)[trans]->FunctionDefined()){
                        basicModel=false;
                        break;
                    }
                    else if ((*transitions)[trans]->getTransitionType()!=STANDARD || (*transitions)[trans]->getTransitionType()!=LEXICAL){
                        basicModel=false;
                        break;
                    }
                }
                
                if (!basicModel){
                    break;
                }
            }
            
            finalized=true;
        }
        return;
        
    }
    
    void model::_addStateToFromTransition(state* st){
        std::vector<transition*>* trans;
        
        //Process Initial State
        trans = st->getTransitions();
        for(size_t i=0;i<trans->size();i++){
            state* temp;
            temp=this->getState((*trans)[i]->getName());
            if (temp){
                st->addToState(temp); //Also add the ptr to state vector::to
                (*trans)[i]->setState(temp);
                if (st!=initial){
                    temp->addFromState(st);
                }
            }
        }
        
        if (st->endi){
            ending->addFromState(st);
        }
    }
    
    
    //!Get pointer to track by track name
    //!\param name track name
    //!\return pointer to track
    //!\return NULL if track with given name doesn't exist in the model
    track* model::getTrack(const std::string& name){
        return trcks.getTrack(name);
    }
    
    
    //!Get string representation of model
    //!\return std::string
    std::string model::stringify(){
        std::string model;
        std::string lnSep(50,'=');
        model+="#STOCHHMM MODEL FILE\n\n";
        
        model+=_stringifyHeader();
        model+=_stringifyTracks();
        model+=_stringifyScaling();
        model+=_stringifyStates();
        return model;
    }
    
    std::string model::_stringifyHeader(){
        std::string headerString;
        std::string lnSep(50,'=');
        headerString+="MODEL INFORMATION\n" + lnSep + "\n";
        std::string headers[] = {"MODEL_NAME", "MODEL_DESCRIPTION", "MODEL_CREATION_DATE","MODEL_CREATION_COMMAND", "MODEL_AUTHOR","MODEL_NUM_ATTRIB","MODEL_NUM_ATTRIB_UPPER","MODEL_NUM_ATTRIB_LOWER"};
        std::string* head[]={&name, &desc, &date, &command,&author};
        
        for(int i=0;i<5;i++){
            if (!(*head[i]).empty()){
                headerString+=headers[i] + ":\t" + (*head[i]) + "\n";
            }
        }
        
        
        if (attribTwo){
            headerString+= headers[7] + ":\t" + double_to_string(range[1]) + "\n";
            headerString+= headers[8] + ":\t" + double_to_string(range[0]) + "\n";
        }
        headerString+="\n";
        
        return headerString;
    }
    
    std::string model::_stringifyTracks(){
        std::string trackString;
        trackString+=trcks.stringify();
        trackString+="\n";
        return trackString;
    }
    
    std::string model::_stringifyScaling(){
        std::string scaleString;
        if (scaling){
            scaleString+=scaling->stringify();
            scaleString+="\n";
        }
        
        return scaleString;
    }
    
    std::string model::_stringifyStates(){
        std::string stateString;
        std::string lnSep(50,'#');
        stateString+="STATE DEFINITIONS\n" + lnSep + "\n";
        
        if (initial){
            stateString+=initial->stringify();
        }
        
        stateString+= lnSep + "\n";
        
        for(size_t i=0;i<states.size();i++){
            if (states[i]){
                stateString+=states[i]->stringify();
                stateString+=lnSep + "\n";
            }
        }
        stateString+= "//END";
        return stateString;
    }
    
    
    //!Print the string representation to stdout
    void model::print(){
        std::cout << stringify() << std::endl;
        return;
    }
    
    
    //To improve we could group into those connected groups
    //Or we could allow user to define groups and clusters;
    //!Create graphvis file for the state of the model
    //!\param filepath  complete path to use for graphviz file
    //!\param q0 TRUE will draw the initial state
    void model::writeGraphViz(std::string filepath,bool q0){
        std::ofstream gv (filepath.c_str());
        if (gv.is_open()){
            std::string modelName=name;
            replaceChar(modelName, ' ', '_');
            
            gv << "digraph " << modelName << " {\n";
            gv << "\tratio=LR\n";
            gv << "\tsize=\"8,5\"\n";
            gv << "\tnode [shape = circle];\n";
            
//            std::bitset<STATE_MAX>* temp;
//            
//            if (q0){
//                temp=initial->getTo();
//                for(size_t temp_iter=0; temp_iter < temp->size(); temp_iter++){
//                    //gv << "\tINIT -> " << (*temp)[temp_iter]->getName() << ";\n";
//                }
//            }
//            
//            
//            for(size_t state_iter=0; state_iter < states.size(); state_iter++){
//                temp = states[state_iter]->getTo();
//                std::string stName=states[state_iter]->getName();
//                
//                if (stName.find_first_of("0123456789")==0){
//                    stName.insert(0,"_");
//                }
//                
//                
//                for(size_t temp_iter=0; temp_iter < temp->size(); temp_iter++){
//                    //std::string tempName=(*temp)[temp_iter]->getName();
//                    
//                    if (tempName.find_first_of("0123456789")==0){
//                        tempName.insert(0,"_");
//                    }
//                    
//                    gv << "\t" << stName << " -> " << tempName << ";\n";
//                }
//                
//                
//                if (q0 && states[state_iter]->endi){
//                    gv << "\t" << stName << " -> END;\n";
//                }
//            }
            
            gv << "}\n";
        }
        gv.close();
        
    }
    
    //!Get pointer to weight by name of weight
    //! \param name Name of weight
    //! \return pointer to weight
    weight* model::getScalingFactor(std::string& name){
        if (scaling==NULL){
            return NULL;
        }
        else{
            return (*scaling)[name];
        }
    }
    
    
    //!Get Distance from the model attrib values to user defined value
    double model::getDistanceToAttrib(double val){
        double mid=range[0];
        if (attribTwo){
            mid+=(range[1]-range[0])/2;
        }
        
        return abs(mid-val);
    }
    
    
    void print_vec (std::vector<std::vector<double> > &x){
        size_t sze=x.size();
        size_t internal_sze=x[0].size();
        
        for(int j=0;j<sze;j++){
            //cout << j << "\t";
            for(int i=0;i<internal_sze;i++){
                std::cout << x[j][i] << '\t' ;
            }
            std::cout << std::endl;
        }
    }
    
}
