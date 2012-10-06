//
//  stateTemplate.cpp
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

#include "modelTemplate.h"
namespace StochHMM{
    
    stencil::stencil(){};
    
    stencil::stencil(std::string& txt){
        parse(txt);
    }
    
    /*!Takes a template string and finds all user defined variables (userParameters).  Adds these to a set of strings */
    
    bool stencil::parse(std::string& txt){
        
        states = txt;
        size_t opening=txt.find("<<");
        size_t closing=0;
        //!Find all userParameters Marks
        while (opening!=std::string::npos){
            closing=txt.find(">>",opening+1);
            std::string tagName=txt.substr(opening+2,closing-opening-2);
            UserParameterKeys.insert(tagName);
            opening=txt.find("<<",closing+1);
        }
        return true;
    }
    
    
    
    
    //!Returns a string with template filled out
    /*!Takes the template map and fills it out according to the specific Template you are looking at and the ID number that you want to use
     \param Template Name of template you want to return
     \param ID The number iterated in this return
     \param map_template map string of strings filled out in other subroutines
     */
    std::string stencil::getTemplate (std::string& Template, std::string& ID, std::map<std::string,std::string>& UserValues){
        
        std::string FilledOutModel = states;
        size_t beginposition=FilledOutModel.find("<<");
        size_t endposition=0;
        //!Find all userParameters and replace
        while (beginposition!=std::string::npos){
            endposition=FilledOutModel.find(">>",beginposition+1);
            /*!Check the coordinates and replace with proper coordinats
             The name should be only the tagName without brackets
             But the start should be the position of first bracket and stop
             should be the position of the last bracket */
            std::string tagName=FilledOutModel.substr(beginposition+2,endposition-beginposition-2);
            std::string ValueToFill = UserValues[tagName];
            FilledOutModel.replace(beginposition,(endposition+2-beginposition),ValueToFill);
            beginposition=FilledOutModel.find("<<",endposition+1);
        }
        
        beginposition=FilledOutModel.find("((");
        endposition=0;
        while (beginposition!=std::string::npos){
            endposition=FilledOutModel.find("))",beginposition+1);
            /*! */
            std::string tagName=FilledOutModel.substr(beginposition+2,endposition-beginposition-2);
            std::string ValueToFill = tagName + ID;
            FilledOutModel.replace(beginposition,(endposition+2-beginposition),ValueToFill);
            beginposition=FilledOutModel.find("((",endposition+1);
        }
        return FilledOutModel;
    }
    
    templates::templates(){}
    
    templates::templates(std::string& template_string){
        templates::parse(template_string);
    }
    
    templates::~templates(){
        std::map<std::string,stencil*>::iterator it;
        for(it=temps.begin();it!=temps.end();it++){
            std::string key = (*it).first;
            delete temps[key];
        }
    }
    
    //!Parses the model templates string
    /*!Takes the string passed to it and parses out each template to then be added to the Templates structure using stensil on each one */
    bool templates::parse(std::string& model_temp){
        
        std::string copyofpassedstring = model_temp;
        size_t model_end_position=model_temp.find("TEMPLATE: ",10);
        size_t model_start_position=0;
        std::string modelstring = "";
        std::string stenName;
        size_t stenNameStart = 10;
        size_t stenNameEnd;
        while (model_start_position!=std::string::npos){
            stencil* stencil_to_push = new(std::nothrow) stencil();
            if (stencil_to_push==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            model_end_position=copyofpassedstring.find("TEMPLATE: ",(model_start_position + 10)) - 1;        
            
            modelstring=model_temp.substr(model_start_position,model_end_position - model_start_position);
            
            stenNameEnd=modelstring.find("\n");
            
            stenName=model_temp.substr(stenNameStart,stenNameEnd-stenNameStart);
            
            stencil_to_push->parse(modelstring);
            //Need to correctly push the stencil to the map
            temps[stenName] = stencil_to_push;
            model_start_position=model_end_position + 1;
        }
        return true;
    }
    
    
    //!Adds more templates when passed a template string
    void templates::addTemplate(std::string& template_string){
        templates::parse(template_string);
    }
    
    
    
    //!Returns a string with template filled out
    /*!Takes the template map and fills it out according to the specific Template you are looking at and the ID number that you want to use
     \param Template Name of template you want to return
     \param ID The number iterated in this return
     */
    std::string templates::getTemplate (std::string& Template, std::string& ID,std::map<std::string,std::string>& UserVals){
        
        std::string FilledOutModel = temps[Template]->getTemplate(Template,ID,UserVals);
        return FilledOutModel;
    }
    
    
    
    
    
//    //!Returns all stencils in templates as a single string
//    std::string templates::returnAllFilledTemplates(std::string& ID,std::map<std::string,std::string>& UserVals){
//        std::map<std::string,stencil*>::iterator it;
//        std::string FilledOutModel = "";
//        for(it=temps.begin();it!=temps.end();it++){
//            std::string key = (*it).first;
//            FilledOutModel = FilledOutModel + returnFilledTemplate (key,ID,UserVals);
//        }
//        return FilledOutModel;
//    }
}

