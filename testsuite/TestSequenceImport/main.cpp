//
//  main.cpp
//  TestSequenceImport
//
//  Created by Paul Lott on 4/16/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#include <iostream>
#include <string>
#include "hmm.h"
#include "seqTracks.h"
using namespace StochHMM;



int main(int argc, const char * argv[])
{
    
    //Single Tracks
    seqTracks fastaTrks;
    seqTracks fastqTrks;
    
    model test;
    std::string modelFile("TestSeqs/SingleTrack_Dice/DiceNew.hmm");

    //Test Single Track
    std::string fastaFiles[]={  "Dice1.fa", "Dice2.fa", "Dice3.fa", "Dice4.fa", "Dice5.fa", "Dice6.fa", "Dice7.fa",
                                "Dice8.fa", "Dice9.fa", "Dice10.fa", "Dice11.fa", "Dice12.fa"};
    
    std::string fastqFiles[]={"Dice1.fastq","Dice2.fastq","Dice3.fastq","Dice4.fastq","Dice5.fastq","Dice6.fastq",
        "Dice7.fastq", "Dice8.fastq", "Dice9.fastq", "Dice10.fastq", "Dice11.fastq", "Dice12.fastq",};


    if (test.import(modelFile)){
//        for (int i=0;i<6;i++){
//            std::string fileTemp = "TestSeqs/Dice/" + files[i];
//            trks.loadSeqs(test, fileTemp, FASTA, NULL);
//            for (int j=0;j<trks.size();j++){
//                seqJob* temp = trks.getJob();
//                temp->printSeq();
//            }
//        }
    }
    
//Test Single Track, Multiple Files
    std::vector<std::string> fastaFileVec(fastaFiles, fastaFiles+12);
    std::vector<std::string> fastqFileVec(fastqFiles, fastqFiles+12);
    
//    for (int i=0;i<12;i++){
//        fastaFileVec[i] = "TestSeqs/SingleTrack_Dice/" + fastaFileVec[i];
//        fastqFileVec[i] = "TestSeqs/SingleTrack_Dice/" + fastqFileVec[i];
//    }
//    
//    fastaTrks.loadSeqs(test, fastaFileVec, FASTA, SINGLE_TRACK);
    seqJob* temp;
//    while(fastaTrks.remainingSeqs()!=0){
//        temp = fastaTrks.getJob();
//        if (temp == NULL){
//            continue;
//        }
//        temp->printSeq();
//    }
//    
//    fastqTrks.loadSeqs(test, fastqFileVec, FASTQ, SINGLE_TRACK);
//    while(fastqTrks.remainingSeqs()!=0){
//        temp = fastqTrks.getJob();
//        if (temp == NULL){
//            continue;
//        }
//        temp->printSeq();
//    }
    
//Multiple Track, Multiple Files
    modelFile ="TestSeqs/MultiTrack_Dice/DiceNew.hmm";
    model multiple;
    multiple.import(modelFile);
    
    fastaFileVec.clear();
    fastaFileVec.insert(fastaFileVec.begin(), fastaFiles, fastaFiles+12);
    
    fastqFileVec.clear();
    fastqFileVec.insert(fastqFileVec.begin(), fastqFiles, fastqFiles+12);
    
    
    for (int i=0;i<12;i++){
        fastaFileVec[i] = "TestSeqs/MultiTrack_Dice/" + fastaFileVec[i];
        fastqFileVec[i] = "TestSeqs/MultiTrack_Dice/" + fastqFileVec[i];
    }
    
//    fastaTrks.loadSeqs(multiple, fastaFileVec, FASTA, MULTI_TRACK);
//    while(fastaTrks.remainingSeqs()!=0){
//        temp = fastaTrks.getJob();
//        if (temp == NULL){
//            continue;
//        }
//        temp->printSeq();
//    }
    
    fastqTrks.loadSeqs(multiple, fastqFileVec, FASTQ, MULTI_TRACK);
    
    while(fastqTrks.remainingSeqs()!=0){
        temp = fastqTrks.getJob();
        if (temp == NULL){
            continue;
        }
        temp->printSeq();
        
        temp->printFilenames();
    }

    
    return 0;
}