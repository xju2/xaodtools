#ifndef __MYXAODTOOLS_CANDIDATE_H__
#define __MYXAODTOOLS_CANDIDATE_H__

#include <TLorentzVector.h>

class Candidate{
private:
    int index;
    TLorentzVector P;
    float charge;

public:
    Candidate(int index,TLorentzVector pp){
        this->index = index;
        this->P = pp;
    }
    void setCharge(float charge){
        this->charge = charge;
    }
    float getCharge(){
        return this->charge;
    }

    TLorentzVector getFourVector(){
    	return this->P;
    }
};

#endif
