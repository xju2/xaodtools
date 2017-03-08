#include "MyXAODTools/HZZ4lHelper.h"

#include "xAODTruth/TruthVertex.h"
#include "xAODBase/IParticleHelpers.h"

#include <memory>
#include <math.h>
#include <algorithm>

#include <TTree.h>
using namespace std;

HZZ4lHelper::HZZ4lHelper():
    kZMASS(91187.6),
    m_tree(NULL)
{
}

HZZ4lHelper::~HZZ4lHelper(){
}

bool HZZ4lHelper::Is_close2Z(vector<Candidate*>* lep_4vec, double& m12, int& idL1, int& idL2)
{
    double tmpZ1 = -10.0;
    for(int i=0; i < (int)lep_4vec->size(); i++){
        for(int j=i+1; j<(int)lep_4vec->size(); ++j){
            float ch1 = lep_4vec->at(i)->getCharge();
            float ch2 = lep_4vec->at(j)->getCharge();
            if(ch1 * ch2 == 1)continue;
            double dimass = (lep_4vec->at(i)->getFourVector()+
                    lep_4vec->at(j)->getFourVector()).M();
            if(fabs(dimass - kZMASS) < fabs(tmpZ1 - kZMASS)){
                    idL1 = i; 
                    idL2 = j;
                    tmpZ1 = dimass;
            }
        }
    }
    if(tmpZ1 > 0) {
        m12 = tmpZ1;
        return true;
    }else{ 
        return false;
    }

}

bool HZZ4lHelper::Is_4mu4e(vector<Candidate*>* lep_4vec )
{
    double tmpZ1 = 0.0;
    int idL1 = -1, idL2 = -1, idL3 = -1, idL4 = -1;
    if( Is_close2Z(lep_4vec, tmpZ1, idL1, idL2) ){
        unique_ptr<vector<Candidate*> > newlep(new vector<Candidate*>());
        for(int i=0; i < (int) lep_4vec->size(); i++){
            if(i==idL1 || i==idL2) continue;
            newlep->push_back(lep_4vec->at(i));
        }
        double tmpZ2 = 0.0;
        if( Is_close2Z(newlep.get(), tmpZ2, idL3, idL4) )
        {
            m_mZ1 = (float)tmpZ1; 
            m_mZ2 = (float)tmpZ2;
            TLorentzVector Higgs_tlv(lep_4vec->at(idL1)->getFourVector());
            Higgs_tlv += lep_4vec->at(idL2)->getFourVector();
            Higgs_tlv += newlep->at(idL3)->getFourVector();
            Higgs_tlv += newlep->at(idL4)->getFourVector();
            m_m4l = (float)Higgs_tlv.M();
            m_Hpt_ = (float) Higgs_tlv.Pt();
            m_Hphi_ = (float) Higgs_tlv.Phi();

            if(lep_4vec->at(idL1)->getCharge() > 0){
                m_Z1_lepplus_pt  = lep_4vec->at(idL1)->getFourVector().Pt();
                m_Z1_lepminus_pt = lep_4vec->at(idL2)->getFourVector().Pt();
                m_Z1_lepplus_eta  = lep_4vec->at(idL1)->getFourVector().Eta();
                m_Z1_lepminus_eta = lep_4vec->at(idL2)->getFourVector().Eta();
            }else{
                m_Z1_lepplus_pt  = lep_4vec->at(idL2)->getFourVector().Pt();
                m_Z1_lepminus_pt = lep_4vec->at(idL1)->getFourVector().Pt();
                m_Z1_lepplus_eta  = lep_4vec->at(idL2)->getFourVector().Eta();
                m_Z1_lepminus_eta = lep_4vec->at(idL1)->getFourVector().Eta();
            }
            if(newlep ->at(idL3)->getCharge() > 0){
                m_Z2_lepplus_pt  = newlep->at(idL3)->getFourVector().Pt();
                m_Z2_lepminus_pt = newlep->at(idL4)->getFourVector().Pt();
                m_Z2_lepplus_eta  = newlep->at(idL3)->getFourVector().Eta();
                m_Z2_lepminus_eta = newlep->at(idL4)->getFourVector().Eta();
            }else{
                m_Z2_lepplus_pt  = newlep->at(idL4)->getFourVector().Pt();
                m_Z2_lepminus_pt = newlep->at(idL3)->getFourVector().Pt();
                m_Z2_lepplus_eta  = newlep->at(idL4)->getFourVector().Eta();
                m_Z2_lepminus_eta = newlep->at(idL3)->getFourVector().Eta();
            }
            return true;
        }
    }
    return false;
}

bool HZZ4lHelper::Is_2e2mu(
        vector<Candidate*>* ele_4vec,
        vector<Candidate*>* muon_4vec, int& type)
{
    double tmpZ1 = 0, tmpZ2 = 0;
    int idL1 = -1, idL2 = -1, idL3 = -1, idL4 = -1;
    if( Is_close2Z(ele_4vec, tmpZ1, idL1, idL2) )
    {
        if( Is_close2Z(muon_4vec, tmpZ2, idL3, idL4) )
        {
            TLorentzVector Higgs_tlv(ele_4vec->at(idL1)->getFourVector());
            Higgs_tlv += ele_4vec->at(idL2)->getFourVector();
            Higgs_tlv += muon_4vec->at(idL3)->getFourVector();
            Higgs_tlv += muon_4vec->at(idL4)->getFourVector();
            m_m4l =  (float) Higgs_tlv.M();
            m_Hpt_ = (float) Higgs_tlv.Pt();
            m_Hphi_ = (float) Higgs_tlv.Phi();

            if(fabs(tmpZ1 - kZMASS) < fabs(tmpZ2 - kZMASS)){
                type = 2; //2e2mu
                m_mZ1 = tmpZ1; 
                m_mZ2 = tmpZ2;
                if(ele_4vec->at(idL1)->getCharge() > 0){
                    m_Z1_lepplus_pt   = ele_4vec->at(idL1)->getFourVector().Pt();
                    m_Z1_lepminus_pt  = ele_4vec->at(idL2)->getFourVector().Pt();
                    m_Z1_lepplus_eta   = ele_4vec->at(idL1)->getFourVector().Eta();
                    m_Z1_lepminus_eta  = ele_4vec->at(idL2)->getFourVector().Eta();
                }else{
                    m_Z1_lepplus_pt   = ele_4vec->at(idL2)->getFourVector().Pt();
                    m_Z1_lepminus_pt  = ele_4vec->at(idL1)->getFourVector().Pt();
                    m_Z1_lepplus_eta   = ele_4vec->at(idL2)->getFourVector().Eta();
                    m_Z1_lepminus_eta  = ele_4vec->at(idL1)->getFourVector().Eta();
                }
                if(muon_4vec->at(idL3)->getCharge() > 0){
                    m_Z2_lepplus_pt   = muon_4vec->at(idL3)->getFourVector().Pt();
                    m_Z2_lepminus_pt  = muon_4vec->at(idL4)->getFourVector().Pt();
                    m_Z2_lepplus_eta   = muon_4vec->at(idL3)->getFourVector().Eta();
                    m_Z2_lepminus_eta  = muon_4vec->at(idL4)->getFourVector().Eta();
                }else{
                    m_Z2_lepplus_pt   = muon_4vec->at(idL4)->getFourVector().Pt();
                    m_Z2_lepminus_pt  = muon_4vec->at(idL3)->getFourVector().Pt();
                    m_Z2_lepplus_eta   = muon_4vec->at(idL4)->getFourVector().Eta();
                    m_Z2_lepminus_eta  = muon_4vec->at(idL3)->getFourVector().Eta();
                }
            }else{
                type = 3; //2mu2e
                m_mZ1 = tmpZ2; m_mZ2 = tmpZ1;
                if(ele_4vec->at(idL1)->getCharge() > 0){
                    m_Z2_lepplus_pt   = ele_4vec->at(idL1)->getFourVector().Pt();
                    m_Z2_lepminus_pt  = ele_4vec->at(idL2)->getFourVector().Pt();
                    m_Z2_lepplus_eta   = ele_4vec->at(idL1)->getFourVector().Eta();
                    m_Z2_lepminus_eta  = ele_4vec->at(idL2)->getFourVector().Eta();
                }else{
                    m_Z2_lepplus_pt   = ele_4vec->at(idL2)->getFourVector().Pt();
                    m_Z2_lepminus_pt  = ele_4vec->at(idL1)->getFourVector().Pt();
                    m_Z2_lepplus_eta   = ele_4vec->at(idL2)->getFourVector().Eta();
                    m_Z2_lepminus_eta  = ele_4vec->at(idL1)->getFourVector().Eta();
                }
                if(muon_4vec->at(idL3)->getCharge() > 0){
                    m_Z1_lepplus_pt   = muon_4vec->at(idL3)->getFourVector().Pt();
                    m_Z1_lepminus_pt  = muon_4vec->at(idL4)->getFourVector().Pt();
                    m_Z1_lepplus_eta   = muon_4vec->at(idL3)->getFourVector().Eta();
                    m_Z1_lepminus_eta  = muon_4vec->at(idL4)->getFourVector().Eta();
                }else{
                    m_Z1_lepplus_pt   = muon_4vec->at(idL4)->getFourVector().Pt();
                    m_Z1_lepminus_pt  = muon_4vec->at(idL3)->getFourVector().Pt();
                    m_Z1_lepplus_eta   = muon_4vec->at(idL4)->getFourVector().Eta();
                    m_Z1_lepminus_eta  = muon_4vec->at(idL3)->getFourVector().Eta();
                }
            }
            return true;
        }
    }
    return false;
}

bool HZZ4lHelper::MakeOutputTree(TTree& MyTree)
{
    MyTree.Branch("m4l", &m_m4l, "m4l/F");
    MyTree.Branch("Hpt", &m_Hpt_, "Hpt/F");
    MyTree.Branch("Hphi", &m_Hphi_, "Hphi/F");
    MyTree.Branch("mZ1", &m_mZ1, "mZ1/F");
    MyTree.Branch("mZ2", &m_mZ2, "mZ2/F");
    MyTree.Branch("Z1_lepplus_pt", &m_Z1_lepplus_pt, "Z1_lepplus_pt/F");
    MyTree.Branch("Z1_lepminus_pt", &m_Z1_lepminus_pt, "Z1_lepminus_pt/F");
    MyTree.Branch("Z2_lepplus_pt", &m_Z2_lepplus_pt, "Z2_lepplus_pt/F");
    MyTree.Branch("Z2_lepminus_pt", &m_Z2_lepminus_pt, "Z2_lepminus_pt/F");
    MyTree.Branch("Z1_lepplus_eta", &m_Z1_lepplus_eta, "Z1_lepplus_eta/F");
    MyTree.Branch("Z1_lepminus_eta", &m_Z1_lepminus_eta, "Z1_lepminus_eta/F");
    MyTree.Branch("Z2_lepplus_eta", &m_Z2_lepplus_eta, "Z2_lepplus_eta/F");
    MyTree.Branch("Z2_lepminus_eta", &m_Z2_lepminus_eta, "Z2_lepminus_eta/F");
    return true;
}

bool HZZ4lHelper::MakeTruthTree(TTree& MyTree)
{
    MyTree.Branch("truth_h_mass", & truth_h_mass, "truth_h_mass/F");
    MyTree.Branch("truth_z1_mass", & truth_z1_mass, "truth_z1_mass/F");
    MyTree.Branch("truth_z2_mass", & truth_z2_mass, "truth_z2_mass/F");

    MyTree.Branch("truth_h_pt", & truth_h_pt, "truth_h_pt/F");
    MyTree.Branch("truth_z1_pt", & truth_z1_pt, "truth_z1_pt/F");
    MyTree.Branch("truth_z2_pt", & truth_z2_pt, "truth_z2_pt/F");

    MyTree.Branch("truth_l1_pt", & truth_l1_pt, "truth_l1_pt/F");
    MyTree.Branch("truth_l2_pt", & truth_l2_pt, "truth_l2_pt/F");
    MyTree.Branch("truth_l3_pt", & truth_l3_pt, "truth_l3_pt/F");
    MyTree.Branch("truth_l4_pt", & truth_l4_pt, "truth_l4_pt/F");
    MyTree.Branch("truth_l1_eta", & truth_l1_eta, "truth_l1_eta/F");
    MyTree.Branch("truth_l2_eta", & truth_l2_eta, "truth_l2_eta/F");
    MyTree.Branch("truth_l3_eta", & truth_l3_eta, "truth_l3_eta/F");
    MyTree.Branch("truth_l4_eta", & truth_l4_eta, "truth_l4_eta/F");
    return true;
}

bool HZZ4lHelper::GetTruthInfo(const xAOD::TruthParticleContainer& particles)
{
    // Get truth info from the Truth container.

    auto par_itr = particles.begin();
    auto par_end = particles.end();

    bool m_debug = false;
    int n_ele = 0;
    int n_muon = 0;
    int n_tau = 0;
    for(; par_itr != par_end; par_itr++){
        const xAOD::TruthParticle* p = (*par_itr);
        int pdgId = p->pdgId();
        // cout <<"PDGID: "<< pdgId << endl;
        if(pdgId == 25) {
            truth_h_mass = p->m();
            truth_h_pt = sqrt(p->px()*p->px()+p->py()*p->py());
            if(p->hasDecayVtx()) { // find decay vertex of Higgs
                const xAOD::TruthVertex* hzz_vtx = p->decayVtx();
                int h_daugthers = hzz_vtx->nOutgoingParticles();
                if(h_daugthers != 2) {
                    if(m_debug){
                      cerr << "Higgs decays to " << h_daugthers << " daughters!" << endl;
                      const xAOD::TruthParticle* tmp = hzz_vtx->outgoingParticle(0);
                      cerr << "\t pdgID: " << tmp->pdgId() << endl; 
                      cerr << "\t mass: " << tmp->m() << endl;
                    }
                    continue;
                }
                
                for(int i_h_daughter = 0; i_h_daughter < h_daugthers; i_h_daughter++){
                    const xAOD::TruthParticle* z = hzz_vtx->outgoingParticle(i_h_daughter);
                    bool is_leading_z = true;
                    double z_pt = sqrt(z->px()*z->px()+z->py()*z->py());
                    if(z->pdgId() != 23) continue;
                    if(truth_z1_mass < 0) {
                        truth_z1_mass = z->m();
                        truth_z1_pt = z_pt;
                    } else {
                        is_leading_z = false;
                        truth_z2_mass = z->m();
                        truth_z2_pt = z_pt;
                    }
                    const xAOD::TruthVertex* zqq_vtx = z->decayVtx();
                    int z_daughters = zqq_vtx->nOutgoingParticles();
                    if(z_daughters != 2) {
                        cerr << "Z decays to " << z_daughters << " daughters!" << endl;
                        continue;
                    }
                    for(int i_z_daughter = 0; i_z_daughter < z_daughters; i_z_daughter++){
                        const xAOD::TruthParticle* quark = zqq_vtx->outgoingParticle(i_z_daughter);
                        // TVector2 quark_xy(quark->px(), quark->py());
                        // double quark_pt = quark_xy.Mod();
                        // double quark_phi = TVector2::Phi_mpi_pi(quark_xy.Phi());
                        float quark_pt = (float)quark->pt();
                        float quark_eta = (float)quark->eta();

                        int pdg = quark->absPdgId();
                        if(pdg == 11) n_ele ++;
                        else if(pdg == 13) n_muon ++;
                        else if(pdg == 15) n_tau ++;
                        else ;

                        if(is_leading_z){
                            if (truth_l1_pt < 0) {
                                truth_l1_pt = quark_pt;
                                truth_l1_eta = quark_eta;
                            }else{ 
                                truth_l2_pt = quark_pt;
                                truth_l2_eta = quark_eta;
                            }
                        }else{
                            if (truth_l3_pt < 0) {
                                truth_l3_pt = quark_pt;
                                truth_l3_eta = quark_eta;
                            } else {
                                truth_l4_pt = quark_pt;
                                truth_l4_eta = quark_eta;
                            }
                        }
                    }
                }
            }
            break;
        }
    }
    if(n_muon == 4) event_type_truth = 0;
    if(n_ele == 4) event_type_truth = 1;
    if(n_muon == 2 && n_ele == 2) event_type_truth = 2;

    if( fabs(truth_z1_mass-kZMASS) > fabs(truth_z2_mass-kZMASS)){
        swap(truth_z1_mass, truth_z2_mass);
        swap(truth_l1_pt, truth_l3_pt);
        swap(truth_l2_pt, truth_l4_pt);
        swap(truth_l1_eta, truth_l3_eta);
        swap(truth_l2_eta, truth_l4_eta);
    }
    if(truth_l1_pt < truth_l2_pt){ 
        swap(truth_l1_pt, truth_l2_pt);
        swap(truth_l1_eta, truth_l2_eta);
    }
    if(truth_l3_pt < truth_l4_pt){  
        swap(truth_l3_pt, truth_l4_pt);
        swap(truth_l3_eta, truth_l4_eta);
    }
    return true;
}

void HZZ4lHelper::Clear(){
    m_m4l = -999;
    m_Hpt_ = -999;
    m_Hphi_ = -999;
    m_mZ1   = -999;
    m_mZ2   = -999;
    m_Z1_lepplus_pt = -999;
    m_Z1_lepminus_pt    = -999;
    m_Z2_lepplus_pt = -999;
    m_Z2_lepminus_pt    = -999;
    m_Z1_lepplus_eta = -999;
    m_Z1_lepminus_eta    = -999;
    m_Z2_lepplus_eta = -999;
    m_Z2_lepminus_eta    = -999;

    truth_h_mass    = -999;
    truth_z1_mass   = -999;
    truth_z2_mass   = -999;
    truth_h_pt      = -999;
    truth_z1_pt     = -999;
    truth_z2_pt     = -999;

    truth_l1_pt     = -999;
    truth_l2_pt     = -999;
    truth_l3_pt     = -999;
    truth_l4_pt     = -999;
    truth_l1_eta     = -999;
    truth_l2_eta     = -999;
    truth_l3_eta     = -999;
    truth_l4_eta     = -999;

    event_type_truth    = -999;
}

int HZZ4lHelper::passFiducial(){
    if (event_type < 0) return 0;

    vector<float> lep_pt;
    vector<float> lep_eta;
    lep_pt.push_back(m_Z1_lepplus_pt/1E3);
    lep_pt.push_back(m_Z1_lepminus_pt/1E3);
    lep_pt.push_back(m_Z2_lepplus_pt/1E3);
    lep_pt.push_back(m_Z2_lepminus_pt/1E3);
    lep_eta.push_back(m_Z1_lepplus_eta);
    lep_eta.push_back(m_Z1_lepminus_eta);
    lep_eta.push_back(m_Z2_lepplus_eta);
    lep_eta.push_back(m_Z2_lepminus_eta);

    if(! passKinemitic(lep_pt, lep_eta, event_type)){ 
        //cout<<"ERROR:"<<m_Z1_lepplus_pt<<" " << m_Z1_lepminus_pt<<" " <<m_Z2_lepplus_pt<<" " << m_Z2_lepminus_pt << endl;
        //cout<<"\t:"<<m_Z1_lepplus_eta<<" " << m_Z1_lepminus_eta<<" " <<m_Z2_lepplus_eta<<" " << m_Z2_lepminus_eta << endl;
        return 0;
    }

    if(! passZMass(m_mZ1/1E3, m_mZ2/1E3, m_m4l/1E3)) return 0;
    return 1;
}

int HZZ4lHelper::passFiducialTruth()
{
    vector<float> lep_pt;
    vector<float> lep_eta;
    lep_pt.push_back(truth_l1_pt/1E3);
    lep_pt.push_back(truth_l2_pt/1E3);
    lep_pt.push_back(truth_l3_pt/1E3);
    lep_pt.push_back(truth_l4_pt/1E3);
    lep_eta.push_back(truth_l1_eta);
    lep_eta.push_back(truth_l2_eta);
    lep_eta.push_back(truth_l3_eta);
    lep_eta.push_back(truth_l4_eta);

    if(! passKinemitic(lep_pt, lep_eta, event_type_truth)) return 0;
    if(! passZMass(truth_z1_mass/1E3, truth_z2_mass/1E3, truth_h_mass/1E3)) return 0;
    return 1;
}

bool HZZ4lHelper::passElectronFid(float pt, float eta){
    return pt > 7 && abs(eta) < 2.47;
}
bool HZZ4lHelper::passMuonFid(float pt, float eta){
    return pt > 5 && abs(eta) < 2.7;
}
bool HZZ4lHelper::passLepPt(vector<float>& lep_pt)
{
    sort(lep_pt.begin(), lep_pt.end(), descend_sort);
    if(lep_pt.at(0) < 20) return false;
    if(lep_pt.at(1) < 15) return false;
    if(lep_pt.at(1) < 10) return false;
    return true;
}
bool HZZ4lHelper::passZMass(float z1, float z2, float m4l)
{
    // mZ1 cuts
    if(z1 >= 106 || z1 <= 50) return false; 

    // mZ2 cuts
    Double_t M_low = 0.0;
    Double_t M_high = 115.0;
    const int npoints = 2;
    Double_t mass[npoints]={140,190};
    Double_t cut[npoints]={12.0,50.};
    //loop to check the low mass cut for Z*
    for (UInt_t j = 0; j < npoints -1; j++){
        UInt_t k = j + 1;
        if (mass[k]>m4l) {
            M_low = cut[j]+(cut[k]-cut[j])/(mass[k]-mass[j])*(m4l-mass[j]);
            break;
        }
    } 
    if(m4l <= mass[0]) M_low = cut[0];
    if(m4l >= mass[1]) M_low = cut[1];
    bool _passZ2 = z2 >  M_low && z2 < M_high;
    if(!_passZ2) return false;

    return true;
}

bool HZZ4lHelper::passKinemitic(vector<float>& lep_pt, vector<float>& lep_eta, int event_type)
{
    if(event_type == 0) { // 4mu
        for(int i = 0; i < 4; i++){
            float pt = lep_pt.at(i);
            float eta = lep_eta.at(i);
            if(! passMuonFid(pt, eta)) return false;
        }
    }else if(event_type == 1){ // 4e
        for(int i = 0; i < 4; i++){
            float pt = lep_pt.at(i);
            float eta = lep_eta.at(i);
            if(! passElectronFid(pt, eta)) return false;
        }
    }else if(event_type == 2){ // 2e2mu
        for(int i = 0; i < 2; i++){
            float pt = lep_pt.at(i);
            float eta = lep_eta.at(i);
            if(! passElectronFid(pt, eta)) return false;
        }
        for(int i = 2; i < 4; i++){
            float pt = lep_pt.at(i);
            float eta = lep_eta.at(i);
            if(! passMuonFid(pt, eta)) return false;
        }
    }else if(event_type == 3){ // 2mu2e
        for(int i = 0; i < 2; i++){
            float pt = lep_pt.at(i);
            float eta = lep_eta.at(i);
            if(! passMuonFid(pt, eta)) return false;
        }
        for(int i = 2; i < 4; i++){
            float pt = lep_pt.at(i);
            float eta = lep_eta.at(i);
            if(! passElectronFid(pt, eta)) return false;
        }
    } else {
        return false;
    }
    // check pT
    if(! passLepPt(lep_pt)) return false;

    return true;
}
