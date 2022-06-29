// -*- C++ -*-
//
// Package:    AnalyzerTT
// Class:      AnalyzerTT
// 
/**\class AnalyzerTT AnalyzerTT.cc Demo/AnalyzerTT/src/AnalyzerTT.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alberto Orso Maria Iorio,40 1-B15,+41227671651,
//         Created:  dom 13 gen 2019, 19.13.15, CET
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
 
#include "DataFormats/JetReco/interface/PFJet.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"



#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"



#include "TTree.h"
#include "TFile.h"
#include "TSpline.h"
#include "TMath.h"

//
// class declaration
//

using namespace edm;
using namespace std;
using namespace reco;

class AnalyzerTT : public edm::EDAnalyzer {
public:
  explicit AnalyzerTT(const edm::ParameterSet&);
  ~AnalyzerTT();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  map< string , vector<float> > vfloats_values;
  map< string , int > sizes;
  map<string, TTree * > trees;

  bool doMCMatch_;
  
  //EDGetTokenT< std::vector< pat::Jet > > jLabel_;
  edm::InputTag m_Muons, m_gsfElectrons, m_globalMuons, theVertexLabel_, jetsLabel_, metLabel_, lhes_;
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
AnalyzerTT::AnalyzerTT(const edm::ParameterSet& iConfig)

{
  m_Muons = iConfig.getParameter<edm::InputTag>("muons");
  m_gsfElectrons = iConfig.getParameter<edm::InputTag>("electrons");
  m_globalMuons = iConfig.getParameter<edm::InputTag>("globalMuons");

  theVertexLabel_ = iConfig.getParameter<edm::InputTag>("vertexLabel");

  jetsLabel_ = iConfig.getParameter<edm::InputTag>("jets");
  metLabel_ = iConfig.getParameter<edm::InputTag>("met");

  lhes_ = iConfig.getParameter<edm::InputTag>( "lhes" );

  doMCMatch_ = iConfig.getUntrackedParameter<bool>("doMCMatch",true);

  //jLabel_             (consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetLabel"))), 
   //now do what ever initialization is needed

  
  Service<TFileService> fs;
  
  TFileDirectory ANTrees = fs->mkdir("event_trees");

  trees["events"] = new TTree("events","events");

  trees["events"]->Branch("muons_size", &sizes["muons"]);
  trees["events"]->Branch("muons_pt", &vfloats_values["muons_pt"]);
  trees["events"]->Branch("muons_eta", &vfloats_values["muons_eta"]);
  trees["events"]->Branch("muons_phi", &vfloats_values["muons_phi"]);
  trees["events"]->Branch("muons_e", &vfloats_values["muons_e"]);
  trees["events"]->Branch("muons_charge", &vfloats_values["muons_charge"]);
  trees["events"]->Branch("muons_chi2", &vfloats_values["muons_chi2"]);
  trees["events"]->Branch("muons_mctruthmatch", &vfloats_values["muons_mctruthmatch"]);
  trees["events"]->Branch("muons_minDR", &vfloats_values["muons_minDR"]);

  trees["events"]->Branch("muons_dB", &vfloats_values["muons_dB"]);
  trees["events"]->Branch("muons_dz", &vfloats_values["muons_dz"]);
  trees["events"]->Branch("muons_reliso", &vfloats_values["muons_reliso"]);

  trees["events"]->Branch("muons_isLooseMuon", &vfloats_values["muons_isLooseMuon"]);
  trees["events"]->Branch("muons_isTightMuon", &vfloats_values["muons_isTightMuon"]);

  trees["events"]->Branch("muons_nmatchstat", &vfloats_values["muons_nmatchstat"]);
  trees["events"]->Branch("muons_nmuonhits", &vfloats_values["muons_nmuonhits"]);
  trees["events"]->Branch("muons_npixhits", &vfloats_values["muons_npixhits"]);
  trees["events"]->Branch("muons_nlaywithmeas", &vfloats_values["muons_nlaywithmeas"]);


  trees["events"]->Branch("electrons_size", &sizes["electrons"]);
  trees["events"]->Branch("electrons_pt", &vfloats_values["electrons_pt"]);
  trees["events"]->Branch("electrons_eta", &vfloats_values["electrons_eta"]);
  trees["events"]->Branch("electrons_phi", &vfloats_values["electrons_phi"]);
  trees["events"]->Branch("electrons_e", &vfloats_values["electrons_e"]);
  trees["events"]->Branch("electrons_mctruthmatch", &vfloats_values["electrons_mctruthmatch"]);
  trees["events"]->Branch("electrons_minDR", &vfloats_values["electrons_minDR"]);

  trees["events"]->Branch("electrons_isTightElectron", &vfloats_values["electrons_isTightElectron"]);
  trees["events"]->Branch("electrons_reliso", &vfloats_values["electrons_reliso"]);
  trees["events"]->Branch("electrons_hOverE", &vfloats_values["electrons_hOverE"]);
  trees["events"]->Branch("electrons_dB", &vfloats_values["electrons_dB"]);
  trees["events"]->Branch("electrons_dz", &vfloats_values["electrons_dz"]);
  trees["events"]->Branch("electrons_dEtaIn", &vfloats_values["electrons_dEtaIn"]);
  trees["events"]->Branch("electrons_dPhiIn", &vfloats_values["electrons_dPhiIn"]);
  trees["events"]->Branch("electrons_ninner", &vfloats_values["electrons_ninner"]);
  trees["events"]->Branch("electrons_sigmaieie", &vfloats_values["electrons_sigmaieie"]);
  trees["events"]->Branch("electrons_ooemoop", &vfloats_values["electrons_ooemoop"]);
 
 trees["events"]->Branch("met_pt", &vfloats_values["met_pt"]);
  trees["events"]->Branch("met_phi", &vfloats_values["met_phi"]);

  trees["events"]->Branch("jets_size", &sizes["jets"]);
  trees["events"]->Branch("jets_pt", &vfloats_values["jets_pt"]);
  trees["events"]->Branch("jets_eta", &vfloats_values["jets_eta"]);
  trees["events"]->Branch("jets_phi", &vfloats_values["jets_phi"]);
  trees["events"]->Branch("jets_e", &vfloats_values["jets_e"]);
  trees["events"]->Branch("muons_c", &vfloats_values["muons_c"]);

  trees["events"]->Branch("genlep_size", &sizes["genleps"]);
  trees["events"]->Branch("genlep_pt", &vfloats_values["genlep_pt"]);
  trees["events"]->Branch("genlep_eta", &vfloats_values["genlep_eta"]);
  trees["events"]->Branch("genlep_phi", &vfloats_values["genlep_phi"]);
  trees["events"]->Branch("genlep_e", &vfloats_values["genlep_e"]);
  trees["events"]->Branch("genlep_id", &vfloats_values["genlep_id"]);
  //  trees["events"]->Branch("genlep_isFromZ", &vfloats_values["genlep_isFromZ"]);

  trees["events"]->Branch("electrons_charge", &vfloats_values["electrons_charge"]);
  
  //TFileDirectory DMTrees;// = fs->mkdir( "systematics_trees" );                                                                                                                                                     
}




AnalyzerTT::~AnalyzerTT()
{
 
  
  

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
AnalyzerTT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   reco::Vertex::Point posVtx;
   reco::Vertex::Error errVtx;
   unsigned int theIndexOfThePrimaryVertex = 999.;

   
   edm::Handle<LHEEventProduct > lhes;
   if(doMCMatch_){
     iEvent.getByLabel(lhes_, lhes);
   }

   edm::Handle<reco::VertexCollection> vertex;
   iEvent.getByLabel(theVertexLabel_, vertex);
   if (vertex.isValid()){
     for (unsigned int ind=0; ind<vertex->size(); ++ind) {
       if ( (*vertex)[ind].isValid() && !((*vertex)[ind].isFake()) ) {
	 theIndexOfThePrimaryVertex = ind;
	 break;
       }
     }
   }

   vfloats_values["muons_pt"].clear();
   vfloats_values["muons_eta"].clear();
   vfloats_values["muons_phi"].clear();
   vfloats_values["muons_e"].clear();
   vfloats_values["muons_charge"].clear();
   vfloats_values["muons_dB"].clear();
   vfloats_values["muons_dz"].clear();
   vfloats_values["muons_chi2"].clear();
   vfloats_values["muons_reliso"].clear();
   vfloats_values["muons_mctruthmatch"].clear();
   vfloats_values["muons_minDR"].clear();
   
   
   vfloats_values["muons_isLooseMuon"].clear();
   vfloats_values["muons_isTightMuon"].clear();

   vfloats_values["muons_nmatchstat"].clear();
   vfloats_values["muons_nmuonhits"].clear();
   vfloats_values["muons_npixhits"].clear();
   vfloats_values["muons_nlaywithmeas"].clear();
   
   vfloats_values["electrons_pt"].clear();
   vfloats_values["electrons_eta"].clear();
   vfloats_values["electrons_phi"].clear();
   vfloats_values["electrons_e"].clear();
   vfloats_values["electrons_minDR"].clear();
   vfloats_values["electrons_charge"].clear();
   vfloats_values["electrons_mctruthmatch"].clear();

   vfloats_values["electrons_isTightElectron"].clear();
   vfloats_values["electrons_reliso"].clear();
   vfloats_values["electrons_hOverE"].clear();
   vfloats_values["electrons_dB"].clear();
   vfloats_values["electrons_dz"].clear();

   vfloats_values["electrons_dEtaIn"].clear();
   vfloats_values["electrons_dPhiIn"].clear();
   vfloats_values["electrons_ninner"].clear();
   vfloats_values["electrons_sigmaieie"].clear();
   vfloats_values["electrons_ooemoop"].clear();

/*  trees["events"]->Branch("electrons_dz", &vfloats_values["electrons_dz"]);
  trees["events"]->Branch("electrons_dEtaIn", &vfloats_values["electrons_dEtaIn"]);
  trees["events"]->Branch("electrons_dPhiIn", &vfloats_values["electrons_dPhiIn"]);
  trees["events"]->Branch("electrons_ninner", &vfloats_values["electrons_ninner"]);
  trees["events"]->Branch("electrons_sigmaieie", &vfloats_values["electrons_sigmaieie"]);
  trees["events"]->Branch("electrons_ooemoop", &vfloats_values["electrons_ooemoop"]);
*/

   vfloats_values["met_pt"].clear();
   vfloats_values["met_phi"].clear();
   
   vfloats_values["jets_pt"].clear();
   vfloats_values["jets_eta"].clear();
   vfloats_values["jets_phi"].clear();
   vfloats_values["jets_e"].clear();
   
   vfloats_values["genlep_pt"].clear();
   vfloats_values["genlep_eta"].clear();
   vfloats_values["genlep_phi"].clear();
   vfloats_values["genlep_e"].clear();
   vfloats_values["genlep_id"].clear();
   vfloats_values["genlep_status"].clear();
   vfloats_values["genlep_isFromZ"].clear();
   vfloats_values["genlep_isFromW"].clear();
   


   edm::Handle<reco::MuonCollection> Muons;
   iEvent.getByLabel(m_Muons, Muons);

   edm::Handle<GsfElectronCollection> gsfElectrons;
   iEvent.getByLabel(m_gsfElectrons,gsfElectrons);
   

   edm::Handle<reco::PFJetCollection> jets;
   iEvent.getByLabel(jetsLabel_, jets);

   edm::Handle<reco::PFMETCollection> pfmet;
   iEvent.getByLabel(metLabel_, pfmet);
   
   //   std::cout << " muons size "<< Muons->size()<<std::endl;
   //std::cout << " electrons size "<< gsfElectrons->size()<<std::endl;
   sizes["muons"]=0;
   sizes["electrons"]=0;
   sizes["jets"]=0;
   sizes["genleps"]=0;
   //   size_t mu=0;gsfElectrons->oisize()
   size_t nup=0;
   if(doMCMatch_){

     nup=lhes->hepeup().NUP;
     sizes["lhes"]=nup;

    //cout << " mark 2.2 ; nup"<<nup<<endl; 
    for( size_t i=0;i<nup;++i){
      //cout << "  " << i +1 << endl;
      int id = lhes->hepeup().IDUP[i];
      float px = lhes->hepeup().PUP[i][0];
      float py = lhes->hepeup().PUP[i][1];
      float pz = lhes->hepeup().PUP[i][2];
      float energy = lhes->hepeup().PUP[i][3];

      int istatus = lhes->hepeup().ISTUP[i];
      int imoth1 =lhes->hepeup().MOTHUP[i].first;
      //      int imoth2 =lhes->hepeup().MOTHUP[i].second;

      float isFromWChain=false;
      float isFromZChain=false;

      if(abs(lhes->hepeup().IDUP[imoth1-1])==24 && istatus==1){
	isFromWChain=true;
      }
      if(abs(lhes->hepeup().IDUP[imoth1-1])==23 && istatus==1){
	isFromZChain=true;
      }
      
      
      math::XYZTLorentzVector part = math::XYZTLorentzVector(px, py, pz, energy);
      float pt = part.pt();
      float phi = part.phi();
      float eta = part.eta();
      
      cout << " id "<<id << " status "<< istatus << " eta "<< eta << " phi "<< phi <<endl;
      if(abs(id)==11 || abs(id)==13 || abs(id)==15){
	vfloats_values["genlep_pt"].push_back(pt);
	vfloats_values["genlep_eta"].push_back(eta);
	vfloats_values["genlep_phi"].push_back(phi);
	vfloats_values["genlep_e"].push_back(energy);
	vfloats_values["genlep_status"].push_back(istatus);
	vfloats_values["genlep_id"].push_back(id);
	vfloats_values["genlep_isFromZ"].push_back(isFromZChain);      
	vfloats_values["genlep_isFromW"].push_back(isFromWChain);
	sizes["genleps"]+=1;
      }
    }
    
    
   }
   bool passes=true;
     
   for (reco::MuonCollection::const_iterator globalMuon = Muons->begin();  globalMuon != Muons->end();  ++globalMuon) {
     //     std::cout << "muon pt "<< globalMuon->combinedMuon()->pt() << " eta "<<  globalMuon->combinedMuon()->eta() << " phi "<<  globalMuon->combinedMuon()->phi()<<std::endl;
     //     std::cout << "muon pt "<< globalMuon->pt() << " eta "<<  globalMuon->eta() << " phi "<<  globalMuon->phi()<<std::endl;
     //  if(globalMuon->pt()<4)continue;

     bool isGlobalMuon = globalMuon-> isGlobalMuon();
     bool isTrackerMuon = globalMuon-> isTrackerMuon();
     if(isGlobalMuon && false)cout<<"isglobalmuon" <<endl;
     // selection part
     bool isTightMuon=true;
     bool isLooseMuon=true;
     isTightMuon=isTightMuon && globalMuon->pt()>20;
     isTightMuon=isTightMuon && fabs(globalMuon->eta())<2.4;
     isTightMuon=isTightMuon&& (isGlobalMuon && isTrackerMuon); //globalmuon condition 
     isLooseMuon=isTightMuon;
     if(!isTightMuon)continue;// all muons have pt > 20 
     
     float db=fabs(globalMuon->innerTrack()->dxy(vertex->at(theIndexOfThePrimaryVertex).position()));
     float dz= fabs(globalMuon->innerTrack()->dz(vertex->at(theIndexOfThePrimaryVertex).position()));
     
     float iso = (globalMuon->pfIsolationR03().sumChargedHadronPt + std::max(0.,globalMuon->pfIsolationR03().sumNeutralHadronEt + globalMuon->pfIsolationR03().sumPhotonEt - (0.5 * globalMuon->pfIsolationR03().sumPUPt)));
     
     float reliso = iso/globalMuon->pt();
     float chi2 = globalMuon->globalTrack()->normalizedChi2();     

     float matchstat = globalMuon->numberOfMatchedStations();
     float nmuonhits = globalMuon->globalTrack()->hitPattern().numberOfValidMuonHits();
     float npixhits = globalMuon->innerTrack()->hitPattern().numberOfValidPixelHits() ;
     float nlaywithmeas = globalMuon->track()->hitPattern().trackerLayersWithMeasurement();

     //isTightMuon=isTightMuon&& (isGlobalMuon && isTrackerMuon); //globalmuon condition 
     isTightMuon=isTightMuon && db<0.02; //dxy condition 
     isTightMuon=isTightMuon && dz<0.1; //dz condition
     isTightMuon=isTightMuon && chi2<10; //chi2 condition
     isTightMuon=isTightMuon && reliso<0.15; //chi2 condition

     isTightMuon=isTightMuon && matchstat >1 && nmuonhits>0 && npixhits>0 &&nlaywithmeas>8;
     //end conditions
    
     sizes["muons"]++;

    
     float isMCMatch=0.;
    
     vfloats_values["muons_pt"].push_back(globalMuon->pt());
     vfloats_values["muons_eta"].push_back(globalMuon->eta());
     vfloats_values["muons_phi"].push_back(globalMuon->phi());
     vfloats_values["muons_e"].push_back(globalMuon->energy());
     vfloats_values["muons_charge"].push_back(globalMuon->charge());
     //     cout << " muon # "<< sizes["muons"] << " pt is "<< globalMuon->pt() <<" chi2 is "<< globalMuon->muonBestTrack()->normalizedChi2()<< endl;
 
     float pt=globalMuon->pt();
     float eta=globalMuon->eta();
     float phi=globalMuon->phi();
     float e=globalMuon->energy();

     float minDR =999.;
     if(doMCMatch_){
       for(size_t mc=0;mc< vfloats_values["genlep_pt"].size();++mc){
	 
	 float ptg = vfloats_values["genlep_pt"].at(mc);
	 float etag = vfloats_values["genlep_eta"].at(mc);
	 float phig = vfloats_values["genlep_phi"].at(mc);
	 float eg = vfloats_values["genlep_e"].at(mc);
	 float idg = vfloats_values["genlep_id"].at(mc);
	 float statusg = vfloats_values["genlep_status"].at(mc);
	 //	 float fromZ = vfloats_values["genlep_isFroMZ"].at(mc);
	 
	 float dr =0.; 
	 dr=deltaR(math::PtEtaPhiELorentzVector(pt, eta, phi, e),math::PtEtaPhiELorentzVector(ptg, etag, phig, eg));
	 
	 if(abs(idg)==13. && statusg==1.0){
	   cout << " muon phig "<< phig <<" phim "<< phi <<" etag "<< etag << " etam "<< eta<<  endl;
	   minDR=min(minDR,dr);
	 }
	 if(dr<0.4 && abs(idg)==13. && statusg==1.0 ){
	   cout << " muon drmatch ! dr " << dr << "phig "<< phig <<" phim "<< phi <<" etag "<< etag << " etam "<< eta<<  endl;
	   isMCMatch=true;
	 }       
       }
     }
     if(minDR>10)minDR=10.;
     vfloats_values["muons_minDR"].push_back(minDR);
    
     
     
     vfloats_values["muons_chi2"].push_back(chi2);
     vfloats_values["muons_dB"].push_back( db);
     vfloats_values["muons_dz"].push_back( dz);
     vfloats_values["muons_reliso"].push_back(reliso);

     //float matchstat = globalMuon->numberOfMatchedStations();
     //float nmuonhits = globalMuon->globalTrack()->hitPattern().numberOfValidMuonHits();
     //float npixhit = globalMuon->innerTrack()->hitPattern().numberOfValidPixelHits() ;
     //float nlaywithmeas = globalMuon->track()->hitPattern().trackerLayersWithMeasurement();
     
     vfloats_values["muons_matchstat"].push_back(matchstat);
     vfloats_values["muons_nmuonhits"].push_back(nmuonhits);
     vfloats_values["muons_npixhits"].push_back(npixhits);
     vfloats_values["muons_nlaywithmeas"].push_back(nlaywithmeas);

     vfloats_values["muons_isLooseMuon"].push_back(isLooseMuon);
     vfloats_values["muons_isTightMuon"].push_back(isTightMuon);
    
     vfloats_values["muons_mctruthmatch"].push_back(isMCMatch);
     //     fabs(globalMuon->muonBestTrack()->dxy(vertex->at(theIndexOfThePrimaryVertex).position()));
     
     
     //     cout << " iso is "<< ((globalMuon->pfIsolationR03().sumChargedHadronPt + std::max(0.,globalMuon->pfIsolationR03().sumNeutralHadronEt + globalMuon->pfIsolationR03().sumPhotonEt - (0.5 * globalMuon->pfIsolationR03().sumPUPt))) / globalMuon->pt())<< endl;
        
   

   }

   for (size_t i =0; i <vfloats_values["muons_c"].size(); ++i ){
     cout << i <<" i " << vfloats_values["muons_c"].at(i) << endl;
       }
   
   for (reco::GsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin();gsfIter!=gsfElectrons->end(); gsfIter++){
     //int mishits = gsfIter->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
     float ninner=10.;
     float db=10.,dz=10.;
     //     if(ack() && false){
     ninner = gsfIter->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
     db=fabs(gsfIter->gsfTrack()->dxy(vertex->at(theIndexOfThePrimaryVertex).position()));
     dz=fabs(gsfIter->gsfTrack()->dz(vertex->at(theIndexOfThePrimaryVertex).position()));
     
       ///     }

     double tkIso = gsfIter->dr03TkSumPt();
     double ecalIso = gsfIter->dr03EcalRecHitSumEt();
     double hcalIso = gsfIter->dr03HcalTowerSumEt();
     //ouble eOverP = gsfIter->eSuperClusterOverP();
     //    double fBrem = gsfIter->fbrem();
     double hOverE = gsfIter->hadronicOverEm();
     double sigmaieie = gsfIter->sigmaIetaIeta();
     // preselect electrons
     // cout << " iso "<<(tkIso+hcalIso+ecalIso)/gsfIter->pt()<< " eop "<<eOverP <<" hoe "<< hOverE << " sieie "<< sigmaieie<<endl;

     bool isTightElectron = true;
     isTightElectron = isTightElectron && gsfIter->pt()>20; 
     isTightElectron = isTightElectron && fabs(gsfIter->eta())<(2.5); 
     if(!isTightElectron) continue;
     
     
     float dPhiIn=fabs(gsfIter->deltaPhiSuperClusterTrackAtVtx());
     float dEtaIn=fabs(gsfIter->deltaEtaSuperClusterTrackAtVtx());
       //     if(eOverP>0.132)continue;
     float EIn=gsfIter->ecalEnergy();
     float PIn=EIn/gsfIter->eSuperClusterOverP();
     //float eop=EIn/PIn;
     float ooemoop= fabs(1/EIn -1/PIn);
     

     //selection part
     isTightElectron=isTightElectron && db<0.04;
     isTightElectron=isTightElectron && dz<0.1;


     isTightElectron=isTightElectron && hOverE<0.12;
     isTightElectron=isTightElectron && ooemoop<0.05;
     isTightElectron=isTightElectron && sigmaieie<0.01;
     

     isTightElectron=isTightElectron && dEtaIn<0.007;
     isTightElectron=isTightElectron && dPhiIn<0.15;
     isTightElectron=isTightElectron && ninner==0;
     
     
     float reliso = ((tkIso+hcalIso+ecalIso)/gsfIter->pt());
     isTightElectron=isTightElectron && reliso < 0.15;
     
     sizes["electrons"]++;
     
     float isMCMatch=0.;
     //Cout () << " Gsf () pt "<<gsfIter->pt()<< "scpt" <<gsfIter->superCluster()->energy()<<endl;
     vfloats_values["electrons_pt"].push_back(gsfIter->pt());
     vfloats_values["electrons_eta"].push_back(gsfIter->eta());
     vfloats_values["electrons_phi"].push_back(gsfIter->phi());
     vfloats_values["electrons_e"].push_back(gsfIter->energy());
     vfloats_values["electrons_charge"].push_back(gsfIter->charge());

     vfloats_values["electrons_isTightElectron"].push_back(isTightElectron);

     vfloats_values["electrons_reliso"].push_back(reliso);
     vfloats_values["electrons_ooemoop"].push_back(ooemoop);
     vfloats_values["electrons_hOverE"].push_back(hOverE);
     vfloats_values["electrons_sigmaieie"].push_back(sigmaieie);

     vfloats_values["electrons_dB"].push_back(db);
     vfloats_values["electrons_dz"].push_back(dz);
     vfloats_values["electrons_dEtaIn"].push_back(dEtaIn);
     vfloats_values["electrons_dPhiIn"].push_back(dPhiIn);
     vfloats_values["electrons_ninner"].push_back(ninner);
 
     
     float pt=gsfIter->pt();
     float eta=gsfIter->eta();
     float phi=gsfIter->phi();
     float e=gsfIter->energy();
     
     float minDR =999.;
     if(doMCMatch_){
       for(size_t mc=0;mc< vfloats_values["genlep_pt"].size();++mc){
	 
	 float ptg = vfloats_values["genlep_pt"].at(mc);
	 float etag = vfloats_values["genlep_eta"].at(mc);
	 float phig = vfloats_values["genlep_phi"].at(mc);
	 float eg = vfloats_values["genlep_e"].at(mc);
	 float idg = vfloats_values["genlep_id"].at(mc);
	 float statusg = vfloats_values["genlep_status"].at(mc);
	 //	 float fromZ = vfloats_values["genlep_isFroMZ"].at(mc);
	 
	 float dr =0.; 
	 dr=deltaR(math::PtEtaPhiELorentzVector(pt, eta, phi, e),math::PtEtaPhiELorentzVector(ptg, etag, phig, eg));
	 
	 if(abs(idg)==11. && statusg==1.0){
	   cout << " ele phig "<< phig <<" phim "<< phi <<" etag "<< etag << " etam "<< eta<<  endl;
	   minDR=min(minDR,dr);
	 }
	 if(dr<0.4 && abs(idg)==11. && statusg==1.0 ){
	   cout << " ele drmatch ! dr " << dr << "phig "<< phig <<" phim "<< phi <<" etag "<< etag << " etam "<< eta<<  endl;
	   isMCMatch=true;
	 }       
       }
     }
     if(minDR>10)minDR=10.;
     vfloats_values["electrons_minDR"].push_back(minDR);
     
     vfloats_values["electrons_mctruthmatch"].push_back(isMCMatch);
     
     //     vfloats_values["electrons_iso"].push_back();
     //vfloats_values["electrons_sigmaieie"].push_back(sigmaieie);
     //vfloats_values["electrons_energyOverP"].push_back(eOverP);
     //vfloats_values["electrons_hadronicOverEM"].push_back(hOverE);
     
   }
   edm::Handle<reco::TrackCollection> globalMuons;
   iEvent.getByLabel(m_globalMuons, globalMuons);
   
   //edm::Handle<reco::PFJetCollection> jets;

   for (reco::PFJetCollection::const_iterator jet = jets->begin();  jet != jets->end();  ++jet) {
     // std::cout << "jet pt "<< jet->pt() << " eta "<<  jet->eta() << " phi "<<  jet->phi()<<std::endl;

     if(jet->pt()<20)continue;
     if(abs(jet->eta())>2.5)continue;
     sizes["jets"]++;
     vfloats_values["jets_pt"].push_back(jet->pt());
     vfloats_values["jets_eta"].push_back(jet->eta());
     vfloats_values["jets_phi"].push_back(jet->phi());
     vfloats_values["jets_e"].push_back(jet->energy());
     
     //std::cout << "muon pt "<< globalMuon->pt() << " eta "<<  globalMuon->eta() << " phi "<<  globalMuon->phi()<<std::endl;
   }

   vfloats_values["met_pt"].push_back(pfmet->front().pt());
   vfloats_values["met_phi"].push_back(pfmet->front().phi());
   
   if((sizes["muons"]+sizes["electrons"])<2)passes=false; 
   passes=true;
   if(passes)trees["events"]->Fill();
   //  edm::Handle<std::vector<pat::Jet> > jetHandle, packedjetHandle;
   //iEvent.getByToken(jLabel_, jetHandle);
   //std::unique_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );
   //for (size_t i = 0; i< jetColl->size(); i++){
   // pat::Jet & jet = (*jetColl)[i];
   //}

}


// ------------ method called once each job just before starting event loop  ------------
void 
AnalyzerTT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AnalyzerTT::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
AnalyzerTT::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
AnalyzerTT::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
AnalyzerTT::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
AnalyzerTT::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnalyzerTT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzerTT);
