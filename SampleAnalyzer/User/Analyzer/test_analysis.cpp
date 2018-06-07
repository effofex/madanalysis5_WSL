#include "SampleAnalyzer/User/Analyzer/test_analysis.h"
#include <cmath>     //min() and abs()
#include <vector>
#include <algorithm> //std::set_difference
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool test_analysis::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;
  // initialize variables, histos
  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void test_analysis::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  // saving histos
  cout << "END   Finalization" << endl;
}

void test_analysis::filterBaseLine(const vector<RecLeptonFormat> &from, vector<RecLeptonFormat> &to, 
                    double min_pt, double maxAbsEta){
  copy_if(from.begin(),from.end(),
	back_inserter(to),
        [](const RecLeptonFormat& l){return (l.pt() > 5) && (abs(l.eta())<2.47);});
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool test_analysis::Execute(SampleFormat& sample, const EventFormat& event)
{
  if(0 != event.rec()){    
    vector<RecLeptonFormat> baseLineElecs;
    vector<RecLeptonFormat> baseLineMuons;
    vector<RecJetFormat> baseLineJets;    
    vector<RecLeptonFormat> signalMuons;
    vector<RecLeptonFormat> signalElecs;
    vector<RecJetFormat> signalJets;
    vector<MAuint32> removeJetIndex; 
    vector<MAuint32> removeElecIndex;
    vector<MAuint32> removeMuonIndex;

/***************************************************
Create baseline vectors
***************************************************/
    filterBaseLine(event.rec()->electrons(),baseLineElecs,5,2.47);
    //Example of looping thru with auto
    /*for( auto &i: baseLineElecs)
    {
      cout << i.pt() << endl;
    }*/
    
    // Looking through the reconstructed muon collection
    for (MAuint32 i=0;i<event.rec()->muons().size();i++)
    {
      const RecLeptonFormat& muon = event.rec()->muons()[i];
      
      const double momentumBound  = 4; //GeV
      const double psuedoRapidityBound = 2.7; 
      if(momentumBound < muon.pt() && psuedoRapidityBound >= abs(muon.eta())){
        baseLineMuons.push_back(muon);
      }
    }
 
    //Looking through reconstructed jet collection
    for(MAuint32 i=0;i<event.rec()->jets().size();i++)
    {
      const RecJetFormat& jet = event.rec()->jets()[i];

      const double momentumBound = 25; //GeV
      const double psuedoRapidityBound = 2.5;

      if(momentumBound < jet.pt() && psuedoRapidityBound >= abs(jet.eta())){
        baseLineJets.push_back(jet);
      }
    }

/***********************************
Determine overlapping elements to remove
***********************************/
    
    //Overlap removal between electrons and jets, partial implementation of col
    //#2 in Table 3 of CERN-ERP-2017-246 https://arxiv.org/pdf/1711.11520.pdf
    for(MAuint32 i=0;i<baseLineJets.size();i++)
    {
      
      for(MAuint32 j=0;j<baseLineElecs.size();j++)
      {
        const RecJetFormat& myjet = baseLineJets[i];
        const RecLeptonFormat& myelec = baseLineElecs[j];
        if(myjet.dr(myelec) < 0.2)
        {
          removeJetIndex.push_back(i);
        }

      }
    }
    
    //Overlap removal between electrons and jets according to criterion 3
    for(MAuint32 i=0;i<baseLineJets.size();i++)
    {
      for(MAuint32 j=0;j<baseLineElecs.size();j++)
      {
        const RecJetFormat& myjet = baseLineJets[i];
        const RecLeptonFormat& myelec = baseLineElecs[j];
        double minDist = min(0.4,0.04 + 10/(myelec.pt()));
	if(myjet.dr(myelec) < minDist)
        {
          removeElecIndex.push_back(j);
	}
         
      }
    }
    
    for(MAuint32 i=0;i<baseLineJets.size();i++)
    {
      for(MAuint32 j=0;j<baseLineMuons.size();j++)
      {
        const RecJetFormat& myjet = baseLineJets[i];
        const RecLeptonFormat& mymuon = baseLineMuons[j];
        double minDist = min(0.4,0.04 + 10/(mymuon.pt()));
	if(myjet.dr(mymuon) < minDist)
        {
          removeMuonIndex.push_back(j);
	}
         
      }
    }
   
/*************************
Create the signal vectors based on items flagged for removal.
Note DOES NOT actually remove anything from the baseline
**************************/
    for(MAuint32 i=0;i<baseLineJets.size();i++)
    {
      bool removeMe = false;
      for(MAuint32 j=0;j<removeJetIndex.size();j++)
      {
	removeMe = removeMe || (i==removeJetIndex[j]);
      }
      if(!removeMe){
          signalJets.push_back(baseLineJets[i]);
      }
    }
    
    for(MAuint32 i=0;i<baseLineElecs.size();i++)
    {
      bool removeMe = false;
      for(MAuint32 j=0;j<removeElecIndex.size();j++)
      {
	removeMe = removeMe || (i==removeElecIndex[j]);
      }
      if(!removeMe){
          signalElecs.push_back(baseLineElecs[i]);
      }
   }
   
   for(MAuint32 i=0;i<baseLineMuons.size();i++)
    {
      bool removeMe = false;
      for(MAuint32 j=0;j<removeMuonIndex.size();j++)
      {
	removeMe = removeMe || (i==removeMuonIndex[j]);
      }
      if(!removeMe){
          signalMuons.push_back(baseLineMuons[i]);
      }
   }

/*******************
Output the results
********************/ 
    cout << "Base and signal vector sizes" << endl;
    cout << "Particle\tBase\tRemoved\tSignal" << endl;
    cout << "\tJets:\t" << baseLineJets.size() << "\t" << removeJetIndex.size() << "\t" << signalJets.size() << endl;
    cout << "\tElecs:\t" << baseLineElecs.size() << "\t" << removeElecIndex.size() << "\t" << signalElecs.size() << endl;
    cout << "\tMuons:\t" << baseLineMuons.size() << "\t" << removeMuonIndex.size() << "\t" << signalMuons.size() << endl;

}


  return true; 
// ***************************************************************************
  // Example of analysis with generated particles
  // Concerned samples : LHE/STDHEP/HEPMC
  // ***************************************************************************
  /* 
  if (event.mc()!=0)
  {

    cout << "---------------NEW EVENT-------------------" << endl;

    // Event weight
    double myWeight=1.;
    if (!Configuration().IsNoEventWeight()) myWeight=event.mc()->weight();

    // Initial state
    for (MAuint32 i=0;i<event.mc()->particles().size();i++)
    {
      const MCParticleFormat& part = event.mc()->particles()[i];

      cout << "----------------------------------" << endl;
      cout << "MC particle" << endl;
      cout << "----------------------------------" << endl;

      // display index particle
      cout << "index=" << i+1;

      // display the status code
      cout << "Status Code=" << part.statuscode() << endl;
      if (PHYSICS->Id->IsInitialState(part)) cout << " (Initial state) ";
      else if (PHYSICS->Id->IsFinalState(part)) cout << " (Final state) ";
      else cout << " (Intermediate state) ";
      cout << endl;

      // pdgid
      cout << "pdg id=" << part.pdgid() << endl;
      if (PHYSICS->Id->IsInvisible(part)) cout << " (invisible particle) ";
      else cout << " (visible particle) ";
      cout << endl;

      // display kinematics information
      cout << "px=" << part.px()
                << " py=" << part.py()
                << " pz=" << part.pz()
                << " e="  << part.e()
                << " m="  << part.m() << endl;
      cout << "pt=" << part.pt() 
                << " eta=" << part.eta() 
                << " phi=" << part.phi() << endl;

      // display particle mother id
      if (part.mothers().empty()) 
      {
        cout << "particle with no mother." << endl;
      }
      else
      {
        std::cout << "particle coming from the decay of";
        for(MAuint32 j=0;j<part.mothers().size();j++)
        {
          const MCParticleFormat* mother = part.mothers()[j];
          cout << " " << mother->pdgid();
        }
        std::cout << "." << endl;
      }
    }

    // Transverse missing energy (MET)
    cout << "MET pt=" << event.mc()->MET().pt()
         << " phi=" << event.mc()->MET().phi() << endl;
    cout << endl;

    // Transverse missing hadronic energy (MHT)
    cout << "MHT pt=" << event.mc()->MHT().pt()
              << " phi=" << event.mc()->MHT().phi() << endl;
    cout << endl;

    // Total transverse energy (TET)
    cout << "TET=" << event.mc()->TET() << endl;
    cout << endl;

    // Total transverse hadronic energy (THT)
    cout << "THT=" << event.mc()->THT() << endl; 
   cout << endl;

  return true;
  }
  */ 


  // ***************************************************************************
  // Example of analysis with reconstructed objects
  // Concerned samples : 
  //   - LHCO samples
  //   - LHE/STDHEP/HEPMC samples after applying jet-clustering algorithm
  // ***************************************************************************
  /*
  // Event weight
  double myWeight=1.;
  if (!Configuration().IsNoEventWeight() && event.mc()!=0) myWeight=event.mc()->weight();

  Manager()->InitializeForNewEvent(myWeight);

  if (event.rec()!=0)
  {
    cout << "---------------NEW EVENT-------------------" << endl;

    // Looking through the reconstructed electron collection
    for (MAuint32 i=0;i<event.rec()->electrons().size();i++)
    {
      const RecLeptonFormat& elec = event.rec()->electrons()[i];
      cout << "----------------------------------" << endl;
      cout << "Electron" << endl;
      cout << "----------------------------------" << endl;
      cout << "index=" << i+1 
                << " charge=" << elec.charge() << endl;
      cout << "px=" << elec.px()
                << " py=" << elec.py()
                << " pz=" << elec.pz()
                << " e="  << elec.e()
                << " m="  << elec.m() << endl;
      cout << "pt=" << elec.pt() 
                << " eta=" << elec.eta() 
                << " phi=" << elec.phi() << endl;
      cout << "pointer address to the matching MC particle: " 
                << elec.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed muon collection
    for (MAuint32 i=0;i<event.rec()->muons().size();i++)
    {
      const RecLeptonFormat& mu = event.rec()->muons()[i];
      cout << "----------------------------------" << endl;
      cout << "Muon" << endl;
      cout << "----------------------------------" << endl;
      cout << "index=" << i+1 
                << " charge=" << mu.charge() << endl;
      cout << "px=" << mu.px()
                << " py=" << mu.py()
                << " pz=" << mu.pz()
                << " e="  << mu.e()
                << " m="  << mu.m() << endl;
      cout << "pt=" << mu.pt() 
                << " eta=" << mu.eta() 
                << " phi=" << mu.phi() << endl;
      cout << "ET/PT isolation criterion =" << mu.ET_PT_isol() << endl;
      cout << "pointer address to the matching MC particle: " 
           << mu.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed hadronic tau collection
    for (MAuint32 i=0;i<event.rec()->taus().size();i++)
    {
      const RecTauFormat& tau = event.rec()->taus()[i];
      cout << "----------------------------------" << endl;
      cout << "Tau" << endl;
      cout << "----------------------------------" << endl;
      cout << "tau: index=" << i+1 
                << " charge=" << tau.charge() << endl;
      cout << "px=" << tau.px()
                << " py=" << tau.py()
                << " pz=" << tau.pz()
                << " e="  << tau.e()
                << " m="  << tau.m() << endl;
      cout << "pt=" << tau.pt() 
                << " eta=" << tau.eta() 
                << " phi=" << tau.phi() << endl;
      cout << "pointer address to the matching MC particle: " 
           << tau.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed jet collection
    for (MAuint32 i=0;i<event.rec()->jets().size();i++)
    {
      const RecJetFormat& jet = event.rec()->jets()[i];
      cout << "----------------------------------" << endl;
      cout << "Jet" << endl;
      cout << "----------------------------------" << endl;
      cout << "jet: index=" << i+1 
           << " charge=" << jet.charge() << endl;
      cout << "px=" << jet.px()
           << " py=" << jet.py()
           << " pz=" << jet.pz()
           << " e="  << jet.e()
           << " m="  << jet.m() << endl;
      cout << "pt=" << jet.pt() 
           << " eta=" << jet.eta() 
           << " phi=" << jet.phi() << endl;
      cout << "b-tag=" << jet.btag()
           << " true b-tag (before eventual efficiency)=" 
           << jet.true_btag() << endl;
      cout << "EE/HE=" << jet.EEoverHE()
           << " ntracks=" << jet.ntracks() << endl;
      cout << endl;
    }

    // Transverse missing energy (MET)
    cout << "MET pt=" << event.rec()->MET().pt()
         << " phi=" << event.rec()->MET().phi() << endl;
    cout << endl;

    // Transverse missing hadronic energy (MHT)
    cout << "MHT pt=" << event.rec()->MHT().pt()
              << " phi=" << event.rec()->MHT().phi() << endl;
    cout << endl;

    // Total transverse energy (TET)
    cout << "TET=" << event.rec()->TET() << endl;
    cout << endl;

    // Total transverse hadronic energy (THT)
    cout << "THT=" << event.rec()->THT() << endl;
    cout << endl;
  }
  */
  return true;
}

