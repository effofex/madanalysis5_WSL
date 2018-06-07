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
        [&min_pt,&maxAbsEta](const RecLeptonFormat& l){return (l.pt() > min_pt ) && (abs(l.eta())< maxAbsEta);});
}

void test_analysis::filterBaseLine(const vector<RecJetFormat> &from, vector<RecJetFormat> &to, 
                    double min_pt, double maxAbsEta){
  copy_if(from.begin(),from.end(),
	back_inserter(to),
        [&min_pt,&maxAbsEta](const RecJetFormat& l){return (l.pt() > min_pt) && (abs(l.eta())< maxAbsEta);});
}
// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool test_analysis::Execute(SampleFormat& sample, const EventFormat& event)
{
  if(0 != event.rec()){    
    vector<RecLeptonFormat> elecs;
    vector<RecLeptonFormat> muons;
    vector<RecJetFormat> jets;    
    int blSizeElec = 0;
    int blSizeMuon = 0;
    int blSizeJet  = 0;

    /***************************************************
    * Create baseline vectors
    ***************************************************/
    filterBaseLine(event.rec()->electrons(),elecs,5,2.47);
    filterBaseLine(event.rec()->muons(),muons,4,2.7);
    filterBaseLine(event.rec()->jets(),jets,25,2.5);
 

    /***********************************
    * Determine overlapping elements to remove
    **********************************/
      
    //Before doing any overlap removal, record the original size of the 
    //baseline vectors, so that they can be reported later
    blSizeJet = jets.size();
    blSizeElec = elecs.size();
    blSizeMuon = muons.size();
    
    //Overlap removal between electrons and jets, partial implementation of col
    //#2 in Table 3 of CERN-ERP-2017-246 https://arxiv.org/pdf/1711.11520.pdf
     for(MAuint32 j=0;j<elecs.size();j++)
    {
      const RecLeptonFormat& myelec = elecs[j];
      jets.erase(remove_if(jets.begin(),jets.end(),
        [&myelec](RecJetFormat& j){return j.dr(myelec);}),
        jets.end());
    }
    
    //Overlap removal between jets and leptons, partialimplmentation of col 
    //#4 in Table 3 of CERN-ERP-2017-246 https://arxiv.org/pdf/1711.11520.pdf
    //Performed in two steps, since both elecs and muons are leptons
    for(MAuint32 i=0;i<jets.size();i++)
    {
      const RecJetFormat& myjet = jets[i];
      elecs.erase(remove_if(elecs.begin(),elecs.end(),
        [&myjet](RecLeptonFormat& l){return myjet.dr(l)<min(0.4,0.4+10/l.pt());}),
        elecs.end());
    }
         
    //Repeating the removal for muons.  The duplication here is a bit of code
    //itch, but I'm holding off on extracting a function until we do the 
    //full overlap algorithm
    for(MAuint32 i=0;i<jets.size();i++)
    {
      const RecJetFormat& myjet = jets[i];
      muons.erase(remove_if(muons.begin(),muons.end(),
        [&myjet](RecLeptonFormat& l){return myjet.dr(l)<min(0.4,0.4+10/l.pt());}),
        muons.end());
    }
 
   
   /*******************
   * Output the results
   ********************/ 
    cout << "Base and signal vector sizes" << endl;
    cout << "Particle\tBase\tRemoved\tSignal" << endl;
    cout << "\tJets:\t" << blSizeJet << "\t" << blSizeJet-jets.size() << "\t" << jets.size() << endl;
    cout << "\tElecs:\t" << blSizeElec << "\t" << blSizeElec-elecs.size() << "\t" << elecs.size() << endl;
    cout << "\tMuons:\t" << blSizeMuon << "\t" << blSizeMuon-muons.size() << "\t" << muons.size() << endl;

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

