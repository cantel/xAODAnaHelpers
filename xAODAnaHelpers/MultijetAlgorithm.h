#ifndef MultijetAlgo_MultijetAlgorithm_H
#define MultijetAlgo_MultijetAlgorithm_H

#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>

//algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

#ifndef __CINT__
  #include "xAODJet/JetContainer.h"
  #include "xAODTruth/TruthParticleContainer.h"
  #include "xAODEventInfo/EventInfo.h"
  #include <xAODAnaHelpers/JetHists.h>
#endif

// ROOT include(s):
#include "TH1D.h"

// external includes
#include "TrackVertexAssociationTool/LooseTrackVertexAssociationTool.h"
#include "TrackVertexAssociationTool/TightTrackVertexAssociationTool.h"
#include "TrackVertexAssociationTool/BaseTrackVertexAssociationTool.h"
#include "TrackVertexAssociationTool/ITrackVertexAssociationTool.h"
#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"

#include <sstream>

static float GeV = 1000.;

class MultijetAlgorithm : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  public:
    bool m_writeTree = false;                     // true will write out a TTree
    std::string m_metadata_filename = ""; // metadata file containg cross-sections, k-factors etc - see README on details to obtain this file.
    std::string m_inContainerName = "";        // input container name
    std::string m_secondJetContainerName = "AntiKt4TruthJets"; // Name of jet container for resolution reference
    std::string m_vertexContainerName = "PrimaryVertices";
    std::string m_trackContainerName = "InDetTrackParticles";
    bool m_useCutFlow = true;
    bool m_truthLevelOnly = false;                // truthLevelOnly info
    std::string m_inputAlgo = "";              // input algo for when running systs
    MSG::Level m_msgLevel = MSG::INFO;
    //bool m_debug = true;                         // set verbose mode
    bool m_isMC = false;                          // Is MC
    bool m_doCompleteJVT = false;
    // float cutValue;
    xAOD::TEvent *m_event;  //!
    xAOD::TStore *m_store;  //!
    unsigned int m_eventCounter;     //!
    std::string m_treeStream;
    std::string m_TrackSelectionType = "Loose"; // choose between Loose or Tight


    std::string m_name;
    TH1D* m_cutflowHist;    //!
    TH1D* m_cutflowHistW;   //!
    int m_cutflowFirst;     //!
    int m_iCutflow;         //!
    float m_mcEvtWeight;  //!
    std::string m_comEnergy; //!

  private:
    //configuration variables
    int m_runNumber; //!
    std::stringstream m_ss; //!

    void passCut();

    asg::AnaToolHandle<CP::ITrackVertexAssociationTool> m_TVATool_handle { "TVATool", this }; //!
    InDet::InDetTrackSelectionTool * m_trkSelTool;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!




  // this is a standard constructor
  MultijetAlgorithm ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // these are the functions not inherited from Algorithm
#ifndef __CINT__
  bool executeAnalysis( const xAOD::EventInfo* eventInfo, 
      const xAOD::JetContainer* signalJets, 
      bool count);
#endif // not __CINT__
  void AddTree( std::string );
  virtual EL::StatusCode getDSWeights( const xAOD::EventInfo* eventInfo, float& weight_xs );
  unsigned int doPUchecks ( const xAOD::Jet * j, const xAOD::JetContainer * truthjCont  );
  bool isMatch(const xAOD::Jet_v1 * j1, const xAOD::Jet_v1 * j2, float cut );
  // this is needed to distribute the algorithm to the workers
  ClassDef(MultijetAlgorithm, 1);
};

#endif
