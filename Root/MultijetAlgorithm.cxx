#include <EventLoop/Job.h> 
#include <EventLoop/Worker.h> 
#include "EventLoop/OutputStream.h" 
#include "AthContainers/ConstDataVector.h"

#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include <xAODJet/JetContainer.h>
#include <xAODTrigger/JetRoIContainer.h>
#include "xAODEventInfo/EventInfo.h"
#include <xAODAnaHelpers/MultijetAlgorithm.h>
#include <xAODAnaHelpers/HelperFunctions.h>
#include "AsgTools/MessageCheck.h"

#include "TFile.h"

#include "TSystem.h"

#include <iostream>
#include <fstream>

using namespace std;

// this is needed to distribute the algorithm to the workers
ClassImp(MultijetAlgorithm)

MultijetAlgorithm :: MultijetAlgorithm () :
  m_cutflowHist(0),
  m_cutflowHistW(0),
  m_treeStream("extradecorationstream"),
  m_trkSelTool(nullptr)
{}

EL::StatusCode MultijetAlgorithm :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  std::cout<< "Calling setupJob in MultijetAlgorithm."<<std::endl;
  job.useXAOD();
  xAOD::Init( "MultijetAlgorithm" ).ignore(); // call before opening first file

  EL::OutputStream outForTree( m_treeStream );
  job.outputAdd (outForTree);
  std::cout<<"in setupJob: all went well"<<std::endl;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MultijetAlgorithm :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  std::cout<< "Calling histInitialize"<<std::endl;
  ANA_CHECK( xAH::Algorithm::algInitialize());

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MultijetAlgorithm :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  std::cout<<"Calling fileExecute"<<std::endl;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MultijetAlgorithm :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  std::cout<< "Calling changeInput"<<std::endl;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MultijetAlgorithm :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  std::cout<< "Calling initialize"<<std::endl;

  m_eventCounter = 0;

  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  const xAOD::EventInfo* eventInfo(nullptr);
  ANA_CHECK( HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, msg()) );
  m_isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );


  // does not work on the grid
  TString fileName = TString(wk()->inputFile()->GetName());
  if( fileName.Contains("13TeV") ) {
    m_comEnergy = "13TeV";
  }
  else if( fileName.Contains("8TeV") ) {
    m_comEnergy = "8TeV";
  }
  else {
    std::cout << "No COM Energy could be found from the string name" << std::endl;
    std::cout << "Default to 13TeV cross section files" << std::endl;
    m_comEnergy = "13TeV";
  }
  std::cout << "Setting m_comEnergy to "<< m_comEnergy << std::endl;

   if(m_useCutFlow) {

    TFile *file = wk()->getOutputFile ("cutflow");
    m_cutflowHist  = (TH1D*)file->Get("cutflow");
    m_cutflowHistW = (TH1D*)file->Get("cutflow_weighted");

    m_cutflowFirst = m_cutflowHist->GetXaxis()->FindBin("okWeights");
    m_cutflowHistW->GetXaxis()->FindBin("okWeights");

  }

  // set up tools

    if ( m_TrackSelectionType == "Loose" ) { 
    	ANA_CHECK( ASG_MAKE_ANA_TOOL(m_TVATool_handle, CP::LooseTrackVertexAssociationTool));
    }
    else if ( m_TrackSelectionType == "Tight" ) { 
    	ANA_CHECK( ASG_MAKE_ANA_TOOL(m_TVATool_handle, CP::TightTrackVertexAssociationTool));
    }
    else {
	ANA_MSG_ERROR("invalid TrackSelection Type : "<<m_TrackSelectionType<<". Set to either 'Loose' or 'Tight'.");
	return EL::StatusCode::FAILURE;
    }
 
    ANA_CHECK( m_TVATool_handle.retrieve());
    ANA_MSG_INFO("Retrieved tool: "<< m_TVATool_handle);
    ANA_MSG_DEBUG("Tool info: ");
    m_TVATool_handle->print();

   m_trkSelTool = new InDet::InDetTrackSelectionTool( "JetTrackSelection", m_TrackSelectionType.c_str() ); // Type: Loose or Tight
   m_trkSelTool->initialize();


  Info("initialize()", "Succesfully initialized! \n");
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MultijetAlgorithm :: execute ()
{
  // Here you do everything that needs to be done on every single
  // event, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  if( m_msgLevel <= MSG::DEBUG) Info("execute()", "Applying selection \n");
  ++m_eventCounter;

  // retrieve event
  const xAOD::EventInfo* eventInfo(nullptr);
  ANA_CHECK( HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, msg()) );

  m_iCutflow = m_cutflowFirst;

  bool pass(false);

  const xAOD::JetContainer *inputJets = nullptr;
  ANA_CHECK( HelperFunctions::retrieve(inputJets, m_inContainerName, m_event, m_store, msg()) ); // for pile-up check.
 
  pass = true;
  passCut();

  
  // do truth matching for isHS, isPU labelling. 
  if ( m_isMC ) {
  const xAOD::JetContainer *truthJets = nullptr;
  if ( m_isMC ) ANA_CHECK( HelperFunctions::retrieve(truthJets, "AntiKt4TruthJets", m_event, m_store, msg()) ); // for pile-up check.
  // deal with large weighting due to pile-up overlay. I think this is noticeable only for higher slices, so maybe limit to just those?
  if( m_isMC && eventInfo->mcChannelNumber() > 361025){ 
      if(inputJets->size() > 1){
      float pTAvg = ( inputJets->at(0)->pt() + inputJets->at(1)->pt() ) /2.0;
      if( truthJets->size() > 0 && truthJets->at(0)->pt() > 0.0 ) { 
	float ratio = pTAvg/truthJets->at(0)->pt();
        if ( ratio < 0.6 || ratio > 1.4 ) {
      	wk()->skipEvent();  return EL::StatusCode::SUCCESS;
      	}
     }
    }
  }
  for ( auto jet : *inputJets ) {

	bool found_truth_match(false);
	bool is_isolated(true);
	unsigned int closeby_truthjet_count(0);
	
	for ( auto truthjet : *truthJets ) {
		std::cout<<"truthjet->pt() = "<<truthjet->pt()<<std::endl;
		if ( isMatch( jet, truthjet, 0.3 ) && truthjet->pt()/1000. >= 10. ) found_truth_match = true;	
		if ( isMatch( jet, truthjet, 0.6 ) && truthjet->pt()/1000. > 7. ) {
			closeby_truthjet_count+=1;
			is_isolated = false;
		} 
	}

	if ( found_truth_match && closeby_truthjet_count == 1 ) {
		jet->auxdecor< bool >("isHS") = true;
		jet->auxdecor< bool >("isPU") = false;
	}
	else if ( is_isolated ) {
		jet->auxdecor< bool >("isPU") = true;
		jet->auxdecor< bool >("isHS") = false;
	}
	else {
		jet->auxdecor< bool >("isPU") = false;
		jet->auxdecor< bool >("isHS") = false;
	}
  
  }
  }

  // MC event weight
  if ( m_isMC ) {
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    static SG::AuxElement::Accessor< float > mcEvtWeightAcc("mcEventWeight");
    if ( ! mcEvtWeightAcc.isAvailable( *eventInfo ) ) {
      ANA_MSG_ERROR( "mcEventWeight is not available as decoration! Aborting" );
      return EL::StatusCode::FAILURE;
    }
    m_mcEvtWeight = mcEvtWeightAcc( *eventInfo );
  }
  else m_mcEvtWeight = 1.;

  float weight_xs(1.);

  getDSWeights(eventInfo, weight_xs);

  eventInfo->auxdecor< float >("weight_xs") = weight_xs;
  eventInfo->auxdecor< float >("weight") = m_mcEvtWeight * weight_xs;
  }

  if ( m_doCompleteJVT ) {

     for ( auto jet : *inputJets ) {


	  const xAOD::VertexContainer* vertices(nullptr);
	    ANA_CHECK( HelperFunctions::retrieve(vertices, m_vertexContainerName, m_event, m_store, msg()) );
	  const xAOD::TrackParticleContainer* tracks(nullptr);
	    ANA_CHECK( HelperFunctions::retrieve(tracks, m_trackContainerName, m_event, m_store, msg()) );

	std::vector<float> JVFCorr_allPV;

	std::vector<float> SumPtTrkPt500 = jet->getAttribute<std::vector<float> >("SumPtTrkPt500");
	
	if ( SumPtTrkPt500.empty() ) {
		jet->auxdecor< std::vector<float> >("JVFCorrAllPV") = JVFCorr_allPV;
		continue;
	}

	for ( unsigned int PVindex = 0; PVindex < vertices->size(); PVindex++ ) {

		if ( PVindex >= SumPtTrkPt500.size() ) continue;

		float SumPtTrkPt500PV = SumPtTrkPt500.at(PVindex); 

		float sumTracknotPV = 0.;
	
		std::vector<const xAOD::TrackParticle*> ghostTracks;
		jet->getAssociatedObjects("GhostTrack", ghostTracks);
		for ( auto gtrklink: ghostTracks ){
			//const xAOD::TrackParticle* gtrk = ghostTracks.at(igTrk);
			const xAOD::TrackParticle* gtrk = dynamic_cast<const xAOD::TrackParticle*>(gtrklink);
			if( !m_trkSelTool->accept(*gtrk,0) ) continue;
			ElementLink< xAOD::VertexContainer> matchedVtxLink = m_TVATool_handle->getUniqueMatchVertexLink( *gtrk, *vertices );
			if ( !matchedVtxLink.isValid() ) continue;
			if ( (*matchedVtxLink)->index() == PVindex ) continue;
			sumTracknotPV+= gtrk->pt();
		}
		//float SumPtTrkPt500_nonPV = 0.;
		//for ( auto sumpttrkpt : SumPtTrkPt500 ) {
		//	SumPtTrkPt500_nonPV+= sumpttrkpt;
		//}
		//SumPtTrkPt500_nonPV-= SumPtTrkPt500PV;
		float jvfcorr(-999.);
	
		if ( SumPtTrkPt500PV + sumTracknotPV  <= 0 ) jvfcorr = -1.;
		else {
		
			xAOD::TrackVertexAssociationMap trktovxmap =m_TVATool_handle->getUniqueMatchMap( *tracks, *vertices );
		
			// count pile-up tracks for JVFCorr correction
			int n_pileuptrackcount = 0;
			for ( size_t iVertex = 0; iVertex < vertices->size(); iVertex++ ) {
				const xAOD::Vertex * vtx = vertices->at(iVertex);
				if ( vtx->index() == PVindex) continue;
				for(size_t iTrack = 0; iTrack < trktovxmap[vtx].size(); ++iTrack) {
					const xAOD::TrackParticle * track = trktovxmap[vtx].at(iTrack);
					if( !m_trkSelTool->accept(*track,0) ) continue;
					if ( (track->pt() < 30000. ) ) ++n_pileuptrackcount;
				 } // end of loop over tracks per vertex
			} // end of loop over vertices
		
			// first test: can recompute JVFCorr for HS?
			//if ( SumPtTrkPt500PV + SumPtTrkPt500_nonPV > 0 ) jvfcorr = SumPtTrkPt500PV/ ( SumPtTrkPt500PV +  (SumPtTrkPt500_nonPV / (0.01 * std::max(n_pileuptrackcount, 1 ) ) ));
			jvfcorr = SumPtTrkPt500PV/ ( SumPtTrkPt500PV +  (sumTracknotPV / (0.01 * std::max(n_pileuptrackcount, 1 ) ) ));
			
	
		} // *if* SumPts are not ZERO.

		JVFCorr_allPV.push_back( jvfcorr );
	
	
		//if ( SumPtTrkPt500PV > 0 ) {
		//	float rederivedCorr = SumPtTrkPt500_nonPV/(SumPtTrkPt500PV/realJVFCorr - SumPtTrkPt500PV);
		//	ANA_MSG_INFO("the rederived corr = "<<rederivedCorr );
		//}

	} // end of loop over vertices
	ANA_MSG_DEBUG("computed JVFCorr ="<<JVFCorr_allPV.at(0));
	float realJVFCorr = jet->getAttribute< float >("JVFCorr");
	ANA_MSG_DEBUG("real JVFCorr ="<<realJVFCorr);
        
    	jet->auxdecor< std::vector<float> >("JVFCorrAllPV") = JVFCorr_allPV;
   } // end of loop over jets	


  }

 // !! don't want to cut on pile-up jets just label them.
 // will use simple geomtreic matching reco <-> truth jets.

/*
  if ( m_isMC ) {

	if ( !m_truthLevelOnly ) {
		unsigned int matchedrecocnt(0);
		unsigned int recocnt(0);

		for ( auto jet : *inJets ) {
			recocnt++;
			unsigned int matchedtruthcnt = doPUchecks( jet, truthJets ); // returns # of truth jets matched to this jet for the moment.
			if ( matchedtruthcnt != 0 ) matchedrecocnt++;
			if (  m_msgLevel <= MSG::DEBUG) std::cout<< "no. of matched truth jets to reco jet # "<<recocnt<<": "<<matchedtruthcnt<<std::endl; 
		}	
		if ( m_msgLevel <= MSG::DEBUG) std::cout<<"total # of matched reco jet count for event# "<<m_eventCounter<<": "<<matchedrecocnt<<std::endl;
  	}


	if ( m_truthLevelOnly ) {

  		const xAOD::TruthParticleContainer *truthParticles = nullptr;
  		if ( m_isMC && m_truthLevelOnly ) ANA_CHECK( HelperFunctions::retrieve(truthParticles, "TruthBSM", m_event, m_store, msg()) ); // attempt to associate to original BSM particles.

		std::cout<<"size of truthParticle container = "<< truthParticles->size()<<std::endl;

		setOriginPDGID( inJets, truthParticles );
	}

  }	*/

  // look what we have in TStore
  if(msgLvl(MSG::VERBOSE)) m_store->print();

  if ( !pass ) {
    wk()->skipEvent();
  }

  std::cout<< "Leaving MultijetAlgo... "<<std::endl;

  return EL::StatusCode::SUCCESS;


}

// grab cross-section, k-factor and filter efficiency from getMetadata-generated text file.

EL::StatusCode MultijetAlgorithm::getDSWeights(const xAOD::EventInfo* eventInfo, float& weight_xs) {

  float xs(1.);
  float filtEff(1.);
  float kfac(1.);

  if(m_isMC){
  
    if(eventInfo->mcChannelNumber()==0) m_runNumber = eventInfo->runNumber();
    else m_runNumber = eventInfo->mcChannelNumber();

    std::string datafilename = PathResolverFindCalibFile("xAODAnaHelpers/"+m_metadata_filename);

    TTree t; t.ReadFile((datafilename).c_str()); // use PathResolver Here? 
    //TTree t; t.ReadFile(("/data/barn01/cantel/ATLAS/analysis/TLA/ntuplemaker/src/xAODAnaHelpers/data/"+m_metadata_filename).c_str()); // use PathResolver Here? 
	
    int channel=0; double crossSection=-999; double kFactor=-999; double genFiltEff=-999; 
    t.SetBranchAddress("dataset_number",&channel);
    t.SetBranchAddress("crossSection",&crossSection);
    t.SetBranchAddress("kFactor",&kFactor);
    t.SetBranchAddress("genFiltEff",&genFiltEff);
 

    for(int i=0;i<t.GetEntries();i++) {
  
      t.GetEntry(i);
      if (channel != m_runNumber) continue;
            xs = crossSection;
            kfac = kFactor;
      filtEff = genFiltEff;     
    }
  }

  std::cout<<"in getLumiWeights: "<<std::endl;
  std::cout<<"\t is MC? "<<m_isMC<<std::endl;
  std::cout<<"\t channel # = "<<m_runNumber<<std::endl;
  std::cout<<"\t\t xSec = "<<xs<<" "<<"kfac = "<<kfac<<" "<<"filtEff = "<<filtEff<<std::endl;

  weight_xs = xs * filtEff * kfac;

  return EL::StatusCode::SUCCESS;

}


EL::StatusCode MultijetAlgorithm :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.

  std::cout<<"Calling postExecute"<<std::endl;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MultijetAlgorithm :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  std::cout<<m_name<<std::endl;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MultijetAlgorithm :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.

  std::cout<<"Calling histFinalize"<<std::endl;
  
  if( m_writeTree ) {
    std::string thisName;
    m_ss.str( std::string() );
    m_ss << m_runNumber;
    TFile * treeFile = wk()->getOutputFile( m_treeStream );
    if(m_useCutFlow) {
      TH1F* thisCutflowHist = (TH1F*) m_cutflowHist->Clone();
      thisName = thisCutflowHist->GetName();
      thisCutflowHist->SetName( (thisName+"_"+m_ss.str()).c_str() );
      thisCutflowHist->SetDirectory( treeFile );

      TH1F* thisCutflowHistW = (TH1F*) m_cutflowHistW->Clone();
      thisName = thisCutflowHistW->GetName();
      thisCutflowHistW->SetName( (thisName+"_"+m_ss.str()).c_str() );
      thisCutflowHistW->SetDirectory( treeFile );
    }
  }

  ANA_CHECK( xAH::Algorithm::algFinalize());
  return EL::StatusCode::SUCCESS;
}

//
////Easy method for automatically filling cutflow and incrementing counter
////
void MultijetAlgorithm::passCut(){
  m_cutflowHist ->Fill(m_iCutflow, 1);
  m_cutflowHistW->Fill(m_iCutflow, m_mcEvtWeight);
  m_iCutflow++;
}

bool MultijetAlgorithm::isMatch(const xAOD::Jet_v1 * j1, const xAOD::Jet_v1 * j2, float cut ) {

	float eta1 = j1->eta();
	float phi1 = j1->phi();

	float eta2 = j2->eta();
	float phi2 = j2->phi();

	float dphi = TVector2::Phi_mpi_pi( phi1 - phi2 );

	float dR = sqrt((eta2-eta1)*(eta2-eta1) + (dphi)*(dphi));
	if (dR <= cut ) return true;	
	else return false;

}

unsigned int MultijetAlgorithm :: doPUchecks ( const xAOD::Jet * j, const xAOD::JetContainer * truthjCont  )
{

	TLorentzVector jvec;
	jvec.SetPtEtaPhiE( j->pt(), j->eta(), j->phi(), j->e() );

	bool isPU(true);
	double mindR(1e3);
	float mintruthpT(0);
	unsigned int matchedtruthcnt(0);
        float truthfrac(0);
        float truthres(-1);

	if ( m_truthLevelOnly ) isPU = false;
	else {
	
		for ( auto truthj : *truthjCont ) {
		
			TLorentzVector truthjvec;
	                truthjvec.SetPtEtaPhiE( truthj->pt(), truthj->eta(), truthj->phi(), truthj->e() );
		
			double dR = fabs(jvec.DeltaR( truthjvec ));
			if ( dR < 0.4 ) {
				isPU = false;
				matchedtruthcnt++;
				if ( dR< mindR ) { mindR = dR; mintruthpT = truthj->pt(); }
			} 		
		}
	}

	j->auxdecor< bool >("isPU") = isPU;
	j->auxdecor< unsigned int >("matchedtruthcnt") = matchedtruthcnt;
	
	if ( !isPU ) {
		truthfrac = j->pt()/mintruthpT;
		truthres = (j->pt() - mintruthpT )/mintruthpT;
	} 

        j->auxdecor< float >("truthfrac") = truthfrac;
        j->auxdecor< float >("truthres") = truthres;

	return matchedtruthcnt;
}

