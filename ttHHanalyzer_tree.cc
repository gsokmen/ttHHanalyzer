//---------------------------------------------------------------------------
// File:        ttHHanalyzer.cc
// Description: Analyzer for simple ntuples, such as those created by
//              TheNtupleMaker
// Created:     Tue Dec  1 00:49:09 2020 by mkanalyzer.py v2.0.2 15-Apr-2019
// Author:      Shakespeare's ghost
//----------------------------------------------------------------------------
#include "tnm.h"
#include <cmath> 
#include <algorithm>
#include <TLorentzVector.h>
#include <TTree.h>

using namespace std;
//----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // If you want canvases to be visible during program execution, just
  // uncomment the line below
  //TApplication app("ttHHanalyzer", &argc, argv);

  // Get command line arguments
  commandLine cl(argc, argv);
    
  // Get names of ntuple files to be processed
  vector<string> filenames = fileNames(cl.filelist);

  // Get global weight 
  double weight = cl.externalweight;

  // Create tree reader
  itreestream stream(filenames, "Delphes");
  if ( !stream.good() ) error("can't read root input files");

  // Create a buffer to receive events from the stream
  // The default is to select all branches
  // Use second argument to select specific branches
  // Example:
  //   varlist = 'Jet_PT Jet_Eta Jet_Phi'
  //   ev = eventBuffer(stream, varlist)

  eventBuffer ev(stream);
  
  int nevents = ev.size();
  cout << "number of events: " << nevents << endl;

  // Create output file for histograms; see notes in header 
  outputFile of(cl.outputfilename);

  // -------------------------------------------------------------------------
  // Define histograms
  // -------------------------------------------------------------------------
  //setStyle();


  of.file->cd();
  TDirectory *jet = of.file->mkdir("jet");
  jet->cd();

  /*  TH1F * hmet = new TH1F("met", "MET", 50, 0, 1000);
  TH1F * hmetPhi = new TH1F("metPhi", "MET #phi", 50, -5, 5);
  TH1F * hmetEta = new TH1F("metEta", "MET #eta", 50, -5, 5);
  TH1F * hjetPT1 = new TH1F("jetPT1", "jet1 p_{T}", 50, 0, 1000);
  TH1F * hjetPT2 = new TH1F("jetPT2", "jet2 p_{T}", 50, 0, 1000);
  TH1F * hjetPT3 = new TH1F("jetPT3", "jet3 p_{T}", 50, 0, 1000);
  TH1F * hjetPT4 = new TH1F("jetPT4", "jet4 p_{T}", 50, 0, 1000); */
  TH1F * hjetPT5 = new TH1F("jetPT5", "jet5 p_{T}", 50, 0, 1000);
  TH1F * hjetPT6 = new TH1F("jetPT6", "jet6 p_{T}", 50, 0, 1000);
  /*  TH1F * hjetEta1 = new TH1F("jetEta1", "jet1 #eta", 50, -5, 5);
  TH1F * hjetEta2 = new TH1F("jetEta2", "jet2 #eta", 50, -5, 5);
  TH1F * hjetEta3 = new TH1F("jetEta3", "jet3 #eta", 50, -5, 5);
  TH1F * hjetEta4 = new TH1F("jetEta4", "jet4 #eta", 50, -5, 5); */
  TH1F * hjetEta5 = new TH1F("jetEta5", "jet5 #eta", 50, -5, 5);
  TH1F * hjetEta6 = new TH1F("jetEta6", "jet6 #eta", 50, -5, 5);
  /*  TH1F * hBjetPT1 = new TH1F("jetBPT1", "bjet1 p_{T}", 50, 0, 1000);
  TH1F * hBjetPT2 = new TH1F("jetBPT2", "bjet2 p_{T}", 50, 0, 1000);
  TH1F * hBjetPT3 = new TH1F("jetBPT3", "bjet3 p_{T}", 50, 0, 1000);
  TH1F * hBjetPT4 = new TH1F("jetBPT4", "bjet4 p_{T}", 50, 0, 1000); */
  TH1F * hBjetPT5 = new TH1F("jetBPT5", "bjet5 p_{T}", 50, 0, 1000);
  TH1F * hBjetPT6 = new TH1F("jetBPT6", "bjet6 p_{T}", 50, 0, 1000);
  /*  TH1F * hBjetEta1 = new TH1F("jetBEta1", "bjet1 #eta", 50, -5, 5);
  TH1F * hBjetEta2 = new TH1F("jetBEta2", "bjet2 #eta", 50, -5, 5);
  TH1F * hBjetEta3 = new TH1F("jetBEta3", "bjet3 #eta", 50, -5, 5);
  TH1F * hBjetEta4 = new TH1F("jetBEta4", "bjet4 #eta", 50, -5, 5); */
  TH1F * hBjetEta5 = new TH1F("jetBEta5", "bjet5 #eta", 50, -5, 5);
  TH1F * hBjetEta6 = new TH1F("jetBEta6", "bjet6 #eta", 50, -5, 5);
  /*  TH1F * hAvgDeltaRjj = new TH1F("deltaRavgjj", "#DeltaR_{jj}^{avg}", 50, 0, 5);
  TH1F * hAvgDeltaRbb = new TH1F("deltaRavgbb", "#DeltaR_{bb}^{avg}", 50, 0, 5);
  TH1F * hAvgDeltaEtajj = new TH1F("deltaEtaavgjj", "#DeltaEta_{jj}^{avg}", 50, 0, 5);
  TH1F * hAvgDeltaEtabb = new TH1F("deltaEtaavgbb", "#DeltaEta_{bb}^{avg}", 50, 0, 5);
  TH1F * hminDeltaRjj = new TH1F("deltaRminjj", "#DeltaR_{jj}^{min}", 50, 0, 5);
  TH1F * hminDeltaRbb = new TH1F("deltaRminbb", "#DeltaR_{bb}^{min}", 50, 0, 5);
  TH1F * hminDeltaRlj = new TH1F("deltaRminlj", "#DeltaR_{lep,j}^{min}", 50, 0, 5);
  TH1F * hminDeltaRlb = new TH1F("deltaRminlb", "#DeltaR_{lep,b}^{min}", 50, 0, 5);
  TH1F * hmaxDeltaEtabb = new TH1F("deltaEtamaxbb", "#DeltaEta_{bb}^{max}", 50, 0, 5);
  TH1F * hjetAverageMass = new TH1F("jetAvgMass", "m_{j}^{avg}", 50, 0, 100);
  TH1F * hBjetAverageMass = new TH1F("jetBAvgMass", "m_{b}^{avg}", 50, 0, 100);
  TH1F * hBjetAverageMassSqr = new TH1F("jetBAvgMassSqr", "(m^{2})_{b}^{avg}", 50, 0, 1000);
  TH1F * hjetHT = new TH1F("jetHT", "H_{T}", 50, 0, 3000);
  TH1F * hBjetHT = new TH1F("jetBHT", "H_{T}^{b}", 50, 0, 3000);
  TH1F * hjetNumber = new TH1F("jetNumber", "N_{jet}", 20, 0, 20);
  TH1F * hBjetNumber = new TH1F("jetBNumber", "N_{bjet}", 20, 0, 20);
  TH1F * hChi2 = new TH1F("chi2", "#chi^{2}", 100, 0, 3000); */
  

  of.file->cd();
  TDirectory *lepton = of.file->mkdir("Lepton");
  lepton->cd();

  TH1F * hLeptonPT = new TH1F("leptonPT", "lepton p_{T}", 50, 0, 1000);
  TH1F * hMuonPT = new TH1F("muonPT", "muon p_{T}", 50, 0, 1000);
  TH1F * hElePT = new TH1F("elePT", "ele p_{T}", 50, 0, 1000);
  TH1F * hLeptonPhi = new TH1F("leptonPhi", "lepton #phi", 50, 0, 1000);
  TH1F * hMuonPhi = new TH1F("muonPhi", "muon #phi", 50, -5, 5);
  TH1F * hElePhi = new TH1F("elePhi", "ele #phi", 50, -5, 5);
  TH1F * hLeptonEta = new TH1F("leptonEta", "lepton #eta", 50, -5, 5);
  TH1F * hMuonEta = new TH1F("muonEta", "muon #eta", 50, -5, 5);
  TH1F * hEleEta = new TH1F("eleEta", "ele #eta", 50, -5, 5);

  // -------------------------------------------------------------------------
  // Define Trees
  // -------------------------------------------------------------------------

  TFile hfile("treeOnly.root","RECREATE","Demo ROOT file with trees");

  TTree tree("T","An example of ROOT tree with a few branches");

  // -------------------------------------------------------------------------
  // Loop over events
  // -------------------------------------------------------------------------
  

  for(int entry = 0; entry < nevents; entry++){ 
      ev.read(entry);  // read an event into event buffer

      // Uncomment the following line if you wish to copy variables into structs. See the file eventBuffer.h to find out what structs are available. Alternatively, you can call individual fill functions,such as, ev.fillElectrons().

      ev.fillObjects();

      // Define shorthand names  // analysis   // Loop over particles:

      std::vector<eventBuffer::Jet_s> jet = ev.Jet;
      std::vector<eventBuffer::Jet_s> selJet;
      std::vector<eventBuffer::Jet_s> selbJet;
      std::vector<eventBuffer::MuonTight_s> muonT = ev.MuonTight;
      std::vector<eventBuffer::MuonTight_s> selMuonT;
      std::vector<eventBuffer::Electron_s> ele = ev.Electron;
      std::vector<eventBuffer::Electron_s> selEle;
      std::vector<eventBuffer::Electron_s> selLept;
      float met = ev.MissingET_MET;
      float metEta = ev.MissingET_Eta;            
      float metPhi = ev.MissingET_Phi;

      //weight = ev.Weight[0].Weight) 
      //float weight = 0.3954; // ttHtobb with 3846848 events
      //float weight = 11.28; //ttbarInc with 1243317 events !! 
      //float weight = 0.2321; // ttZ with 9243300 events 
      //float weight = 0.0113; // ttHH with 193300 events
      //float weight = 0.0092; //ttZZ with 499986 events      // DON'T FORGET TO CHANGE THE WEIGHT // 
      //float weight = 0.007422; //ttZH with 500000 events
      
      float sumJetMass = 0, sumBjetMass = 0, scalarBjetSumPT = 0, avgBjetMass =0,  avgBjetMassSqr = 0, avgJetMass =  0, scalarSumPT = 0;
      float deltaR = 0, deltaPhi = 0, deltaEta = 0, sqrDeltaPhi = 0, sqrDeltaEta = 0, sumDeltaR = 0, sumDeltaEta = 0, avgDeltaR = 0, avgDeltaEta = 0;
      float deltaR_bb = 0, deltaPhi_bb = 0, deltaEta_bb = 0, sqrDeltaPhi_bb = 0, sqrDeltaEta_bb = 0, sumDeltaR_bb = 0, sumDeltaEta_bb = 0, avgDeltaR_bb = 0, avgDeltaEta_bb = 0;
      float deltaR_lj = 0, deltaPhi_lj = 0, deltaEta_lj = 0, sqrDeltaPhi_lj = 0, sqrDeltaEta_lj = 0;
      float deltaR_lb = 0, deltaPhi_lb = 0, deltaEta_lb = 0, sqrDeltaPhi_lb = 0, sqrDeltaEta_lb = 0;
      float minDeltaRjj = 1000000., minDeltaRbb = 1000000., minDeltaRlb = 1000000., minDeltaRlj = 1000000.;
      float maxDeltaEtabb = 0;

      int nJets = 0, nBJets = 0, nMuons = 0, nEle = 0, nNonVetoMuons = 0, nNonVetoEle = 0;
      int nCombJ = 0, nCombB = 0 ;

      float minChi2 = 1000000.;
      float higgsMass = 125.;
      float chi2 = 0;
      TLorentzVector pBB1 , pBB2 ;
      //      double_t resolution = 0.15;

      TLorentzVector bjet1, bjet2, bjet3, bjet4;

      //      tree.Branch("scalarSumPT", &scalarSumPT, "scalarSumPT/f");
      // tree.Branch("scalarBjetSumPT", &scalarBjetSumPT, "scalarBjetSumPT/f");


      for(int i=0; i < jet.size(); i++){
	if(jet[i].PT > 30 && fabs(jet[i].Eta) < 2.4){
	  scalarSumPT += fabs(jet[i].PT);
	  sumJetMass += jet[i].Mass;
	  for( int j = i+1; j < jet.size(); j++){
	    if(jet[j].PT > 30 && fabs(jet[j].Eta) < 2.4){
	      deltaPhi = fabs(jet[i].Phi - jet[j].Phi);
	      deltaEta = fabs(jet[i].Eta - jet[j].Eta);
	      sqrDeltaPhi = deltaPhi * deltaPhi;
	      sqrDeltaEta = deltaEta * deltaEta;
	      deltaR = sqrt(sqrDeltaPhi + sqrDeltaEta);  
	      if( minDeltaRjj > deltaR){
		minDeltaRjj = deltaR;
	      }
	      sumDeltaEta += deltaEta;
	      sumDeltaR += deltaR;
	      nCombJ++;
	    }
	  }
	  nJets++;
	  selJet.push_back(jet[i]);
	  // cout <<  "count: " << i << endl; 
	  if(jet[i].BTag == 1){
	    scalarBjetSumPT += fabs(jet[i].PT);
	    sumBjetMass +=jet[i].Mass;
	    for( int j = i+1; j < jet.size(); j++){	    
	      if(jet[j].PT > 30 && fabs(jet[j].Eta) < 2.4){
		if(jet[j].BTag == 1){
 		  deltaPhi_bb = fabs(jet[i].Phi - jet[j].Phi);
		  deltaEta_bb = fabs(jet[i].Eta - jet[j].Eta);
		  sqrDeltaPhi_bb = deltaPhi_bb * deltaPhi_bb;
		  sqrDeltaEta_bb = deltaEta_bb * deltaEta_bb;
		  deltaR_bb = sqrt(sqrDeltaPhi_bb + sqrDeltaEta_bb);  
		  if( minDeltaRbb > deltaR_bb){
		    minDeltaRbb = deltaR_bb;
		  }
		  if( maxDeltaEtabb < deltaEta_bb){
		    maxDeltaEtabb = deltaEta_bb;
		  }
		  sumDeltaEta_bb += deltaEta_bb;
		  sumDeltaR_bb += deltaR_bb;
		  nCombB++;
		}
	      }
	    }
	    nBJets++;
	    selbJet.push_back(jet[i]);
	  }
	}
      }
      

      avgDeltaR = sumDeltaR/(double)nCombJ;
      avgDeltaR_bb = sumDeltaR_bb/(double)nCombB;
      avgDeltaEta = sumDeltaEta/(double)nCombJ;
      avgDeltaEta_bb = sumDeltaEta_bb/(double)nCombB;
      avgJetMass = sumJetMass/nJets;      
      avgBjetMass = sumBjetMass/nBJets;
      avgBjetMassSqr = (avgBjetMass)*(avgBjetMass);
      

      //
      if( nBJets > 3){
	for(int i=0; i < selbJet.size(); i++){
	  for (int j = i+1; j < selbJet.size(); j++){
	    if( i != j){
	      for(int k = 0; k < selbJet.size(); k++){
		if(j != k && j != i){
		  for(int m = k+1; m < selbJet.size(); m++){
		    if(m != k && m != j && m != i){
		      bjet1.SetPtEtaPhiM(jet[i].PT, jet[i].Eta, jet[i].Phi, jet[i].Mass);
		      bjet2.SetPtEtaPhiM(jet[j].PT, jet[j].Eta, jet[j].Phi, jet[j].Mass);
		      bjet3.SetPtEtaPhiM(jet[k].PT, jet[k].Eta, jet[k].Phi, jet[k].Mass);
		      bjet4.SetPtEtaPhiM(jet[m].PT, jet[m].Eta, jet[m].Phi, jet[m].Mass);
		      pBB1 = bjet1 + bjet2;
		      pBB2 = bjet3 + bjet4;
		      chi2 = (pBB1.M() - higgsMass)*(pBB1.M() - higgsMass) + (pBB2.M() - higgsMass)*(pBB2.M() - higgsMass);
		      if(minChi2 > chi2){
			minChi2 = chi2;
		      }
		    }
		  } 
		}
	      }
	    }
	  }
	}
	//	cout << "minChi2: " << minChi2 << " jet count: " << nBJets << endl; 	
      }
      
 	    //      

      for(int i = 0; i < muonT.size(); i++){
	if(muonT[i].PT > 29 && fabs(muonT[i].Eta) < 2.4){	 
	  for( int j = 0; j < jet.size(); j++){	    
	    if(jet[j].PT > 30 && fabs(jet[j].Eta) < 2.4){
	      deltaPhi_lj = fabs(muonT[i].Phi - jet[j].Phi);
	      deltaEta_lj = fabs(muonT[i].Eta - jet[j].Eta);
	      sqrDeltaPhi_lj = deltaPhi_lj * deltaPhi_lj;
	      sqrDeltaEta_lj = deltaEta_lj * deltaEta_lj;
	      deltaR_lj = sqrt(sqrDeltaPhi_lj + sqrDeltaEta_lj);  
	      if(minDeltaRlj > deltaR_lj){
		minDeltaRlj = deltaR_lj;
	      } 
	      if(jet[j].BTag == 1){
		deltaPhi_lb = fabs(muonT[i].Phi - jet[j].Phi);
		deltaEta_lb = fabs(muonT[i].Eta - jet[j].Eta);
		sqrDeltaPhi_lb = deltaPhi_lb * deltaPhi_lb;
		sqrDeltaEta_lb = deltaEta_lb * deltaEta_lb;
		deltaR_lb = sqrt(sqrDeltaPhi_lb + sqrDeltaEta_lb);  
		if(minDeltaRlb > deltaR_lb){
		  minDeltaRlb = deltaR_lb;
		} 
	      }
	    }
	  }
	  nMuons++;
	  selMuonT.push_back(muonT[i]);
	}
	if(muonT[i].PT >= 15 && muonT[i].PT <= 29){	 
	  nNonVetoMuons++;
	}
      }

     
      for(int i = 0; i < ele.size(); i++){
	if(ele[i].PT > 30 && fabs(ele[i].Eta) < 2.4){	 
	  for( int j = 0; j < jet.size(); j++){	    
	    if(jet[j].PT > 30 && fabs(jet[j].Eta) < 2.4){
	      deltaPhi_lj = fabs(ele[i].Phi - jet[j].Phi);
	      deltaEta_lj = fabs(ele[i].Eta - jet[j].Eta);
	      sqrDeltaPhi_lj = deltaPhi_lj * deltaPhi_lj;
	      sqrDeltaEta_lj = deltaEta_lj * deltaEta_lj;
	      deltaR_lj = sqrt(sqrDeltaPhi_lj + sqrDeltaEta_lj);  
	      if(minDeltaRlj > deltaR_lj){
		minDeltaRlj = deltaR_lj;
	      } 
	      if(jet[j].BTag == 1){
		deltaPhi_lb = fabs(ele[i].Phi - jet[j].Phi);
		deltaEta_lb = fabs(ele[i].Eta - jet[j].Eta);
		sqrDeltaPhi_lb = deltaPhi_lb * deltaPhi_lb;
		sqrDeltaEta_lb = deltaEta_lb * deltaEta_lb;
		deltaR_lb = sqrt(sqrDeltaPhi_lb + sqrDeltaEta_lb);  
		if(minDeltaRlb > deltaR_lb){
		  minDeltaRlb = deltaR_lb;
		} 
	      }
	    }
	  }
	  nEle++;
	  selEle.push_back(ele[i]);
	}
	if(ele[i].PT >= 15 && ele[i].PT <= 30){	 
	  nNonVetoEle++;
	}
      }


      int nLeptons = 0;

      if(nNonVetoEle == 0 && nNonVetoMuons == 0){
	nLeptons = nMuons + nEle;		
      }
      
      float jetPT1, jetPT2, jetPT3, jetPT4, jetPT5, jetPT6, bjetPT1, bjetPT2, bjetPT3, bjetPT4, bjetPT5, bjetPT6, jetEta1, jetEta2, jetEta3, jetEta4, jetEta5, jetEta6, bjetEta1, bjetEta2, bjetEta3, bjetEta4, bjetEta5, bjetEta6;
      float leptonPT, leptonEta, muonPT, muonPhi, muonEta, elePT, elePhi, eleEta;

      //      std::cout<< "starting to initialize the trees" << std::endl;

      tree.Branch("jetPT1", &jetPT1, "jetPT1/f");      
      tree.Branch("jetPT2", &jetPT2, "jetPT2/f");	
      tree.Branch("jetPT3", &jetPT3, "jetPT3/f");	
      tree.Branch("jetPT4", &jetPT4, "jetPT4/f");	
      tree.Branch("jetPT5", &jetPT5, "jetPT5/f");	
      tree.Branch("jetPT6", &jetPT6, "jetPT6/f");	
      tree.Branch("bjetPT1", &bjetPT1, "bjetPT1/f");      
      tree.Branch("bjetPT2", &bjetPT2, "bjetPT2/f");	
      tree.Branch("bjetPT3", &bjetPT3, "bjetPT3/f");	
      tree.Branch("bjetPT4", &bjetPT4, "bjetPT4/f");	
      tree.Branch("bjetPT5", &bjetPT5, "bjetPT5/f");	
      tree.Branch("bjetPT6", &bjetPT6, "bjetPT6/f");	
      tree.Branch("jetEta1", &jetEta1, "jetEta1/f");
      tree.Branch("jetEta2", &jetEta2, "jetEta2/f");		
      tree.Branch("jetEta3", &jetEta3, "jetEta3/f");	
      tree.Branch("jetEta4", &jetEta4, "jetEta4/f");	
      tree.Branch("jetEta5", &jetEta5, "jetEta5/f");	
      tree.Branch("jetEta6", &jetEta6, "jetEta6/f");	
      tree.Branch("bjetEta1", &bjetEta1, "bjetEta1/f");	
      tree.Branch("bjetEta2", &bjetEta2, "bjetEta2/f");	
      tree.Branch("bjetEta3", &bjetEta3, "bjetEta3/f");	
      tree.Branch("bjetEta4", &bjetEta4, "bjetEta4/f");	
      tree.Branch("bjetEta5", &bjetEta5, "bjetEta5/f");	
      tree.Branch("bjetEta6", &bjetEta6, "bjetEta6/f");	
      tree.Branch("scalarSumPT", &scalarSumPT, "scalarSumPT/f");
      tree.Branch("scalarBjetSumPT", &scalarBjetSumPT, "scalarBjetSumPT/f");
      tree.Branch("met", &met, "met/f");      
      tree.Branch("metPhi", &metPhi, "metPhi/f");      
      tree.Branch("metEta", &metEta, "metEta/f");   
      tree.Branch("jetNumber", &nJets, "nJets/i");
      tree.Branch("bjetNumber", &nBJets, "nBJets/i");
      tree.Branch("avgDeltaR", &avgDeltaR, "avgDeltaR/f");
      tree.Branch("avgDeltaR_bb", &avgDeltaR_bb, "avgDeltaR_bb/f");
      tree.Branch("avgDeltaEta", &avgDeltaEta, "avgDeltaEta/f");
      tree.Branch("avgDeltaEta", &avgDeltaEta_bb, "avgDeltaEta_bb/f");
      tree.Branch("minDeltaRjj", &minDeltaRjj, "minDeltaRjj/f");
      tree.Branch("minDeltaRbb", &minDeltaRbb, "minDeltaRbb/f");
      tree.Branch("minDeltaRlj", &minDeltaRlj, "minDeltaRlj/f");
      tree.Branch("minDeltaRlb", &minDeltaRlb, "minDeltaRlb/f");
      tree.Branch("maxDeltaEtabb", &maxDeltaEtabb, "maxDeltaEtabb/f");
      tree.Branch("avgJetMass", &avgJetMass, "avgJetMass/f");
      tree.Branch("avgBjetMass", &avgBjetMass, "avgBjetMass/f");
      tree.Branch("avgBjetMassSqr", &avgBjetMassSqr, "avgBjetMassSqr/f");
      tree.Branch("minchi2", &minChi2, "minChi2/f");
      tree.Branch("leptonPT", &leptonPT, "leptonPT/f");
      tree.Branch("leptonEta", &leptonEta, "leptonEta/f"); 
      tree.Branch("muonPT", &muonPT, "muonPT/f");
      tree.Branch("muonPhi", &muonPhi, "muonPhi/f");
      tree.Branch("muonEta", &muonEta, "muonEta/f");
      tree.Branch("elePT", &elePT, "elePT/f");
      tree.Branch("elePhi", &elePhi, "elePhi/f");
      tree.Branch("eleEta", &eleEta, "eleEta/f"); 

      if(nJets >= 4 && nBJets >= 4  && nLeptons == 1 && met > 20){
	
	/*	hmet->Fill(met, weight);
	hmetEta->Fill(metEta, weight);
	hmetPhi->Fill(metPhi, weight);
	hjetPT1->Fill(selJet[0].PT, weight);
	hjetPT2->Fill(selJet[1].PT, weight);
	hjetPT3->Fill(selJet[2].PT, weight);
	hjetPT4->Fill(selJet[3].PT, weight);
	hjetEta1->Fill(selJet[0].Eta, weight);
	hjetEta2->Fill(selJet[1].Eta, weight);
	hjetEta3->Fill(selJet[2].Eta, weight);
	hjetEta4->Fill(selJet[3].Eta, weight);
	hBjetPT1->Fill(selbJet[0].PT, weight);
	hBjetPT2->Fill(selbJet[1].PT, weight);
	hBjetPT3->Fill(selbJet[2].PT, weight);
	hBjetPT4->Fill(selbJet[3].PT, weight);
	hBjetEta1->Fill(selbJet[0].Eta, weight);
	hBjetEta2->Fill(selbJet[1].Eta, weight);
	hBjetEta3->Fill(selbJet[2].Eta, weight);
	hBjetEta4->Fill(selbJet[3].Eta, weight); */

	jetPT1 = selJet[0].PT;
	jetPT2 = selJet[1].PT;
	jetPT3 = selJet[2].PT;
	jetPT4 = selJet[3].PT;
	bjetPT1 = selbJet[0].PT;
	bjetPT2 = selbJet[1].PT;
	bjetPT3 = selbJet[2].PT;
	bjetPT4 = selbJet[3].PT;
	jetEta1 = selJet[0].Eta;
	jetEta2 = selJet[1].Eta;
	jetEta3 = selJet[2].Eta;
	jetEta4 = selJet[3].Eta;		
	bjetEta1 = selbJet[0].Eta;
	bjetEta2 = selbJet[1].Eta;
	bjetEta3 = selbJet[2].Eta;     
	bjetEta4 = selbJet[3].Eta;     

	if(nJets > 4){
	  hjetPT5->Fill(selJet[4].PT, weight);
	  hjetEta5->Fill(selJet[4].Eta, weight);
	  jetPT5 = selJet[4].PT;
	  jetEta5 = selJet[4].Eta;
	} else{
	  jetPT5 = -99;
	  jetEta5 = -99;
	}

	if(nJets > 5){
	  hjetPT6->Fill(selJet[5].PT, weight);
	  hjetEta6->Fill(selJet[5].Eta, weight);
	  jetPT6 = selJet[5].PT;
	  jetEta6 = selJet[5].Eta;
	} else{
	  jetPT6 = -99;
	  jetEta6 = -99;
	}

	if(nBJets > 4) {
	  hBjetEta5->Fill(selbJet[4].Eta, weight);		
	  hBjetPT5->Fill(selbJet[4].PT, weight);
	  bjetPT5 = selbJet[4].PT;
	  bjetEta5 = selbJet[4].Eta;	
	} else {
	  bjetPT5 = -99;
	  bjetEta5 = -99;
	}

	if(nBJets > 5) {
	  hBjetEta6->Fill(selbJet[5].Eta, weight);		
	  hBjetPT6->Fill(selbJet[5].PT, weight);
	  bjetPT6 = selbJet[5].PT;
	  bjetEta6 = selbJet[5].Eta;	
	} else {
	  bjetPT6 = -99;
	  bjetEta6 = -99;
	}
	
	if(nMuons == 1){
	  hLeptonPT->Fill(selMuonT[0].PT, weight);
	  hLeptonPhi->Fill(selMuonT[0].Phi, weight);
	  hLeptonEta->Fill(selMuonT[0].Eta, weight);
	  leptonPT = selMuonT[0].PT;
	  leptonEta = selMuonT[0].Eta;
	}

	else if(nEle == 1){
	  hLeptonPT->Fill(selEle[0].PT, weight);
	  hLeptonPhi->Fill(selEle[0].Phi, weight);
	  hLeptonEta->Fill(selEle[0].Eta, weight);
	  leptonPT = selEle[0].PT;
	  leptonEta = selEle[0].Eta;
	} else {
	  leptonPT = -99;
	  leptonEta = -99;
	}

	if(nMuons == 1){
	  hMuonPT->Fill(selMuonT[0].PT, weight);
	  hMuonPhi->Fill(selMuonT[0].Phi, weight);
	  hMuonEta->Fill(selMuonT[0].Eta, weight);
	  muonPT = selMuonT[0].PT;
	  muonPhi = selMuonT[0].Phi;
	  muonEta = selMuonT[0].Eta;
	} else {
	  muonPT = -99;
	  muonPhi = -99;
	  muonEta = -99;
	}

	if(nEle == 1){
	  hElePT->Fill(selEle[0].PT, weight);
	  hElePhi->Fill(selEle[0].Phi, weight);
	  hEleEta->Fill(selEle[0].Eta, weight);
	  elePT = selEle[0].PT;
	  elePhi = selEle[0].Phi;
	  eleEta = selEle[0].Eta;
	} else {
	  elePT = -99;
	  elePhi = -99;
	  eleEta = -99;
	}

	/*	hChi2->Fill(minChi2, weight);
	hjetNumber->Fill(nJets, weight);
	hBjetNumber->Fill(nBJets, weight);
	hjetHT -> Fill(scalarSumPT, weight);
	hBjetHT -> Fill(scalarBjetSumPT, weight);
	hAvgDeltaRjj -> Fill(avgDeltaR, weight);
	hAvgDeltaRbb -> Fill(avgDeltaR_bb, weight);
	hAvgDeltaEtajj -> Fill(avgDeltaEta, weight);
	hAvgDeltaEtabb -> Fill(avgDeltaEta_bb, weight);
	hminDeltaRjj->Fill(minDeltaRjj, weight);
	hminDeltaRbb->Fill(minDeltaRbb, weight);
	hminDeltaRlj->Fill(minDeltaRlj, weight);
	hminDeltaRlb->Fill(minDeltaRlb, weight);
	hmaxDeltaEtabb -> Fill(maxDeltaEtabb, weight);
	hjetAverageMass ->Fill(avgJetMass, weight);
	hBjetAverageMass ->Fill(avgBjetMass, weight);
	hBjetAverageMassSqr ->Fill(avgBjetMassSqr, weight); */

	tree.Fill();

      }	      
      
  }

  /*  double meanMet = hmet->GetMean();
  double rmsMet = hmet->GetRMS(); 
  double entriesMet = hmet->GetEntries();
  cout <<"Mean Met: " << meanMet << ", RMS Met:" << rmsMet << ", Entries Met:" << entriesMet << endl; */

 
  jet->cd();

  /*  hmet->Write();
  hmetEta->Write();
  hmetPhi->Write();
  hjetPT1->Write();
  hjetPT2->Write();
  hjetPT3->Write();
  hjetPT4->Write();  
  hjetEta1->Write();
  hjetEta2->Write();
  hjetEta3->Write();
  hjetEta4->Write();
  hBjetPT1->Write();
  hBjetPT2->Write();
  hBjetPT3->Write();
  hBjetPT4->Write();  
  hBjetEta1->Write();
  hBjetEta2->Write();
  hBjetEta3->Write();
  hBjetEta4->Write();
  hAvgDeltaRjj->Write();
  hAvgDeltaRbb->Write();
  hAvgDeltaEtajj->Write();
  hAvgDeltaEtabb->Write();
  hminDeltaRjj->Write();
  hminDeltaRbb->Write();
  hminDeltaRlj->Write();
  hminDeltaRlb->Write();
  hmaxDeltaEtabb->Write();
  hjetAverageMass->Write();
  hBjetAverageMass->Write();
  hBjetAverageMassSqr->Write();
  hjetHT->Write();
  hBjetHT->Write();
  hjetNumber->Write();
  hBjetNumber->Write();
  hChi2->Write();  */

  hjetPT5->Write();
  hjetPT6->Write();
  hBjetPT5->Write();
  hBjetPT6->Write();
  hjetEta5->Write();
  hjetEta6->Write();
  hBjetEta5->Write();
  hBjetEta6->Write();

  lepton->cd();

  hLeptonPT->Write();
  hMuonPT->Write();
  hElePT->Write();
  hLeptonEta->Write();
  hMuonEta->Write();  
  hEleEta->Write();
  hLeptonPhi->Write();
  hMuonPhi->Write();
  hElePhi->Write();

  cout <<"starting to write the root file: " << endl;

  hfile.Write();
  hfile.Close();

  ev.close();
  of.close();
  return 0;
}



