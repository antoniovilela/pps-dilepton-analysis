#define Analysis_cxx

#include "Analysis.h"

#include <TStyle.h>

#include <TCanvas.h>

#include <TH1D.h>
#include <TH2D.h>

#include <TMath.h>
#include <TLorentzVector.h>

#include <iostream>

#define MASS_MU 0.1057 // GeV
#define MASS_E  0.000511 // GeV
#define MASS_P  0.938272029 // GeV
#define ECM 13000.0 // GeV

using namespace std;

void Analysis::Loop() {
   //   In a ROOT session, you can do:
   //      root> .L Analysis.C
   //      root> Analysis t
   //      root> t.GetEntry(12); // Fill t data members with entry number 12
   //      root> t.Show();       // Show values of entry 12
   //      root> t.Show(16);     // Read and show values of entry 16
   //      root> t.Loop();       // Loop on all entries
   //

   //     This is the loop skeleton where:
   //    jentry is the global entry number in the chain
   //    ientry is the entry number in the current Tree
   //  Note that the argument to GetEntry must be:
   //    jentry for TChain::GetEntry
   //    ientry for TTree::GetEntry and TBranch::GetEntry
   //
   //       To read only selected branches, Insert statements like:
   // METHOD1:
   //    fChain->SetBranchStatus("*",0);  // disable all branches
   //    fChain->SetBranchStatus("branchname",1);  // activate branchname
   // METHOD2: replace line
   //    fChain->GetEntry(jentry);       //read all branches
   //by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   map<string,TH1D*> histosTH1D;
   histosTH1D["Acoplanarity"] = new TH1D("Acoplanarity","Acoplanarity",100,0.,1.);
   histosTH1D["nExtraTracks0p5mm"] = new TH1D("nExtraTracks0p5mm","nExtraTracks0p5mm",150,0,150);
   histosTH1D["nExtraTracks1mm"] = new TH1D("nExtraTracks1mm","nExtraTracks1mm",150,0,150);
   histosTH1D["nExtraTracks2mm"] = new TH1D("nExtraTracks2mm","nExtraTracks2mm",150,0,150);
  
   map<string,TH2D*> histosTH2D;
   histosTH2D["Log10AcoplanarityVsLog10DistanceClosestTrack"] = new TH2D("Log10AcoplanarityVsLog10DistanceClosestTrack", "Log10AcoplanarityVsLog10DistanceClosestTrack", 200, -5., 0., 200, -5., 1.);
   histosTH2D["xiCorrelation_SingleHit_Arm0"] = new TH2D("xiCorrelation_SingleHit_Arm0", "xiCorrelation_SingleHit_Arm0", 500, -1., 1., 500, -1., 1.);
   histosTH2D["xiCorrelation_SingleHit_Arm1"] = new TH2D("xiCorrelation_SingleHit_Arm1", "xiCorrelation_SingleHit_Arm1", 500, -1., 1., 500, -1., 1.);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      /* Plot leading and second muon pT, eta, dimuon mass, pT, rapidity, vertex, 
	 number of extra tracks, acoplanbarity, distance of closest track to vertex, etc.*/ 
      if( LeadingMuonPt > 50. && SecondMuonPt > 50. && 
	    TMath::Abs(LeadingMuonEta) < 2.4 && TMath::Abs(SecondMuonEta) < 2.4 &&
	    LeadingMuonTightID && SecondMuonTightID ){

	 if( DimuonMass > 110. ){

	    if( TMath::Abs(FittedVtxZ) < 15. ){

               histosTH1D.at("Acoplanarity")->Fill( Acoplanarity );
               histosTH2D.at("Log10AcoplanarityVsLog10DistanceClosestTrack")->Fill( TMath::Log10(DistanceClosestExtraTrack), TMath::Log10(Acoplanarity) );

               histosTH1D.at("nExtraTracks0p5mm")->Fill(nExtraTracks0p5mm);
               histosTH1D.at("nExtraTracks1mm")->Fill(nExtraTracks1mm);
               histosTH1D.at("nExtraTracks2mm")->Fill(nExtraTracks2mm);

	       if( nExtraTracks1mm == 0){
		  // Other cuts (acoplanarity, distance of closest track to vertex, etc.)
		  TLorentzVector muon1, muon2;
		  muon1.SetPtEtaPhiM(LeadingMuonPt,LeadingMuonEta,LeadingMuonPhi,MASS_MU);
		  muon2.SetPtEtaPhiM(SecondMuonPt,SecondMuonEta,SecondMuonPhi,MASS_MU);

                  cout << "Leading muon pT, eta, phi: " << muon1.Pt() << ", " << muon1.Eta() << ", " << muon1.Phi() << endl; 
                  cout << "Second muon pT, eta, phi:  " << muon2.Pt() << ", " << muon2.Eta() << ", " << muon2.Phi() << endl; 

		  Double_t xill_L = 0., xill_R = 0.;
		  xill_L = ( muon1.Pt()*TMath::Exp( +1.*muon1.Eta() ) ) + ( muon2.Pt()*TMath::Exp( +1.*muon2.Eta() ) );
		  xill_R = ( muon1.Pt()*TMath::Exp( -1.*muon1.Eta() ) ) + ( muon2.Pt()*TMath::Exp( -1.*muon2.Eta() ) );
		  xill_L /= ECM;
		  xill_R /= ECM;

		  // Events with one hit on each side
		  if( nPixelArm0 == 1){
		     Double_t xi_proton_L = ProtonXiArm0[0];

		     cout << "Proton (L)" << " - xi(p): " << xi_proton_L << " - xi(mumu): " << xill_L << endl;

		     histosTH2D.at("xiCorrelation_SingleHit_Arm0")->Fill( xi_proton_L, xill_L );
		  }
		  if( nPixelArm1 == 1){
		     Double_t xi_proton_R = ProtonXiArm1[0];

		     cout << "Proton (R)" << " - xi(p): " << xi_proton_R << " - xi(mumu): " << xill_R << endl;

		     histosTH2D.at("xiCorrelation_SingleHit_Arm1")->Fill( xi_proton_R, xill_R );
		  }

	       } 
	    } 
	 } 
      }
   } // End of event loop

   TFile* output = new TFile("output.root","RECREATE");
   output->cd();
   for( auto& item: histosTH1D ) item.second->Write();
   for( auto& item: histosTH2D ) item.second->Write();
   output->Close();  
}
