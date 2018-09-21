//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 20 16:57:19 2018 by ROOT version 6.10/09
// from TTree Events/Events
// found on file: TTreeB_test.root
//////////////////////////////////////////////////////////

#ifndef Analysis_h
#define Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "Math/GenVector/PxPyPzE4D.h"

class Analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxmuons = 2;

   // Declaration of leaf types
   Int_t           Run;
   Int_t           LumiSection;
   Int_t           EventNum;
   Int_t           Trigger;
   Int_t           muons_;
   Double_t        muons_fCoordinates_fX[kMaxmuons];   //[muons_]
   Double_t        muons_fCoordinates_fY[kMaxmuons];   //[muons_]
   Double_t        muons_fCoordinates_fZ[kMaxmuons];   //[muons_]
   Double_t        muons_fCoordinates_fT[kMaxmuons];   //[muons_]
   Int_t           nMuonCand;
   Double_t        LeadingMuonPt;
   Double_t        LeadingMuonEta;
   Double_t        LeadingMuonPhi;
   Double_t        LeadingMuonVtxZ;
   Int_t           LeadingMuonTightID;
   Double_t        SecondMuonPt;
   Double_t        SecondMuonEta;
   Double_t        SecondMuonPhi;
   Double_t        SecondMuonVtxZ;
   Int_t           SecondMuonTightID;
   Double_t        DimuonMass;
   Double_t        DimuonEta;
   Double_t        DimuonRapidity;
   Double_t        DimuonPhi;
   Double_t        DimuonPt;
   Double_t        Acoplanarity;
   Int_t           ChargeDimuon;
   Double_t        FittedVtxZ;
   Int_t           nExtraTracks0p5mm;
   Int_t           nExtraTracks1mm;
   Int_t           nExtraTracks2mm;
   Double_t        DistanceClosestExtraTrack;
   Int_t           nLocalProtCand;
   Int_t           SelectProtons;
   Double_t        SelectedProtonPzArm0;
   Double_t        SelectedProtonPzArm1;
   Double_t        SelectedProtonXiArm0;
   Double_t        SelectedProtonXiArm1;
   Double_t        ProtonMass;
   Double_t        MissingMass;
   Int_t           CrossingAngle;
   Int_t           nPixelArm0;
   Int_t           nPixelArm1;
   Double_t        PixelXArm0[10];   //[nPixelArm0]
   Double_t        PixelXArm1[10];   //[nPixelArm1]
   Double_t        PixelYArm0[10];   //[nPixelArm0]
   Double_t        PixelYArm1[10];   //[nPixelArm1]
   Double_t        ProtonXiArm0[10];   //[nPixelArm0]
   Double_t        ProtonXiArm1[10];   //[nPixelArm1]
   Double_t        ProtonPzArm0[10];   //[nPixelArm0]
   Double_t        ProtonPzArm1[10];   //[nPixelArm1]
   Float_t         AverageTimeSec45;
   Float_t         AverageTimeSec56;
   Double_t        VtxZFromPPS;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_LumiSection;   //!
   TBranch        *b_EventNum;   //!
   TBranch        *b_Trigger;   //!
   TBranch        *b_muons_;   //!
   TBranch        *b_muons_fCoordinates_fX;   //!
   TBranch        *b_muons_fCoordinates_fY;   //!
   TBranch        *b_muons_fCoordinates_fZ;   //!
   TBranch        *b_muons_fCoordinates_fT;   //!
   TBranch        *b_nMuonCand;   //!
   TBranch        *b_LeadingMuonPt;   //!
   TBranch        *b_LeadingMuonEta;   //!
   TBranch        *b_LeadingMuonPhi;   //!
   TBranch        *b_LeadingMuonVtxZ;   //!
   TBranch        *b_LeadingMuonTightID;   //!
   TBranch        *b_SecondMuonPt;   //!
   TBranch        *b_SecondMuonEta;   //!
   TBranch        *b_SecondMuonPhi;   //!
   TBranch        *b_SecondMuonVtxZ;   //!
   TBranch        *b_SecondMuonTightID;   //!
   TBranch        *b_DimuonMass;   //!
   TBranch        *b_DimuonEta;   //!
   TBranch        *b_DimuonRapidity;   //!
   TBranch        *b_DimuonPhi;   //!
   TBranch        *b_DimuonPt;   //!
   TBranch        *b_Acoplanarity;   //!
   TBranch        *b_ChargeDimuon;   //!
   TBranch        *b_FittedVtxZ;   //!
   TBranch        *b_nExtraTracks0p5mm;   //!
   TBranch        *b_nExtraTracks1mm;   //!
   TBranch        *b_nExtraTracks2mm;   //!
   TBranch        *b_DistanceClosestExtraTrack;   //!
   TBranch        *b_nLocalProtCand;   //!
   TBranch        *b_SelectProtons;   //!
   TBranch        *b_SelectedProtonPzArm0;   //!
   TBranch        *b_SelectedProtonPzArm1;   //!
   TBranch        *b_SelectedProtonXiArm0;   //!
   TBranch        *b_SelectedProtonXiArm1;   //!
   TBranch        *b_ProtonMass;   //!
   TBranch        *b_MissingMass;   //!
   TBranch        *b_CrossingAngle;   //!
   TBranch        *b_nPixelArm0;   //!
   TBranch        *b_nPixelArm1;   //!
   TBranch        *b_PixelXArm0;   //!
   TBranch        *b_PixelXArm1;   //!
   TBranch        *b_PixelYArm0;   //!
   TBranch        *b_PixelYArm1;   //!
   TBranch        *b_ProtonXiArm0;   //!
   TBranch        *b_ProtonXiArm1;   //!
   TBranch        *b_ProtonPzArm0;   //!
   TBranch        *b_ProtonPzArm1;   //!
   TBranch        *b_AverageTimeSec45;   //!
   TBranch        *b_AverageTimeSec56;   //!
   TBranch        *b_VtxZFromPPS;   //!

   Analysis(TTree *tree=0);
   virtual ~Analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analysis_cxx
Analysis::Analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TTreeB_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("TTreeB_test.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

Analysis::~Analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("LumiSection", &LumiSection, &b_LumiSection);
   fChain->SetBranchAddress("EventNum", &EventNum, &b_EventNum);
   fChain->SetBranchAddress("Trigger", &Trigger, &b_Trigger);
   fChain->SetBranchAddress("muons", &muons_, &b_muons_);
   fChain->SetBranchAddress("muons.fCoordinates.fX", muons_fCoordinates_fX, &b_muons_fCoordinates_fX);
   fChain->SetBranchAddress("muons.fCoordinates.fY", muons_fCoordinates_fY, &b_muons_fCoordinates_fY);
   fChain->SetBranchAddress("muons.fCoordinates.fZ", muons_fCoordinates_fZ, &b_muons_fCoordinates_fZ);
   fChain->SetBranchAddress("muons.fCoordinates.fT", muons_fCoordinates_fT, &b_muons_fCoordinates_fT);
   fChain->SetBranchAddress("nMuonCand", &nMuonCand, &b_nMuonCand);
   fChain->SetBranchAddress("LeadingMuonPt", &LeadingMuonPt, &b_LeadingMuonPt);
   fChain->SetBranchAddress("LeadingMuonEta", &LeadingMuonEta, &b_LeadingMuonEta);
   fChain->SetBranchAddress("LeadingMuonPhi", &LeadingMuonPhi, &b_LeadingMuonPhi);
   fChain->SetBranchAddress("LeadingMuonVtxZ", &LeadingMuonVtxZ, &b_LeadingMuonVtxZ);
   fChain->SetBranchAddress("LeadingMuonTightID", &LeadingMuonTightID, &b_LeadingMuonTightID);
   fChain->SetBranchAddress("SecondMuonPt", &SecondMuonPt, &b_SecondMuonPt);
   fChain->SetBranchAddress("SecondMuonEta", &SecondMuonEta, &b_SecondMuonEta);
   fChain->SetBranchAddress("SecondMuonPhi", &SecondMuonPhi, &b_SecondMuonPhi);
   fChain->SetBranchAddress("SecondMuonVtxZ", &SecondMuonVtxZ, &b_SecondMuonVtxZ);
   fChain->SetBranchAddress("SecondMuonTightID", &SecondMuonTightID, &b_SecondMuonTightID);
   fChain->SetBranchAddress("DimuonMass", &DimuonMass, &b_DimuonMass);
   fChain->SetBranchAddress("DimuonEta", &DimuonEta, &b_DimuonEta);
   fChain->SetBranchAddress("DimuonRapidity", &DimuonRapidity, &b_DimuonRapidity);
   fChain->SetBranchAddress("DimuonPhi", &DimuonPhi, &b_DimuonPhi);
   fChain->SetBranchAddress("DimuonPt", &DimuonPt, &b_DimuonPt);
   fChain->SetBranchAddress("Acoplanarity", &Acoplanarity, &b_Acoplanarity);
   fChain->SetBranchAddress("ChargeDimuon", &ChargeDimuon, &b_ChargeDimuon);
   fChain->SetBranchAddress("FittedVtxZ", &FittedVtxZ, &b_FittedVtxZ);
   fChain->SetBranchAddress("nExtraTracks0p5mm", &nExtraTracks0p5mm, &b_nExtraTracks0p5mm);
   fChain->SetBranchAddress("nExtraTracks1mm", &nExtraTracks1mm, &b_nExtraTracks1mm);
   fChain->SetBranchAddress("nExtraTracks2mm", &nExtraTracks2mm, &b_nExtraTracks2mm);
   fChain->SetBranchAddress("DistanceClosestExtraTrack", &DistanceClosestExtraTrack, &b_DistanceClosestExtraTrack);
   fChain->SetBranchAddress("nLocalProtCand", &nLocalProtCand, &b_nLocalProtCand);
   fChain->SetBranchAddress("SelectProtons", &SelectProtons, &b_SelectProtons);
   fChain->SetBranchAddress("SelectedProtonPzArm0", &SelectedProtonPzArm0, &b_SelectedProtonPzArm0);
   fChain->SetBranchAddress("SelectedProtonPzArm1", &SelectedProtonPzArm1, &b_SelectedProtonPzArm1);
   fChain->SetBranchAddress("SelectedProtonXiArm0", &SelectedProtonXiArm0, &b_SelectedProtonXiArm0);
   fChain->SetBranchAddress("SelectedProtonXiArm1", &SelectedProtonXiArm1, &b_SelectedProtonXiArm1);
   fChain->SetBranchAddress("ProtonMass", &ProtonMass, &b_ProtonMass);
   fChain->SetBranchAddress("MissingMass", &MissingMass, &b_MissingMass);
   fChain->SetBranchAddress("CrossingAngle", &CrossingAngle, &b_CrossingAngle);
   fChain->SetBranchAddress("nPixelArm0", &nPixelArm0, &b_nPixelArm0);
   fChain->SetBranchAddress("nPixelArm1", &nPixelArm1, &b_nPixelArm1);
   fChain->SetBranchAddress("PixelXArm0", PixelXArm0, &b_PixelXArm0);
   fChain->SetBranchAddress("PixelXArm1", PixelXArm1, &b_PixelXArm1);
   fChain->SetBranchAddress("PixelYArm0", PixelYArm0, &b_PixelYArm0);
   fChain->SetBranchAddress("PixelYArm1", PixelYArm1, &b_PixelYArm1);
   fChain->SetBranchAddress("ProtonXiArm0", ProtonXiArm0, &b_ProtonXiArm0);
   fChain->SetBranchAddress("ProtonXiArm1", ProtonXiArm1, &b_ProtonXiArm1);
   fChain->SetBranchAddress("ProtonPzArm0", ProtonPzArm0, &b_ProtonPzArm0);
   fChain->SetBranchAddress("ProtonPzArm1", ProtonPzArm1, &b_ProtonPzArm1);
   fChain->SetBranchAddress("AverageTimeSec45", &AverageTimeSec45, &b_AverageTimeSec45);
   fChain->SetBranchAddress("AverageTimeSec56", &AverageTimeSec56, &b_AverageTimeSec56);
   fChain->SetBranchAddress("VtxZFromPPS", &VtxZFromPPS, &b_VtxZFromPPS);
   Notify();
}

Bool_t Analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analysis_cxx
