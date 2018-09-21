//STANDARD ROOT INCLUDES
#include <TROOT.h>

#include <TStyle.h>

#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include <TMath.h>
#include <TLorentzVector.h>

#include <Math/Math.h>
#include <Math/LorentzVector.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

#include "TRandom3.h"

//STANDARD C++ INCLUDES

#include <cmath>

#include <vector>
#include <map>
#include <string>

#include <iostream>
#include <sstream>
#include <fstream>

#include <algorithm>

//#include "statusbar.h"

//#include "pid_tools.h"
//#include "shared_fill_info.h"

#define MASS_MU 0.1057 // GeV
#define MASS_E  0.000511 // GeV
#define MASS_P  0.938272029 // GeV
#define ECM 13000.0 // GeV

const float ns_to_s_ = 1e-9;
const float m_to_cm_ = 1e2;

using namespace std;
using namespace ROOT;

//bool showBar = true;

TFile* fOut;
TTree* trout;

//Particle flow meaning
//0 unknown
//1 charged (trk)
//2 electron (ch+ECAL)
//3 muon (trk+muon chamber)
//4 gamma (ecal and no trk)
//5 neutral hadron (HCal and no tracks)
//6 HF hadronics
//7 HF electromagnetic
///eos/totem/data/cmstotem/2015/90m/Merged/4510/9985

//#define preTS2 0

int preTS2;
int getXangle(int run,int lumi, const char* filename);

struct SortByPt{
   template <class T>
   bool operator()( const T& a, const T& b){
      return a.Pt() > b.Pt();
   }
};

struct ParticleWithIdx{
   unsigned int idx_;
   double pt_;
   ParticleWithIdx( unsigned int idx, double pt ) { idx_ = idx; pt_ = pt; } 
   double Pt() const { return pt_; } 
   double pt() const { return pt_; } 
};

struct FillInfo
{
  unsigned int fillNumber;
  unsigned int runMin, runMax;
  string alignmentTag;

  FillInfo(unsigned int _fn=0, unsigned int _rmin=0, unsigned int _rmax=0, const string &at=""):fillNumber(_fn), runMin(_rmin), runMax(_rmax), alignmentTag(at){}
};

//----------------------------------------------------------------------------------------------------

struct FillInfoCollection : public vector<FillInfo> 
{
  FillInfo FindByRun(unsigned int run) const
  {
    for (const auto it : *this)
    {
      if (it.runMin <= run && it.runMax >= run)
	return it;
    }
    for (const auto it : *this)
    {
      if (it.runMin <= 300000  && it.runMax >= 300000)
      {
	return it;
      }
    }
    printf("THIS INSTEAD SHOULD NOT HAPPEN. EXIT\n");
    exit(1);
    return FillInfo();
  }
};

FillInfoCollection fillInfoCollection;

void InitFillInfoCollection()
{
  fillInfoCollection.push_back(FillInfo(4947, 273725, 273730, "period1_physics_margin/fill_4947/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4953, 274094, 274094, "period1_physics_margin/fill_4953/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4961, 274198, 274200, "period1_physics_margin/fill_4961/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4964, 274240, 274241, "period1_physics_margin/fill_4964/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4964, 274244, 274244, "period1_physics/fill_4964/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4976, 274282, 274286, "period1_physics_margin/fill_4976/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4985, 274387, 274388, "period1_physics/fill_4985/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4988, 274420, 274422, "period1_physics/fill_4988/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4990, 274440, 274443, "period1_physics/fill_4990/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5005, 274954, 274959, "period1_physics/fill_5005/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5013, 274966, 274971, "period1_physics/fill_5013/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5017, 274998, 275001, "period1_physics/fill_5017/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5020, 275059, 275074, "period1_physics/fill_5020/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5021, 275124, 275125, "period1_physics/fill_5021/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5024, 275282, 275293, "period1_physics/fill_5024/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5026, 275309, 275311, "period1_physics/fill_5026/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5027, 275319, 275338, "period1_physics/fill_5027/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5028, 275344, 275345, "period1_physics/fill_5028/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5029, 275370, 275371, "period1_physics/fill_5029/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5030, 275375, 275376, "period1_physics/fill_5030/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5038, 275656, 275659, "period1_physics/fill_5038/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5043, 275757, 275783, "period1_physics/fill_5043/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5045, 275828, 275847, "period1_physics/fill_5045/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5048, 275886, 275890, "period1_physics/fill_5048/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5052, 275911, 275931, "period1_physics/fill_5052/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5261, 279760, 279767, "period1_physics/fill_5261/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5264, 279794, 279794, "period1_physics/fill_5264/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5265, 279823, 279823, "period1_physics/fill_5265/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5266, 279841, 279841, "period1_physics/fill_5266/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5267, 279844, 279865, "period1_physics/fill_5267/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5274, 279931, 279931, "period1_physics/fill_5274/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5275, 279966, 279966, "period1_physics/fill_5275/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5276, 279975, 279975, "period1_physics/fill_5276/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5277, 279993, 280024, "period1_physics/fill_5277/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5279, 280187, 280194, "period1_physics/fill_5279/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5287, 280327, 280364, "period1_physics/fill_5287/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5288, 280383, 280385, "period1_physics/fill_5288/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(8000, 300000, 300001, "period1_physics/fill_8000/2017_01_17"));

  fillInfoCollection.push_back(FillInfo(5393, 282730, 282735, "TODO"));
  fillInfoCollection.push_back(FillInfo(5401, 282920, 282924, "TODO"));
  fillInfoCollection.push_back(FillInfo(5405, 283040, 283043, "TODO"));
  fillInfoCollection.push_back(FillInfo(5406, 283049, 283067, "TODO"));
  fillInfoCollection.push_back(FillInfo(5418, 283305, 283308, "TODO"));
  fillInfoCollection.push_back(FillInfo(5421, 283353, 283359, "TODO"));
  fillInfoCollection.push_back(FillInfo(5423, 283407, 283416, "TODO"));
  fillInfoCollection.push_back(FillInfo(5424, 283453, 283453, "TODO"));
  fillInfoCollection.push_back(FillInfo(5427, 283478, 283481, "TODO"));
  fillInfoCollection.push_back(FillInfo(5433, 283548, 283560, "TODO"));
  fillInfoCollection.push_back(FillInfo(5437, 283672, 283685, "TODO"));
  fillInfoCollection.push_back(FillInfo(5439, 283820, 283835, "TODO"));
  fillInfoCollection.push_back(FillInfo(5441, 283863, 283865, "TODO"));
  fillInfoCollection.push_back(FillInfo(5442, 283876, 283878, "TODO"));
  fillInfoCollection.push_back(FillInfo(5443, 283884, 283885, "TODO"));
  fillInfoCollection.push_back(FillInfo(5446, 283933, 283934, "TODO"));
  fillInfoCollection.push_back(FillInfo(5448, 283946, 283964, "TODO"));
  fillInfoCollection.push_back(FillInfo(5450, 284006, 284014, "TODO"));
  fillInfoCollection.push_back(FillInfo(5451, 284025, 284044, "TODO"));
}

struct AlignmentResult
{
  double sh_x, sh_x_unc;		// mm
  double sh_y, sh_y_unc;		// mm

  AlignmentResult(double _sh_x=0., double _sh_x_unc=0., double _sh_y=0., double _sh_y_unc=0.) :
    sh_x(_sh_x), sh_x_unc(_sh_x_unc), sh_y(_sh_y), sh_y_unc(_sh_y_unc)
  {
  }

  void Write(FILE *f) const
  {
    fprintf(f, "sh_x=%.3f,sh_x_unc=%.3f,sh_y=%.3f,sh_y_unc=%.3f\n", sh_x, sh_x_unc, sh_y, sh_y_unc);
  }
};

struct AlignmentResults : public map<unsigned int, AlignmentResult>
{
  void Write(FILE *f) const
  {
    for (auto &p : *this)
    {
      fprintf(f, "id=%u,", p.first);
      p.second.Write(f);
    }
  }

  int Add(char *line)
  {
    bool idSet = false;
    unsigned int id = 0;
    AlignmentResult result;

    // loop over entries separated by ","
    char *p = strtok(line, ",");
    while (p != NULL)
    {
      // isolate key and value strings
      char *pe = strstr(p, "=");
      if (pe == NULL)
      {
	printf("ERROR in AlignmentResults::Add > entry missing = sign: %s.\n", p);
	return 2;
      }

      char *s_key = p;

      p = strtok(NULL, ",");
      *pe = 0;
      char *s_val = pe+1;

      // interprete keys
      if (strcmp(s_key, "id") == 0)
      {
	idSet = true;
	id = atoi(s_val);
	continue;
      }

      if (strcmp(s_key, "sh_x") == 0)
      {
	result.sh_x = atof(s_val);
	continue;
      }

      if (strcmp(s_key, "sh_x_unc") == 0)
      {
	result.sh_x_unc = atof(s_val);
	continue;
      }

      if (strcmp(s_key, "sh_y") == 0)
      {
	result.sh_y = atof(s_val);
	continue;
      }

      if (strcmp(s_key, "sh_y_unc") == 0)
      {
	result.sh_y_unc = atof(s_val);
	continue;
      }
      printf("ERROR in AlignmentResults::Add > unknown key: %s.\n", s_key);
      return 3;
    }

    if (!idSet)
    {
      printf("ERROR in AlignmentResults::Add > id not set on the following line:\n%s.\n", line);
      return 4;
    }
    insert({id, result});
    return 0;
  }

  double GETX(int rpnumber) const
  {
    auto ait = find(rpnumber);
    double toret=0;
    if (ait == end())
    {
      printf("ERROR in AlignmentResults::Apply > no alignment data for RP %u.\n", rpnumber);
      exit(1);
    }
    else
    {
      toret=ait->second.sh_x;
    }
    return toret;
  }

};

struct AlignmentResultsCollection : public map<string, AlignmentResults>
{
  int Write(const string &fn) const
  {
    FILE *f = fopen(fn.c_str(), "w");
    if (!f)
      return -1;
    Write(f);
    fclose(f);
    return 0;
  }

  void Write(FILE *f) const
  {
    for (auto &p : *this)
    {
      fprintf(f, "\n[%s]\n", p.first.c_str());
      p.second.Write(f);
    }
  }

  int Load(const string &fn)
  {
    FILE *f = fopen(fn.c_str(), "r");
    if (!f)
      return -1;
    return Load(f);
    fclose(f);
  }

  int Load(FILE *f)
  {
    string label = "unknown";
    AlignmentResults block;

    while (!feof(f))
    {
      char line[300];
      char *success = fgets(line, 300, f);
      if (success == NULL)
	break;
      if (line[0] == '\n')
	continue;
      if (line[0] == '[')
      {
	if (block.size() > 0)
	  insert({label, block});
	block.clear();
	char *end = strstr(line, "]");
	if (end == NULL)
	{
	  printf("ERROR in AlignmentResultsCollection::Load > missing closing ].\n");
	  return 1;
	}
	*end = 0;
	label = line+1;
	continue;
      }
      if (block.Add(line) != 0)
	return 2;
    }
    if (block.size() > 0)
      insert({label, block});
    return 0;
  }
};

TRandom3 *Rgenerator;

double Rndm(){
  return ( double(rand())/RAND_MAX );
}

class Particle{

  public:
    Particle();
    Particle(double);
    double   px, py, pz, E, m, pAbs, pt, p[4]; 
    void     p4(double, double, double, double);
    void     boost(Particle);
    void     boost2(double px_, double py_, double pz_, double E_);
    void     twoBodyDecay(Particle&, Particle&);
    void     print();
};

Particle::Particle(){
  px = py = pz = E = m = pAbs = pt = 0.0;
  p[0] = p[1] = p[2] = p[3] = 0.0;
}
Particle::Particle(double mass){
  m = mass;
  px = py = pz = E = pAbs = pt = 0.0;
  p[0] = p[1] = p[2] = p[3] = 0.0;
}

//
//*** Sets components of 4-momentum vector -------------------------------------
//
void Particle::p4(double momx, double momy, double momz, double energy){
  // components of 4-momenta
  px = p[0] = momx;
  py = p[1] = momy;
  pz = p[2] = momz;
  E  = p[3] = energy;
  // transverse momentum and the magnitude of the space momentum
  pt      = sqrt(momx*momx + momy*momy);
  pAbs    = sqrt(momx*momx + momy*momy + momz*momz);
}

//
//*** Prints 4-vector ----------------------------------------------------------
//
void Particle::print(){
  cout << "(" << p[0] <<",\t" << p[1] <<",\t"<< p[2] <<",\t"<< p[3] << ")" << endl;
}

//
//*** Evaluates 4-vector of decay product in the rest frame of parent. ---------
//
void Particle::twoBodyDecay(Particle& prod1, Particle& prod2) {

  double m1 = prod1.m;
  double m2 = prod2.m;

  // check if particle m can decay
  if( m < m1+m2 ){
    cout << "Particle::twoBodyDecay  parent mass is less than sum of products." 
      << endl;
    prod1.p4 ( 0., 0., 0., m1);
    prod2.p4 ( 0., 0., 0., m2);
    return;
  }

  // CM energies and momentum
  double e1 = (m*m + m1*m1 - m2*m2) / (2.0*m);
  double e2 = (m*m - m1*m1 + m2*m2) / (2.0*m);
  double P  = sqrt(e1*e1 - m1*m1);

  // Isotropic random angles
  double theta = acos( 2.0*Rndm() - 1.0 );
  double phi   = 2.0*M_PI *Rndm();

  double pX = P*sin(theta)*cos(phi);
  double pY = P*sin(theta)*sin(phi);
  double pZ = P*cos(theta);

  // set the 4-momenta
  prod1.p4(  pX,  pY,  pZ, e1 );
  prod2.p4( -pX, -pY, -pZ, e2 );
}

//
//*** Lorentz Boost  -----------------------------------------------------------
//
void Particle::boost(Particle parent){

  // beta and gamma values
  double betax = parent.px / parent.E;
  double betay = parent.py / parent.E;
  double betaz = parent.pz / parent.E;
  double beta2 = betax*betax + betay*betay + betaz*betaz;
  double gamma = 1.0/sqrt(1.0-beta2);
  double dot   = betax*px + betay*py + betaz*pz;
  double prod  = gamma*( gamma*dot/(1.0+gamma) + E );

  double pX = px + betax*prod;
  double pY = py + betay*prod;
  double pZ = pz + betaz*prod;
  double e  = gamma*(E + dot);

  p4( pX, pY, pZ, e );
}

void Particle::boost2(double px_, double py_, double pz_, double E_){

  // beta and gamma values
  double betax = px_ / E_;
  double betay = py_ / E_;
  double betaz = pz_ / E_;
  double beta2 = betax*betax + betay*betay + betaz*betaz;
  double gamma = 1.0/sqrt(1.0-beta2);
  double dot   = betax*px + betay*py + betaz*pz;
  double prod  = gamma*( gamma*dot/(1.0+gamma) + E );

  double pX = px + betax*prod;
  double pY = py + betay*prod;
  double pZ = pz + betaz*prod;
  double e  = gamma*(E + dot);

  p4( pX, pY, pZ, e );
}

void Dilepton(string year, vector<string> const& fileNames, string const& outputFileNameTree = "TTreeOutput.root",const Int_t nevt_max = 100)
{ 

  bool verbose = false;
  bool debug = false;

  bool isMC  = false;
  bool saveExtraTracks = false; 

  string treeName;
  //treeName = (!isMC) ? "ntp1" : "ntp1";
  treeName = (!isMC) ? "ggll_aod/ntp1" : "ntp1";
  map<string,TH1F*> histosTH1F;

  /*double energyMin = -10.;
  double energyMax = 190.;
  int nBinsEnergy = 1000;*/

  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;
  vector<TString>* vfiles = new vector<TString>; 

  for(size_t idx_file = 0; idx_file < fileNames.size(); ++idx_file)
  {
    vfiles->push_back( fileNames[idx_file] );
    std::cout<<"Loaded Input file: "<< fileNames[idx_file]<<std::endl;
  }

  // Declaration of tree and its branches variables
  TTree* tree = NULL;

  int nLocalProtCand;
  int therunnumber;

  unsigned short LocalProtCand_arm[100];
  unsigned short LocalProtCand_side[100];
  int HLT_Accept[5];

  unsigned int nPrimVertexCand;

  int const NPAIRMAX = 100;
  unsigned int nPair;
  int    Pair_lepton1[NPAIRMAX];
  int    Pair_lepton2[NPAIRMAX];
  double KalmanVertexCand_z[NPAIRMAX];
  int    Pair_extratracks0p5mm[NPAIRMAX];
  int    Pair_extratracks1mm[NPAIRMAX];
  int    Pair_extratracks2mm[NPAIRMAX];
  double Pair_mass[NPAIRMAX];
  double Pair_eta[NPAIRMAX];
  double Pair_phi[NPAIRMAX];
  double Pair_dphi[NPAIRMAX];
  double Pair_pt[NPAIRMAX];
  double ClosestExtraTrack_vtxdxyz[NPAIRMAX];
  double ClosestHighPurityExtraTrack_vtxdxyz[NPAIRMAX];

  int const NETMAX = 1000; 
  int          nExtraTracks;
  int          ExtraTrack_pair[NETMAX];
  int          ExtraTrack_purity[NETMAX];
  unsigned int ExtraTrack_nhits[NETMAX];
  int          ExtraTrack_charge[NETMAX];
  unsigned int ExtraTrack_ndof[NETMAX];
  double       ExtraTrack_px[NETMAX], 
               ExtraTrack_py[NETMAX], 
               ExtraTrack_pz[NETMAX];
  double       ExtraTrack_chi2[NETMAX];
  double       ExtraTrack_vtxdxyz[NETMAX];
  double       ExtraTrack_vtxT[NETMAX], ExtraTrack_vtxZ[NETMAX];
  double       ExtraTrack_x[NETMAX], 
               ExtraTrack_y[NETMAX], 
               ExtraTrack_z[NETMAX];
  /*std::vector<double> VExtraTrack_z;
  std::vector<double> VExtraTrack_px;
  std::vector<double> VExtraTrack_py;
  std::vector<double> VExtraTrack_pz;*/

  int const NMUMAX = 100;
  unsigned int nMuonCand;
  double MuonCand_pt[NMUMAX];
  double MuonCand_px[NMUMAX]; 
  double MuonCand_py[NMUMAX];
  double MuonCand_pz[NMUMAX];
  //double MuonCand_p[NMUMAX];
  double MuonCand_e[NMUMAX];
  int MuonCand_istight[NMUMAX];
  int MuonCand_isglobal[NMUMAX];
  double MuonCand_eta[NMUMAX];
  double MuonCand_phi[NMUMAX];
  int MuonCand_charge[NMUMAX];
  double MuonCand_vtxz[NMUMAX];

  double LocalProtCand_x[100];
  double LocalProtCand_z[100];

  int const NRPPixelMAX = 100; 
  unsigned int nPixelTrack=0;
  double PixelTrackX[NRPPixelMAX];
  double PixelTrackY[NRPPixelMAX];
  int    PixelTrackArm[NRPPixelMAX];

  unsigned int Run=0;
  unsigned int LumiSection=0;
  unsigned int EventNum=0;

  /*double TimingRecHitT[800];
  double TimingRecHitToT[800];
  int TimingRecHitArm[800];
  int nRecHitsTiming = 0;*/
  int const NRPTTMAX = 100; 
  unsigned int nRPTimingTrack = 0;
  double       RPTimingTrack_Time[NRPTTMAX];
  int          RPTimingTrack_arm[NRPTTMAX];

  int chargeVal=0;
  //int compactRpID=0;

  double pixelcoordX[2][150];
  double pixelcoordY[2][150];

  double LeadingMuonPt=0.;
  double SecondMuonPt=0.;
  double LeadingMuonEta=0.;
  double SecondMuonEta=0.;
  double LeadingMuonPhi=0.;
  double SecondMuonPhi=0.;
  double LeadingMuonVtxZ=0.;
  double SecondMuonVtxZ=0.;
  int LeadingMuonTightID=0;
  int SecondMuonTightID=0;

  double VtxZ    = 0.;
  double VtxZPPS = 0.;
  int nExtraTracks0p5mm = -1; 
  int nExtraTracks1mm   = -1; 
  int nExtraTracks2mm   = -1; 

  double Muons_pxtot=0.;
  double Muons_pytot=0.;
  double Muons_pztot=0.;
  double Muons_pttot=0.;
  double pairMass=0.;
  double pairPhi=0.;
  double pairEta=0.;
  double pairRapidity=0.;
  double pairpT=0.;
  double acoplan=0.;
  double closestExtraTrack=0.;
  //int nExtraTrack0pt5mm = 0;
  //int nExtraTrack1mm = 0;
  double Muons_Etot=0.;
  double align_trk0=0.; double align_trk1=0.;
  double corrected_X_0=0.;
  double corrected_X_1=0.;
  double proton1_pz=0.;
  double proton2_pz=0.;
  double protonmass=0.;
  double missingmass=0.;
  double xi1=0.;
  double xi2=0.;
  int angle = -1;
  int trigger = 0;
  double pixelXarm0=0.;
  double pixelXarm1=0.;
  double pixelYarm0=0.;
  double pixelYarm1=0.;

  /*double AverageLeadingSec45 = 0.;
  double AverageToTSec45 = 0.;
  double AverageLeadingSec56 = 0.;
  double AverageToTSec56 = 0.;*/
  double AverageTimeSec45 = 0.;
  double AverageTimeSec56 = 0.;

  int selectProtons = 0;

  int const NPMAX = 20;
  int nPixel_L = 0;
  int nPixel_R = 0;
  double xi_L[NPMAX];
  double xi_R[NPMAX];
  double pixelX_L[NPMAX];
  double pixelX_R[NPMAX];
  double pixelY_L[NPMAX];
  double pixelY_R[NPMAX];
  double proton_pz_L[NPMAX];
  double proton_pz_R[NPMAX];
  for(int k=0; k < NPMAX; ++k) { xi_L[k] = 0.; xi_R[k] = 0.; pixelX_L[k] = 0.; pixelX_R[k] = 0.; pixelY_L[k] = 0.; pixelY_R[k] = 0.; proton_pz_L[k] = 0; proton_pz_R[k] = 0;}

  std::vector<ROOT::Math::XYZTVector> muons;
  std::vector<ROOT::Math::XYZTVector>* muons_ptr = &muons;

  /*TH1D *hxi_mumu_L = new TH1D("xi_mumu_L", "xi_mumu_L", 500, -1,1);
  TH1D *hxi_mumu_R = new TH1D("xi_mumu_R", "xi_mumu_R", 500, -1,1);
  TH1D *hxi_L = new TH1D("xi_L", "xi_L", 500, -1,1);
  TH1D *hxi_R = new TH1D("xi_R", "xi_R", 500, -1,1);
  TH2D* hxi_correlation = new TH2D("xiCorrelation", "xiCorrelation", 500, -1.,1., 500,-1.,1.);*/

  TString ttreefilename = outputFileNameTree;
  fOut = new TFile(ttreefilename, "RECREATE");
  fOut->cd();

  trout = new TTree("Events", "Events");
  // Event info
  trout->Branch("Run"         , &Run         , "Run/I");
  trout->Branch("LumiSection" , &LumiSection , "LumiSection/I");
  trout->Branch("EventNum"    , &EventNum    , "EventNum/I");
  trout->Branch("Trigger"     , &trigger     , "Trigger/I");
  // Muon info
  trout->Branch("muons","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&muons_ptr);
  trout->Branch("nMuonCand"          , &nMuonCand          , "nMuonCand/I");
  trout->Branch("LeadingMuonPt"      , &LeadingMuonPt      , "LeadingMuonPt/D");
  trout->Branch("LeadingMuonEta"     , &LeadingMuonEta     , "LeadingMuonEta/D");
  trout->Branch("LeadingMuonPhi"     , &LeadingMuonPhi     , "LeadingMuonPhi/D");
  trout->Branch("LeadingMuonVtxZ"    , &LeadingMuonVtxZ    , "LeadingMuonVtxZ/D");
  trout->Branch("LeadingMuonTightID" , &LeadingMuonTightID , "LeadingMuonTightID/I");
  trout->Branch("SecondMuonPt"       , &SecondMuonPt       , "SecondMuonPt/D");
  trout->Branch("SecondMuonEta"      , &SecondMuonEta      , "SecondMuonEta/D");
  trout->Branch("SecondMuonPhi"      , &SecondMuonPhi      , "SecondMuonPhi/D");
  trout->Branch("SecondMuonVtxZ"     , &SecondMuonVtxZ     , "SecondMuonVtxZ/D");
  trout->Branch("SecondMuonTightID"  , &SecondMuonTightID  , "SecondMuonTightID/I");
  trout->Branch("DimuonMass"         , &pairMass           , "DimuonMass/D");
  trout->Branch("DimuonEta"          , &pairEta            , "DimuonEta/D");
  trout->Branch("DimuonRapidity"     , &pairRapidity       , "DimuonRapidity/D");
  trout->Branch("DimuonPhi"          , &pairPhi            , "DimuonPhi/D");
  trout->Branch("DimuonPt"           , &pairpT             , "DimuonPt/D");
  trout->Branch("Acoplanarity"       , &acoplan            , "Acoplanarity/D");
  trout->Branch("ChargeDimuon"       , &chargeVal          , "ChargeDimuon/I");
  trout->Branch("FittedVtxZ"         , &VtxZ               , "FittedVtxZ/D");
  trout->Branch("nExtraTracks0p5mm"  , &nExtraTracks0p5mm  , "nExtraTracks0p5mm/I");
  trout->Branch("nExtraTracks1mm"    , &nExtraTracks1mm    , "nExtraTracks1mm/I");
  trout->Branch("nExtraTracks2mm"    , &nExtraTracks2mm    , "nExtraTracks2mm/I");
  // Extra track info 
  if(saveExtraTracks){
     trout->Branch( "nExtraTracks", &nExtraTracks, "nExtraTracks/I" );
     trout->Branch( "ExtraTrack_pair", ExtraTrack_pair, "ExtraTrack_pair[nExtraTracks]/I" );
     trout->Branch( "ExtraTrack_purity", ExtraTrack_purity, "ExtraTrack_purity[nExtraTracks]/I" );
     trout->Branch( "ExtraTrack_nhits", ExtraTrack_nhits, "ExtraTrack_nhits[nExtraTracks]/I" );
     trout->Branch( "ExtraTrack_charge", ExtraTrack_charge, "ExtraTrack_charge[nExtraTracks]/I" );
     trout->Branch( "ExtraTrack_ndof", ExtraTrack_ndof, "ExtraTrack_ndof[nExtraTracks]/I" );
     trout->Branch( "ExtraTrack_px", ExtraTrack_px, "ExtraTrack_px[nExtraTracks]/D" );
     trout->Branch( "ExtraTrack_py", ExtraTrack_py, "ExtraTrack_py[nExtraTracks]/D" );
     trout->Branch( "ExtraTrack_pz", ExtraTrack_pz, "ExtraTrack_pz[nExtraTracks]/D" );
     trout->Branch( "ExtraTrack_chi2", ExtraTrack_chi2, "ExtraTrack_chi2[nExtraTracks]/D" );
     trout->Branch( "ExtraTrack_vtxdxyz", ExtraTrack_vtxdxyz, "ExtraTrack_vtxdxyz[nExtraTracks]/D" );
     trout->Branch( "ExtraTrack_vtxT", ExtraTrack_vtxT, "ExtraTrack_vtxT[nExtraTracks]/D" );
     trout->Branch( "ExtraTrack_vtxZ", ExtraTrack_vtxZ, "ExtraTrack_vtxZ[nExtraTracks]/D" );
     trout->Branch( "ExtraTrack_x", ExtraTrack_x, "ExtraTrack_x[nExtraTracks]/D" );
     trout->Branch( "ExtraTrack_y", ExtraTrack_y, "ExtraTrack_y[nExtraTracks]/D" );
     trout->Branch( "ExtraTrack_z", ExtraTrack_z, "ExtraTrack_z[nExtraTracks]/D" );
  } 
  trout->Branch("DistanceClosestExtraTrack",&closestExtraTrack,"DistanceClosestExtraTrack/D");
  // Proton/RP info 
  trout->Branch("nLocalProtCand",&nLocalProtCand,"nLocalProtCand/I");
  trout->Branch("SelectProtons",&selectProtons,"SelectProtons/I");
  trout->Branch("SelectedProtonPzArm0",&proton1_pz,"SelectedProtonPzArm0/D");
  trout->Branch("SelectedProtonPzArm1",&proton2_pz,"SelectedProtonPzArm1/D");
  trout->Branch("SelectedProtonXiArm0",&xi1,"SelectedProtonXiArm0/D");
  trout->Branch("SelectedProtonXiArm1",&xi2,"SelectedProtonXiArm1/D");
  trout->Branch("ProtonMass",&protonmass,"ProtonMass/D");
  // Missing mass 
  trout->Branch("MissingMass",&missingmass,"MissingMass/D");
  if(year=="2017"){
     // Event info
     trout->Branch("CrossingAngle",&angle,"CrossingAngle/I");
     // Proton/RP info 
     trout->Branch("nPixelArm0",&nPixel_L,"nPixelArm0/I");
     trout->Branch("nPixelArm1",&nPixel_R,"nPixelArm1/I");
     trout->Branch("PixelXArm0",&pixelX_L,"PixelXArm0[nPixelArm0]/D");
     trout->Branch("PixelXArm1",&pixelX_R,"PixelXArm1[nPixelArm1]/D");
     trout->Branch("PixelYArm0",&pixelY_L,"PixelYArm0[nPixelArm0]/D");
     trout->Branch("PixelYArm1",&pixelY_R,"PixelYArm1[nPixelArm1]/D");
     trout->Branch("ProtonXiArm0",&xi_L,"ProtonXiArm0[nPixelArm0]/D");
     trout->Branch("ProtonXiArm1",&xi_R,"ProtonXiArm1[nPixelArm1]/D");
     trout->Branch("ProtonPzArm0",&proton_pz_L,"ProtonPzArm0[nPixelArm0]/D");
     trout->Branch("ProtonPzArm1",&proton_pz_R,"ProtonPzArm1[nPixelArm1]/D");
     /*trout->Branch("nRecHitsTiming",&nRecHitsTiming, "nRecHitsTiming/I");
     trout->Branch("AverageLeadingSec45",&AverageLeadingSec45, "AverageLeadingSec45/D");
     trout->Branch("AverageLeadingSec56",&AverageLeadingSec56, "AverageLeadingSec56/D");
     trout->Branch("AverageToTSec45",&AverageToTSec45, "AverageToTSec45/D");
     trout->Branch("AverageToTSec56",&AverageToTSec56, "AverageToTSec56/D");*/
     trout->Branch("AverageTimeSec45",&AverageTimeSec45, "AverageTimeSec45");
     trout->Branch("AverageTimeSec56",&AverageTimeSec56, "AverageTimeSec56");
     trout->Branch("VtxZFromPPS",&VtxZPPS,"VtxZFromPPS/D");
  }

  int i_tot = 0 , nevt_tot = 0; int indicatorefile=-1;

  for(vector<TString>::iterator itfiles = vfiles->begin() ; itfiles != vfiles->end() && i_tot < nevt_max_corr ; ++itfiles){
     std::cout<<"File loop:"<<indicatorefile<<std::endl;
     TFile* file = TFile::Open(*itfiles,"READ");
     indicatorefile++;

     // Access TTree from current file
     tree = (TTree*) file->Get( treeName.c_str() );
     if( tree == nullptr ) exit (EXIT_FAILURE);
    
     int nev = int(tree->GetEntriesFast());
     nevt_tot += nev;
     cout <<"The current file ("<< indicatorefile <<") has " << nev << " entries : "<< endl << *itfiles << endl;
     std::vector<unsigned int> rp_list;

     //  to TTree ----------------------------------------------------------------------
     tree->SetBranchAddress("Run",&Run);  

     if(year=="2016"){
	tree->SetBranchAddress("nLocalProtCand",&nLocalProtCand);
	tree->SetBranchAddress("LocalProtCand_x",&LocalProtCand_x);
	tree->SetBranchAddress("LocalProtCand_z",&LocalProtCand_z);
     }else if (year=="2017"){
	/*tree->SetBranchAddress("nPixelTracks", &nPixelTracks);
	tree->SetBranchAddress("PixTrackX", &PixTrackX);
	tree->SetBranchAddress("PixTrackY", &PixTrackY);
	tree->SetBranchAddress("PixTrackArm", &PixTrackArm);
	tree->SetBranchAddress("nRecHitsTiming", &nRecHitsTiming);
	tree->SetBranchAddress("TimingRecHitArm",&TimingRecHitArm);
	tree->SetBranchAddress("TimingRecHitT",&TimingRecHitT);
	tree->SetBranchAddress("TimingRecHitToT",&TimingRecHitToT);*/
	tree->SetBranchAddress("nRPPixelTrack"     , &nPixelTrack);
	tree->SetBranchAddress("RPPixelTrack_x"    , &PixelTrackX);
	tree->SetBranchAddress("RPPixelTrack_y"    , &PixelTrackY);
	tree->SetBranchAddress("RPPixelTrack_arm"  , &PixelTrackArm);
	tree->SetBranchAddress("nRPTimingTrack"    , &nRPTimingTrack);
	tree->SetBranchAddress("RPTimingTrack_arm" , &RPTimingTrack_arm);
	tree->SetBranchAddress("RPTimingTrack_Time", &RPTimingTrack_Time);
     }else{
	exit (EXIT_FAILURE);
     }
     tree->SetBranchAddress("EventNum", &EventNum);
     tree->SetBranchAddress("LumiSection", &LumiSection);
     tree->SetBranchAddress("HLT_Accept", &HLT_Accept);
     tree->SetBranchAddress("nMuonCand",&nMuonCand);
     tree->SetBranchAddress("nPrimVertexCand", &nPrimVertexCand);
     tree->SetBranchAddress("MuonCand_vtxz",&MuonCand_vtxz);
     tree->SetBranchAddress("MuonCand_pt",&MuonCand_pt);
     //tree->SetBranchAddress("MuonCand_px",&MuonCand_px);    
     //tree->SetBranchAddress("MuonCand_py",&MuonCand_py);    
     //tree->SetBranchAddress("MuonCand_pz",&MuonCand_pz);    
     tree->SetBranchAddress("MuonCand_eta",&MuonCand_eta);
     tree->SetBranchAddress("MuonCand_phi",&MuonCand_phi);
     //tree->SetBranchAddress("MuonCand_p",&MuonCand_p);
     tree->SetBranchAddress("MuonCand_e",&MuonCand_e);
     tree->SetBranchAddress("MuonCand_charge",&MuonCand_charge);    
     tree->SetBranchAddress("MuonCand_isglobal",&MuonCand_isglobal);    
     tree->SetBranchAddress("MuonCand_istight",&MuonCand_istight);    
     tree->SetBranchAddress("nPair",&nPair);
     tree->SetBranchAddress("Pair_lepton1", &Pair_lepton1);
     tree->SetBranchAddress("Pair_lepton2", &Pair_lepton2);
     tree->SetBranchAddress("Pair_mass", &Pair_mass);
     tree->SetBranchAddress("Pair_pt", &Pair_pt);
     tree->SetBranchAddress("Pair_eta", &Pair_eta);
     tree->SetBranchAddress("Pair_phi", &Pair_phi);
     tree->SetBranchAddress("Pair_dphi", &Pair_dphi);
     tree->SetBranchAddress("KalmanVertexCand_z",&KalmanVertexCand_z);
     tree->SetBranchAddress("Pair_extratracks0p5mm",&Pair_extratracks0p5mm);
     tree->SetBranchAddress("Pair_extratracks1mm",&Pair_extratracks1mm);
     tree->SetBranchAddress("Pair_extratracks2mm",&Pair_extratracks2mm);
     if(saveExtraTracks){
	tree->SetBranchAddress( "nExtraTracks", &nExtraTracks);
	tree->SetBranchAddress( "ExtraTrack_pair", ExtraTrack_pair);
	tree->SetBranchAddress( "ExtraTrack_purity", ExtraTrack_purity);
	tree->SetBranchAddress( "ExtraTrack_nhits", ExtraTrack_nhits);
	tree->SetBranchAddress( "ExtraTrack_charge", ExtraTrack_charge);
	tree->SetBranchAddress( "ExtraTrack_ndof", ExtraTrack_ndof);
	tree->SetBranchAddress( "ExtraTrack_px", ExtraTrack_px);
	tree->SetBranchAddress( "ExtraTrack_py", ExtraTrack_py);
	tree->SetBranchAddress( "ExtraTrack_pz", ExtraTrack_pz);
	tree->SetBranchAddress( "ExtraTrack_chi2", ExtraTrack_chi2);
	tree->SetBranchAddress( "ExtraTrack_vtxdxyz", ExtraTrack_vtxdxyz);
	tree->SetBranchAddress( "ExtraTrack_vtxT", ExtraTrack_vtxT);
	tree->SetBranchAddress( "ExtraTrack_vtxZ", ExtraTrack_vtxZ);
	tree->SetBranchAddress( "ExtraTrack_x", ExtraTrack_x);
	tree->SetBranchAddress( "ExtraTrack_y", ExtraTrack_y);
	tree->SetBranchAddress( "ExtraTrack_z", ExtraTrack_z);
     }
     tree->SetBranchAddress("ClosestExtraTrack_vtxdxyz", &ClosestExtraTrack_vtxdxyz);

     //  to TTree ----------------------------------------------------------------------
     //  load alignment collection
     AlignmentResultsCollection alignmentCollection;
     alignmentCollection.Load("./analysis_ctpps_alignment/shared_alignment/collect_alignments.out");
     AlignmentResults *alignments = NULL;
     InitFillInfoCollection();

     int NpixelLL;
     int NpixelRR;
     int pixelLL;
     int pixelRR;

     double px_1= 0.;
     double py_1 = 0.;
     double pz_1= 0.;
     double pE_1 = 0.;
     double px_2 = 0.;
     double py_2 = 0.;
     double pz_2  = 0.;
     double pE_2= 0.;
     double px_3 = 0.;
     double py_3 = 0.;
     double pz_3  = 0.;
     double pE_3= 0.;
     double px_4 = 0.;
     double py_4 = 0.;
     double pz_4  = 0.;
     double pE_4= 0.;

     double alignX=0.;
     double alignX_2=0.;
     double alignX_3=0.;
     double alignX_102=0.;
     double alignX_103=0.;

     double fakedispersion_mm=80.;

     for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot) {

        if(debug) cout << "Analyzing event " << ( i_tot + 1 ) << endl;

	/*if(showBar){
	  if (i_evt==0) {
	  std::cout << "" << std::endl;
	  std::cout<< "Status Bar" << std::endl;
	  }
	  loadBar(i_evt, nev);
	//tree->GetEntry(i_evt);
	}*/
	if ( i_evt && i_evt % 10000 == 0 )
	   cout << i_evt << " events analyzed." << endl;
	tree->GetEntry(i_evt);

	for(int i=0;i<150;i++){
	   pixelcoordX[0][i]=0;
	   pixelcoordX[1][i]=0;
	   pixelcoordY[0][i]=0;
	   pixelcoordY[1][i]=0;
	}
	// ---
	//compactRpID=0;

	NpixelLL=0;
	NpixelRR=0;
	pixelLL=0;
	pixelRR=0;

	chargeVal=0;

	Muons_pxtot=0.;
	Muons_pytot=0.;
	Muons_pztot=0.;
	Muons_pttot=0.;
	Muons_Etot=0.;

	LeadingMuonPt=0;
	SecondMuonPt=0;

	LeadingMuonEta=0;
	SecondMuonEta=0;

	LeadingMuonPhi=0;
	SecondMuonPhi=0;

	LeadingMuonTightID=0;
	SecondMuonTightID=0;

	LeadingMuonVtxZ=0;
	SecondMuonVtxZ=0;
	VtxZ=0;
	VtxZPPS=0;
        nExtraTracks0p5mm = -1;
        nExtraTracks1mm = -1;
        nExtraTracks2mm = -1;

	pairMass=0.;
	pairPhi=0.;
	pairEta=0.;
	pairRapidity=0.;
	pairpT=0.;
	acoplan=0.;
	closestExtraTrack=0.;

	align_trk0=0.;
	align_trk1=0.;
	corrected_X_0=0.;
	corrected_X_1=0.;

	proton1_pz=0.;
	proton2_pz=0.;
	protonmass=0.;
	missingmass=0.;
	xi1=0;
	xi2=0;

	trigger = 0;
	angle = -1;
	pixelXarm0=0;
	pixelXarm1=0;
	pixelYarm0=0;
	pixelYarm1=0;

	/*AverageLeadingSec45 = 0.;
	AverageToTSec45 = 0.;
	AverageLeadingSec56 = 0.;
	AverageToTSec56 = 0.;*/
        AverageTimeSec45 = 0.;
        AverageTimeSec56 = 0.;

        selectProtons = 0;
 
        nPixel_L = 0;
        nPixel_R = 0;
	for(int k=0; k < NPMAX; ++k) { xi_L[k] = 0.; xi_R[k] = 0.; pixelX_L[k] = 0.; pixelX_R[k] = 0.; pixelY_L[k] = 0.; pixelY_R[k] = 0.; proton_pz_L[k] = 0; proton_pz_R[k] = 0; }

        muons.clear();

        therunnumber = Run;

	if(HLT_Accept[0] || HLT_Accept[1] || HLT_Accept[2] ) trigger = 1;
        bool triggerAccept = false;
        if( trigger == 1 ) triggerAccept = true;

        int nLocalProtCand_L = 0;
        int nLocalProtCand_R = 0;
	if(year=="2017"){
	   nLocalProtCand = 0;
	   nLocalProtCand_L = 0;
	   nLocalProtCand_R = 0;
	   for(int s = 0; s < nPixelTrack; ++s)
	   {
	      if(PixelTrackArm[s]==0){
		 pixelcoordX[0][NpixelLL]=PixelTrackX[s]; 
		 pixelcoordY[0][NpixelLL]=PixelTrackY[s]; 
		 NpixelLL++; 
		 pixelLL=s; 
		 nLocalProtCand++;
		 nLocalProtCand_L++;
	      }
	      if(PixelTrackArm[s]==1){
		 pixelcoordX[1][NpixelRR]=PixelTrackX[s]; 
		 pixelcoordY[1][NpixelRR]=PixelTrackY[s]; 
		 NpixelRR++;
		 pixelRR=s; 
		 nLocalProtCand++;
		 nLocalProtCand_R++;
	      }
	      if(debug){ 
		 std::cout << "pixel X [mm]: " << PixelTrackX[s] << ", ARM: " << PixelTrackArm[s] << std::endl;
		 std::cout << "pixel Y [mm]: " << PixelTrackY[s] << ", ARM: " << PixelTrackArm[s] << std::endl;
	      }
	   }
	} else if(year=="2016"){

	   const auto &fillInfo = fillInfoCollection.FindByRun(therunnumber);
	   if(fillInfo.fillNumber==8000)
	      continue;

	   alignX=0.;
	   alignX_2=0.;
	   alignX_3=0.;
	   alignX_102=0.;
	   alignX_103=0.;

	   const auto alignment_it = alignmentCollection.find(fillInfo.alignmentTag);
	   if (alignment_it == alignmentCollection.end())
	   {
	      printf("ERROR: no alignment for tag '%s'.\n", fillInfo.alignmentTag.c_str());
	      continue;
	   }
	   else
	   {
	      alignments = &alignment_it->second;	      
	      alignX_2=alignments->GETX(2);
	      alignX_3=alignments->GETX(3);
	      alignX_102=alignments->GETX(102);
	      alignX_103=alignments->GETX(103);
	   }
	}

	if(debug) std::cout << "nLocalProtCand: " << nLocalProtCand << std::endl;

	int anyproton=0;

	fakedispersion_mm=80.;

	if(nLocalProtCand > 0){

	   //if(debug) std::cout << "Analyzing event with at least one proton..." << std::endl;
	   anyproton++;

	   //if(nLocalProtCand_L == 1 && nLocalProtCand_R == 1)

	   // Calculate total PT of the system
	   /*LeadingMuonPt=0.;
	   SecondMuonPt=0.;
	   LeadingMuonEta=0.;
	   SecondMuonEta=0.;
	   LeadingMuonPhi=0.;
	   SecondMuonPhi=0.;
	   LeadingMuonVtxZ=0.;
	   SecondMuonVtxZ=0.;
	   VtxZ=0.;
	   VtxZpps=0.;
	   nExtraTracks1mm = 0;
	   LeadingMuonTightID=0.;
	   SecondMuonTightID=0.;
	   Muons_pxtot=0.;
	   Muons_pytot=0.;
	   Muons_pztot=0.;
	   Muons_pttot=0.;
	   Muons_Etot=0.;*/

           //FIXME 
	   align_trk0=0.; align_trk1=0.;
	   if(year=="2016"){
	      if(LocalProtCand_z[0]>0)
		 align_trk0=alignX_2;
	      else
		 align_trk0=alignX_102;
	      if(LocalProtCand_z[1]>0)
		 align_trk1=alignX_2;
	      else
		 align_trk1=alignX_102;
	      corrected_X_0=align_trk0+1000.*LocalProtCand_x[0];
	      corrected_X_1=align_trk1+1000.*LocalProtCand_x[1];
	   }

	   if(year=="2017"){

	      angle = -1;
	      double dispersionL=-9999;
	      double dispersionR=-9999;

	      double Dispersion[2][200];
	      Dispersion[0][120]=9.145*10;
	      Dispersion[0][130]=8.591*10;
	      Dispersion[0][140]=8.226*10;
	      Dispersion[1][120]=7.291*10;
	      Dispersion[1][130]=6.621*10;
	      Dispersion[1][140]=6.191*10;

	      double relcorr=0;
	      preTS2 ? relcorr=(42.05):relcorr=(42.2);

	      if(preTS2==1){
		 angle=getXangle(therunnumber,  LumiSection , "./inp/new/xangle_tillTS2_STABLEBEAMS_CLEANUP.csv");
	      }

	      if(preTS2==0) angle=getXangle(therunnumber,  LumiSection , "./inp/xangle_afterTS2_cleanup.csv");
	      if(angle==120 || angle ==130 || angle==140){dispersionL=Dispersion[0][angle];dispersionR=Dispersion[1][angle];}
	      if(angle !=120 && angle != 130 && angle != 140){
		 dispersionL = Dispersion[0][140] - 4.6*((angle-140.)/10.) ;
		 dispersionR = Dispersion[1][140] - 4.6*((angle-140.)/10.) ;
	      }

	      if(angle > 0){
		 /*if(NpixelLL==1){
		   for(int k=0;k<NpixelLL;k++){
		   xi1=(pixelcoordX[0][k]-relcorr)/dispersionL;
		   proton1_pz=(1-(xi1))*(ECM/2.);
		   pixelXarm0=pixelcoordX[0][k];
		   pixelYarm0=pixelcoordY[0][k];
		   }
		   }else{
		   pixelXarm0=0.;
		   pixelYarm0=0.;
		   }*/
		 if(NpixelLL > 0){
		    for(int k=0;k<NpixelLL && k<NPMAX;k++){
		       xi_L[k]        = (pixelcoordX[0][k]-relcorr)/dispersionL;
		       proton_pz_L[k] = (1-(xi_L[k]))*(ECM/2.);
		       pixelX_L[k]    = pixelcoordX[0][k];
		       pixelY_L[k]    = pixelcoordY[0][k];
		       ++nPixel_L; 
		    }
		 }
		 /*if(NpixelRR==1){
		   for(int k=0;k<NpixelRR;k++){
		   xi2=(pixelcoordX[1][k]-relcorr)/dispersionR;
		   proton2_pz=-(1-(xi2))*(ECM/2.);
		   pixelXarm1=pixelcoordX[1][k];
		   pixelYarm1=pixelcoordY[1][k];
		   }
		   }else{
		   pixelXarm1=0.;
		   pixelYarm1=0.;
		   }*/
		 if(NpixelRR > 0){
		    for(int k=0;k<NpixelRR && k<NPMAX;k++){
		       xi_R[k]        = (pixelcoordX[1][k]-relcorr)/dispersionR;
		       proton_pz_R[k] = (1-(xi_R[k]))*(ECM/2.);
		       pixelX_R[k]    = pixelcoordX[1][k];
		       pixelY_R[k]    = pixelcoordY[1][k];
		       ++nPixel_R;
		    }
		 }
		 //protonmass=ECM*sqrt((  (xi1)*(xi2) ));
	      }
	   }
	   /* 2016
	      {
	      proton1_pz=(1-(corrected_X_0/fakedispersion_mm))*(ECM/2.);
	      proton2_pz=-(1-(corrected_X_1/fakedispersion_mm))*(ECM/2.);
	      xi1 = corrected_X_0/fakedispersion_mm;
	      xi2 = corrected_X_1/fakedispersion_mm;
	      protonmass=ECM*sqrt(( (corrected_X_0/fakedispersion_mm)*(corrected_X_1/fakedispersion_mm) ));
	      }*/
	}

	TLorentzVector muon0;
	TLorentzVector muon1;
	TLorentzVector muonpair;
	TLorentzVector proton1;
	TLorentzVector proton2;

	bool muonsAccept = false;
	if(nMuonCand >= 2) muonsAccept = true;
	//if(protonsAccept)
	if( triggerAccept && muonsAccept ){

           // Index of leading and second leading muons in pT
           vector<ParticleWithIdx> muonsWithIdx;
           muonsWithIdx.clear();
           for(unsigned int idx_muons = 0; idx_muons < nMuonCand; ++idx_muons) muonsWithIdx.push_back( ParticleWithIdx( idx_muons, MuonCand_pt[idx_muons] ) );
           stable_sort( muonsWithIdx.begin(), muonsWithIdx.end(), SortByPt() ); 
           unsigned int idx_mu0 = muonsWithIdx.at(0).idx_;  
           unsigned int idx_mu1 = muonsWithIdx.at(1).idx_;  
           if(debug) cout << "Found muon with index " << idx_mu0 << endl
                          << "Found muon with index " << idx_mu1 << endl;
               
	   LeadingMuonPt=MuonCand_pt[idx_mu0];
	   SecondMuonPt=MuonCand_pt[idx_mu1];

	   LeadingMuonEta=MuonCand_eta[idx_mu0];
	   SecondMuonEta=MuonCand_eta[idx_mu1];

	   LeadingMuonPhi=MuonCand_phi[idx_mu0];
	   SecondMuonPhi=MuonCand_phi[idx_mu1];

	   LeadingMuonTightID=MuonCand_istight[idx_mu0];
	   SecondMuonTightID=MuonCand_istight[idx_mu1];

	   LeadingMuonVtxZ=MuonCand_vtxz[idx_mu0];
	   SecondMuonVtxZ=MuonCand_vtxz[idx_mu1];

           Math::RhoEtaPhiVector vec3D_muon0( MuonCand_pt[idx_mu0], MuonCand_eta[idx_mu0], MuonCand_phi[idx_mu0] );
           Math::XYZTVector vec4D_muon0( vec3D_muon0.X(), vec3D_muon0.Y(), vec3D_muon0.Z(), MuonCand_e[idx_mu0] );

           Math::RhoEtaPhiVector vec3D_muon1( MuonCand_pt[idx_mu1], MuonCand_eta[idx_mu1], MuonCand_phi[idx_mu1] );
           Math::XYZTVector vec4D_muon1( vec3D_muon1.X(), vec3D_muon1.Y(), vec3D_muon1.Z(), MuonCand_e[idx_mu1] );
           muons.resize(2);
           muons.at(0) = vec4D_muon0;
           muons.at(1) = vec4D_muon1;

	   /*VtxZ = KalmanVertexCand_z[0];
	   nExtraTracks1mm = Pair_extratracks1mm[0];
	   closestExtraTrack = ClosestExtraTrack_vtxdxyz[0];*/
   
           // Look for pair from leading muons
           bool foundPair = false;
           for(unsigned int idx_pairs = 0; idx_pairs < nPair && !foundPair; ++idx_pairs){
              if( ( Pair_lepton1[idx_pairs] == idx_mu0 && Pair_lepton2[idx_pairs] == idx_mu1 ) || 
                  ( Pair_lepton1[idx_pairs] == idx_mu1 && Pair_lepton2[idx_pairs] == idx_mu0 ) ){

                 foundPair = true;

                 if(debug) cout << "Found pair with index " << idx_pairs << "." 
                                << " First lepton: " << Pair_lepton1[idx_pairs] << "."
                                << " Second lepton: " << Pair_lepton2[idx_pairs] << "." << endl;
 
	         VtxZ = KalmanVertexCand_z[idx_pairs];
                 nExtraTracks0p5mm = Pair_extratracks0p5mm[idx_pairs];
	         nExtraTracks1mm   = Pair_extratracks1mm[idx_pairs];
	         nExtraTracks2mm   = Pair_extratracks2mm[idx_pairs];
                 
	         closestExtraTrack = ClosestExtraTrack_vtxdxyz[idx_pairs];
              }
           }

	   //muon0.SetPtEtaPhiE(MuonCand_pt[0], MuonCand_eta[0], MuonCand_phi[0], MuonCand_p[0]);
	   //muon1.SetPtEtaPhiE(MuonCand_pt[1], MuonCand_eta[1], MuonCand_phi[1], MuonCand_p[1]);
	   //muon0.SetPtEtaPhiE( MuonCand_pt[idx_mu0], MuonCand_eta[idx_mu0], MuonCand_phi[idx_mu0], TMath::Sqrt(MuonCand_p[idx_mu0]*MuonCand_p[idx_mu0] + MASS_MU*MASS_MU) );
	   //muon1.SetPtEtaPhiE( MuonCand_pt[idx_mu1], MuonCand_eta[idx_mu1], MuonCand_phi[idx_mu1], TMath::Sqrt(MuonCand_p[idx_mu1]*MuonCand_p[idx_mu1] + MASS_MU*MASS_MU) );
	   muon0.SetPtEtaPhiE( MuonCand_pt[idx_mu0], MuonCand_eta[idx_mu0], MuonCand_phi[idx_mu0], MuonCand_e[idx_mu0] );
	   muon1.SetPtEtaPhiE( MuonCand_pt[idx_mu1], MuonCand_eta[idx_mu1], MuonCand_phi[idx_mu1], MuonCand_e[idx_mu1] );
	   muonpair = muon0+muon1;

	   //Muons_Etot=sqrt(MuonCand_p[0]*MuonCand_p[0]+MASS_MU*MASS_MU) + sqrt(MuonCand_p[1]*MuonCand_p[1]+MASS_MU*MASS_MU);
	   Muons_Etot =  muon0.E() + muon1.E();
	   Muons_pxtot = muon0.Px() + muon1.Px();
	   Muons_pytot = muon0.Py() + muon1.Py();
	   Muons_pztot = muon0.Pz() + muon1.Pz();
	   chargeVal=(int)MuonCand_charge[0]+(int)MuonCand_charge[1];

	   pairMass=muonpair.M();
	   pairEta=muonpair.Eta();
	   pairRapidity=muonpair.Rapidity();
	   pairPhi=muonpair.Phi();
	   pairpT=muonpair.Pt();

	   double dphi = ( muon0.Phi() - muon1.Phi() );
           if( dphi >= Math::Pi() ) dphi -= 2*Math::Pi();
           if( dphi < -Math::Pi() ) dphi += 2*Math::Pi();
	   //acoplan = 1 - std::abs(Pair_dphi[0])/pi;
	   acoplan = 1 - std::abs(dphi) / Math::Pi();

	   /*VExtraTrack_z.clear();
	     VExtraTrack_px.clear();
	     VExtraTrack_py.clear();
	     VExtraTrack_pz.clear();
	     for(int k=0; k<21; k++){
	     VExtraTrack_z.push_back(ExtraTrack_z[k]);
	     VExtraTrack_px.push_back(ExtraTrack_px[k]);
	     VExtraTrack_py.push_back(ExtraTrack_py[k]);
	     VExtraTrack_pz.push_back(ExtraTrack_pz[k]);
	     }*/

	   bool protonsAccept = false;
	   // Events with only one proton per side
	   if(year=="2016"){
	      if((LocalProtCand_z[0]>0 && LocalProtCand_z[1]<0) || (LocalProtCand_z[0]<0 && LocalProtCand_z[1]>0)) protonsAccept = true;
	      if(protonsAccept){

		 selectProtons = 1;

		 proton1_pz=(1-(corrected_X_0/fakedispersion_mm))*(ECM/2.);
		 proton2_pz=-(1-(corrected_X_1/fakedispersion_mm))*(ECM/2.);

		 xi1 = corrected_X_0/fakedispersion_mm;
		 xi2 = corrected_X_1/fakedispersion_mm;

		 //protonmass=ECM*sqrt(( (corrected_X_0/fakedispersion_mm)*(corrected_X_1/fakedispersion_mm) ));
		 protonmass = ECM*sqrt( (xi1)*(xi2) );
	      }
	   } else if(year=="2017"){
	      if(NpixelRR==1 && NpixelLL==1) protonsAccept = true;
	      if(protonsAccept){

		 selectProtons = 1;

		 xi1 = xi_L[0];
		 proton1_pz = proton_pz_L[0];
		 pixelXarm0 = pixelX_L[0];
		 pixelYarm0 = pixelY_L[0];

		 xi2 = xi_R[0];
		 proton2_pz = proton_pz_R[0];
		 pixelXarm1 = pixelX_R[0];
		 pixelYarm1 = pixelY_R[0];

		 protonmass = ECM*sqrt( (xi1)*(xi2) );
	      }
	   }

	   if(protonsAccept){
	      double proton1_energy = std::abs( proton1_pz ); 
	      double proton2_energy = std::abs( proton2_pz ); 
	      proton1.SetPxPyPzE(0., 0.,  std::abs( proton1_pz ), proton1_energy);
	      proton2.SetPxPyPzE(0., 0., -std::abs( proton2_pz ), proton2_energy);
	      //missingmass = sqrt((ECM-(Muons_Etot+fabs(proton1_pz)+fabs(proton2_pz)))*(ECM-(Muons_Etot+fabs(proton1_pz)+fabs(proton2_pz)))-(Muons_pxtot)*(Muons_pxtot)-(Muons_pytot)*(Muons_pytot)- (Muons_pztot+proton1_pz+proton2_pz)*(Muons_pztot+proton1_pz+proton2_pz));
	      TLorentzVector TotalMom; TotalMom.SetPxPyPzE(0., 0., 0., 0.);
	      TotalMom += muonpair;
	      TotalMom += proton1;
	      TotalMom += proton2;
	      missingmass = ( -TotalMom ).M(); 

	      if(year=="2017"){

		 /*double TotalLeadingSec45 = 0.;
		 double TotalToTSec45 = 0.;
		 double TotalLeadingSec56 = 0.;
		 double TotalToTSec56 = 0.;*/

		 /*AverageLeadingSec45 = 0.;
		 AverageToTSec45 = 0.;
		 AverageLeadingSec56 = 0.;
		 AverageToTSec56 = 0.;*/

		 // Need to check with Ksenia which kind of container did she used. If DIGI, check with Nicola about 25./1024 factor.
		 /*for(int k=0; k<nRecHitsTiming; k++){
		    if(TimingRecHitArm[k]==0){
		       TotalLeadingSec45 += TimingRecHitT[k];
		       TotalToTSec45 += TimingRecHitToT[k];
		    }
		    if(TimingRecHitArm[k]==1){
		       TotalLeadingSec56 += TimingRecHitT[k];
		       TotalToTSec56 += TimingRecHitToT[k];
		    }
		 }

		 AverageLeadingSec45 = TotalLeadingSec45/nRecHitsTiming;
		 AverageToTSec45 = TotalToTSec45/nRecHitsTiming;
		 AverageLeadingSec56 = TotalLeadingSec56/nRecHitsTiming;
		 AverageToTSec56 = TotalToTSec56/nRecHitsTiming;

		 VtxZpps=(AverageLeadingSec45-AverageLeadingSec56) * ns_to_s_ * TMath::C() * 0.5 * m_to_cm_;*/

		 AverageTimeSec45 = 0.;
		 AverageTimeSec56 = 0.;

		 double TotalTimeSec45 = 0.;
		 double TotalTimeSec56 = 0.;

		 int NSec45 = 0;
		 int NSec56 = 0;
		 for(int k=0; k<nRPTimingTrack; k++){
		    if(RPTimingTrack_arm[k]==0){
		       TotalTimeSec45 += RPTimingTrack_Time[k];
		       NSec45++;
		    }
		    if(RPTimingTrack_arm[k]==1){
		       TotalTimeSec56 += RPTimingTrack_Time[k];
		       NSec56++;
		    }
		 }

		 double deltaT = 0;

		 if(NSec45>0) AverageTimeSec45 = TotalTimeSec45/NSec45;
		 else AverageTimeSec45 = 0.;
		 if(NSec56>0) AverageTimeSec56 = TotalTimeSec56/NSec56;
		 else AverageTimeSec56 = 0.;
		 if(NSec45 > 0 && NSec56 > 0) deltaT = ( AverageTimeSec45 - AverageTimeSec56 );

		 VtxZPPS = deltaT * ns_to_s_ * TMath::C() * 0.5 * m_to_cm_;

	      }

	   }

           // Fill tree
	   trout->Fill();

	} // Dilepton selection
     } // End of loop over events

     file->Close();

  } // End of loop over files

  cout << "Total of " << i_tot << " events analyzed." << endl;

  fOut->cd();
  trout->Write();
  fOut->Close();

}

//----------------------------------------------------------------------------------------------------


int main(int argc, char** argv)
{

  std::cout << "Main is called with " << argc << " arguments:" <<" Expected: int maxevent, string extrastr, era (e.g. B,C,G), string year"<< std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout << argv[i] << std::endl;
  }
  int maxevent=atoi(argv[1]);
  string extrastr=argv[2];
  string era_str=argv[3];
  string year_str=argv[4];


  std::cout<<"CALLING MAIN.."<<std::endl;
  std::vector<string> allnamesWithPrefixFinal;

  if(year_str=="2016"){
    preTS2 = 0;
    std::cout << "Analyzing 2016 data..." << std::endl;
    if(era_str=="B")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/jjhollar/output_mumu_merge_data_2016B_23Sep2016v3_FullFinal_pT40.root");
    if(era_str=="C")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/jjhollar/output_mumu_merge_data_2016C_23Sep2016v1_FullFinal_pT40.root");
    if(era_str=="G")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/jjhollar/output_mumu_merge_data_2016G_23Sep2016v1_FullFinal_pT40.root");
  } else if(year_str=="2017"){
    std::cout << "Analyzing 2017 data..." << std::endl;
    if(era_str=="B"||era_str=="C"||era_str=="D") preTS2 = 1;
    if(era_str=="E"||era_str=="F") preTS2 = 0;
    /*if(era_str=="B")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/kshcheli/NewTriggers/DoubleMuon/DoubleMuRunB2017gt.root");
    if(era_str=="C")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/kshcheli/2017DileptonsAndWW/DoubleMuon/DoubleMu2017_C.root");
    if(era_str=="D")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/kshcheli/2017DileptonsAndWW/DoubleMuon/DoubleMu2017_D.root");
    if(era_str=="E")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/kshcheli/2017DileptonsAndWW/DoubleMuon/DoubleMu2017_E.root");
    if(era_str=="F")allnamesWithPrefixFinal.push_back("/eos/cms/store/user/kshcheli/2017DileptonsAndWW/DoubleMuon/DoubleMu2017_F_fixed.root");*/
    if(era_str=="B")allnamesWithPrefixFinal.push_back("/eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonB.root");
    if(era_str=="C")allnamesWithPrefixFinal.push_back("/eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonC.root");
    if(era_str=="D")allnamesWithPrefixFinal.push_back("/eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonD.root");
    if(era_str=="E")allnamesWithPrefixFinal.push_back("/eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonE.root");
    if(era_str=="F")allnamesWithPrefixFinal.push_back("/eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonF.root");
  }else{
    std::cout << "Please, put the right year." << std::endl;
    return 0;
  }

  string filenameTree = "TTree" + era_str + "_" + extrastr + ".root";
  Dilepton(year_str, allnamesWithPrefixFinal, filenameTree, maxevent);

  return 0;
}

int getXangle(int run, int lumi, const char* filename)
{
  TString drun;
  TString dfill;
  TString dlumi;
  TString temp;
  int Xangle=-1;
  TString runs;runs.Form("%d", run);
  TString lumis;lumis.Form("%d", lumi);
  ifstream F;
  F.open((const char*)filename);
  //F.clear();
  //F.seekg(0, ios::beg);
  int counter=0;
  if(F){
    while (!F.eof())
    {
      F>>drun;
      F>>dfill;
      F>>dlumi;
      F>>temp;
      F>>Xangle;
      //if(runfill==runlumi && )
      if( runs == drun &&  lumis==dlumi )
      {
	//  cout << "Read: "<< run<<" " << lumi<<"    -----> "<<Xangle << endl;   
	break;
	//return Xangle;
      }
      //    prevrunfill=runfill;
      //    prevXangle=Xangle;
    }
  }//endif
  else cout << "[!] gerXangle():Error reading from file" << endl;
  if(F.eof())
  {
    cout << "[getXangle() warning:] No Xangle data for this run found!" <<endl;
    F.close();
    return -1;
  }
  else  {F.close(); //cout << "~~~~~~~ xangle: "<< Xangle << endl;
    return Xangle;
  }
}
