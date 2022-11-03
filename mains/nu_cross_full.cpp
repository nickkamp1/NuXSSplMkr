#include "lhapdf_cross_section.h"

//#define SINGLE_ENERGY_TEST
#define SAVE_PATH "/data/user/lfischer/software/NuXSSplMkr/data/HNL/"

int main(int argc, char* argv[]){

//#ifndef SINGLE_ENERGY_TEST
//  if(argc != 2){
//    cerr << "Argument number not valid! Given: " <<  argc << endl;
//    cerr << "Usage: pdfname " << endl;
//    return 1;
//  }
//#else
//  if(argc != 5){
//    cerr << "Argument number not valid! Given: " <<  argc << endl;
//    cerr << "Usage: pdfname enu x y" << endl;
//    return 1;
//  }
//#endif

  PhysConst * pc = new PhysConst();

  // Get arguments
  std::string pdfname = (string) argv[1];
  double mass_double = atof(argv[2])/1000;  // mass must be given in MeV
  std::ostringstream ss;
  ss << std::setw(4) << std::setfill('0') << (string) argv[2];
  std::string mass_string = ss.str();
  LHAXS xs_obj(pdfname);

  //enum IntType {CC,NC};
  enum NeutrinoType {neutrino,antineutrino};
  enum PDFVar {central,minus,plus};

  std::map<Current,std::string> IntTypeLabel {{CC,"cc"},{NC,"nc"}};
  std::map<NeutrinoType,double> CP_factor {{neutrino,1.},{antineutrino,-1}};
  //std::map<NeutrinoType,std::string> NeutrinoTypeLabel {{neutrino,"numu"},{antineutrino,"numubar"}};
  std::map<NeutrinoType,std::string> NeutrinoTypeLabel {{neutrino,"nutau"},{antineutrino,"nutaubar"}};
  std::map<PDFVar,int> PDFVarIndex {{central,0},{minus,-1},{plus,1}};
  std::map<PDFVar,std::string> PDFVarLabel {{central,"central"},{minus,"minus"},{plus,"plus"}};

  // muon mass
  //xs_obj.Set_M_Lepton(0.105*xs_obj.pc->GeV);
  // HNL mass
  std::cout << "HNL mass: " << mass_double << std::endl;
  xs_obj.Set_M_Lepton(mass_double*xs_obj.pc->GeV);

  double cm2 = SQ(pc->cm);
  double m2 = SQ(pc->meter);

  for (Current IT : {CC,NC}) {
    std::cout << "Interaction Type: " << IT << std::endl;
    xs_obj.Set_InteractionType(IT);
    for (NeutrinoType neutype : {neutrino,antineutrino}){
      std::cout << "Neutrino Type: " << neutype << std::endl;
      xs_obj.Set_CP_factor(CP_factor[neutype]);
      for (PDFVar pdfvar : {central}){
        // neutrino = 1., antineutrino = -1
        //std::string filename_dsdy = static_cast<std::string>(SAVE_PATH) + "dsdy-"+NeutrinoTypeLabel[neutype]+"-N-"+IntTypeLabel[IT]+"-"+pdfname+"_"+PDFVarLabel[pdfvar]+".dat";
        std::string filename_sigma = static_cast<std::string>(SAVE_PATH)+"M_"+mass_string+"MeV/sigma-"+NeutrinoTypeLabel[neutype]+"-N-"+IntTypeLabel[IT]+"-"+pdfname+"_"+PDFVarLabel[pdfvar]+".dat";

        std::cout << "Filename: " << filename_dsdxdy << std::endl;

        //ofstream outputfile_dsdy(filename_dsdy.c_str());
        ofstream outputfile_sigma(filename_sigma.c_str());

        for (double logenu=0.;logenu<=7.;logenu+=0.05){
          double enu = pow(10, logenu);
          xs_obj.Set_Neutrino_Energy(enu*pc->GeV);
          //for (double logx=-5.;logx<0.;logx+=0.025){
          //  double x = pow(10, logx);
          //  // here x is really y :).
          //  //double dsigdy = xs_obj.dsdy(x,PDFVarIndex[pdfvar])/cm2;
          //  //double dsigdy = xs_obj.dsdy(x)/cm2;
          //  double dsigdy = xs_obj.dsdy(x)/m2;
          //  outputfile_dsdy << enu << "\t"<< x << "\t" << dsigdy << std::endl;
          //}
          double sigma = xs_obj.total();
          //outputfile_sigma << enu << "\t"<< sigma/cm2 << std::endl;
          outputfile_sigma << enu << "\t"<< sigma/m2 << std::endl;
        }

        //outputfile_dsdy.close();
        outputfile_sigma.close();
      }
    }
  }

  return 0;
}
