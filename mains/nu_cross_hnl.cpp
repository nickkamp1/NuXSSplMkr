#include "lhapdf_cross_section.h"
#include <filesystem>

//#define SINGLE_ENERGY_TEST

int main(int argc, char* argv[]){

  PhysConst * pc = new PhysConst();

  // Get arguments
  std::string pdfname = (string) argv[1];
  double mass_double = atof(argv[2])/1000.;  // mass must be given in MeV
  std::stringstream ss_bool(argv[3]);
  bool is_hnl;
  ss_bool >> std::boolalpha >> is_hnl;
  std::string SAVE_PATH = (string) argv[4];
  std::ostringstream ss;
  ss << std::setw(7) << std::setfill('0') << (string) argv[2];
  std::string mass_string = ss.str();
  // std::string filename = SAVE_PATH+"M"+mass_string+"MeV/";
  // std::cout << "Creating " << filename << std::endl;
  // std::filesystem::create_directory(filename);
  LHAXS xs_obj(pdfname);
  //LHAXS xs_obj_HNL(pdfname);

  //enum IntType {CC,NC};
  enum NeutrinoType {neutrino,antineutrino};
  enum PDFVar {central,minus,plus};

  std::map<Current,std::string> IntTypeLabel {{CC,"cc"},{NC,"nc"},{EM,"em"}};
  std::map<NeutrinoType,double> CP_factor {{neutrino,1.},{antineutrino,-1}};
  std::map<NeutrinoType,std::string> NeutrinoTypeLabel {{neutrino,"nu"},{antineutrino,"nubar"}};
  //std::map<NeutrinoType,std::string> NeutrinoTypeLabel {{neutrino,"nutau"},{antineutrino,"nutaubar"}};
  std::map<PDFVar,int> PDFVarIndex {{central,0},{minus,-1},{plus,1}};
  std::map<PDFVar,std::string> PDFVarLabel {{central,"central"},{minus,"minus"},{plus,"plus"}};

  // // HNL mass
  std::cout << "HNL mass: " << mass_double << std::endl;
  xs_obj.Set_M_Lepton(mass_double*xs_obj.pc->GeV);
  std::cout << "Threshold: " << xs_obj.Threshold() << std::endl;
  // set bool to use custom cross section
  xs_obj.Set_IS_HNL(is_hnl);



  double cm2 = SQ(pc->cm);
  double m2 = SQ(pc->meter);


  for (Current IT : {NC,EM}) {
    std::cout << "Interaction Type: " << IT << std::endl;
    xs_obj.Set_InteractionType(IT);
    for (NeutrinoType neutype : {neutrino,antineutrino}){
      std::cout << "Neutrino Type: " << neutype << std::endl;
      xs_obj.Set_CP_factor(CP_factor[neutype]);
      for (PDFVar pdfvar : {central}){
        std::string filename_dsdxdy = static_cast<std::string>(SAVE_PATH)+"M_"+mass_string+"MeV/dsdxdy-"+NeutrinoTypeLabel[neutype]+"-N-"+IntTypeLabel[IT]+"-"+pdfname+"_"+PDFVarLabel[pdfvar]+"_v2.dat";
        std::string filename_sigma = static_cast<std::string>(SAVE_PATH)+"M_"+mass_string+"MeV/sigma-"+NeutrinoTypeLabel[neutype]+"-N-"+IntTypeLabel[IT]+"-"+pdfname+"_"+PDFVarLabel[pdfvar]+"_v2.dat";

        std::cout << "Diff filename: " << filename_dsdxdy << std::endl;
        std::cout << "Total filename: " << filename_sigma << std::endl;


        ofstream outputfile_dsdxdy(filename_dsdxdy.c_str());
        ofstream outputfile_sigma(filename_sigma.c_str());
        std::cout << "Threshold: " << xs_obj.Threshold() << std::endl;

        double Emin = std::max(10.,1.2*xs_obj.Threshold());
        for (double logenu=std::log10(Emin);logenu<=6.;logenu+=0.05){
          double enu = pow(10, logenu);
          std::cout << enu << std::endl;
          xs_obj.Set_Neutrino_Energy(enu*xs_obj.pc->GeV);
          double sigma = xs_obj.total(); // this is in eV^-2 if CC/NC, dimensionless if EM
          if (IT == NC || IT == CC) sigma /= m2;
          outputfile_sigma << enu << "\t" << sigma << std::endl;
          for (double logx=-5.;logx<0.;logx+=0.025){
            double x = pow(10, logx);
            for (double logy=-5.;logy<0.;logy+=0.025){
                double y = pow(10, logy);
                double zz[2];
                zz[0] = log(x);
                zz[1] = log(y);

                double dsigdxdy = xs_obj.KernelXS(zz);
                if (IT == NC || IT == CC) dsigdxdy /= m2; // same units as sigma
                // NOW WE NEED TO REMOVE THE JACOBIAN
                dsigdxdy /= (x*y); // convert from d/dlnx d/dlny to d/dx dy
                outputfile_dsdxdy << enu << "\t"<< x <<  "\t" << y << "\t" << dsigdxdy << std::endl;
            }
          }
        }

        outputfile_dsdxdy.close();
        outputfile_sigma.close();
      }
    }
  }

  return 0;
}
