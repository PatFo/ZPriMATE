#ifndef INPUT_H
#define INPUT_H

#include <map>
#include <string>
#include <vector>

// integer index for parameters
typedef std::map<unsigned int, double>  parmap;
typedef std::map<std::string, parmap >  dict;

typedef std::map<std::string, std::vector<std::string> > strmap;




  //**************************************//
  //         CONFIG READER                //  
  //**************************************//

class conf_reader{
  private:
    dict config;
  // Flags for input scheme
  void extract_config(dict &config, std::ifstream &istr);
  bool fullInput;
  // Flag to identify if there are universal couplings
  // e.g. QUARK, LEPTON, LEPTONR
  // QUARK|QUARKL|QUARKR|LEPTON|LEPTONL|LEPTONR
  unsigned int CouplingFlag;
  unsigned int getFlag(std::string blockName);
  unsigned int updateFlag(std::string blockName);
  void checkFlag(unsigned int flag);
  void addConfig(std::string, parmap);
  public:
  // Check if ALL couplings are diagonal and family universal
  bool couplingUniversal();

  // These methods do NOT check for diagonal or family universal but merely
  // offer a tool for a quick readout of the flag
  bool QuarkUni();
  bool QuarkLeftUni();
  bool QuarkRightUni();
  bool LeptonUni();
  bool LeptonLeftUni();
  bool LeptonRightUni();
  dict get_config();
  //Constructor
  conf_reader(const char* filename);
  
};





  //**************************************//
  //         START FILE READER            //
  //**************************************//
  
  
  class settings{
  private:
    bool _use_ssm;
    bool _verb;
    bool _force;
    int _proc_id; 
    double _acc;
    double _ebeam;
    double _luminosity;
    double _mzssm;
    std::string _binning;
    std::string _efficiencies;
    std::string _limdir;
    std::string _model;
    std::string _odir;
    std::string _pdfset;

    void initOptions();
    int parseCmdLine(int argc,char** argv);
    int getOptions(char ** begin, char ** end, const std::string & option);
    int getArguments(char ** begin, char ** end, const std::string & option);
    int setOption(char**,int index);
  public:
    // Options with operands, e.g. '-h'
    std::map<std::string, std::string> options;
    // Arguments of command line without operands
    // position 0: settingsFile
    // position 1: tmpFile (optional)
    std::vector<std::string> arguments;
    bool force();
    bool use_ssm();
    bool verbose();
    int proc_id();    
    double int_acc();
    double ebeam();
    double luminosity();
    double mzssm();
    std::string binning();
    std::string efficiencies();
    std::string limdir();
    std::string model();
    std::string odir();
    std::string pdfset();
    //Constructor
    settings(int argc,char** argv);
  };



#endif
