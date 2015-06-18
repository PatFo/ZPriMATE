#include <config.h>
#include <cstring>
#include <errno.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#if MAC_SYS
        #include <mach-o/dyld.h>
#endif

#include "input.h"



const int MAX_CHARS(100);
const int MAX_ITEMS(30);
const char* const DELIMITERS1 = " $=";
const char* const DELIMITERS2 = " =";
const char* const SEPARATOR = "/";





  //**************************************//
  //         CONFIG READER                //  
  //**************************************//



int split_line(const char** itemlist,  char* line , const char* const DELIMS)
///Split a line into its components separated by the DELIMS
{
//   itemlist[0] = std::strtok(line, DELIMS); //Split the string
//   std::cout<<itemlist[0]<<std::endl;  // ###################################### DEBUG ####################### 
  
 
    char *token = std::strtok(line, DELIMS); //Split the string
    if(token ==NULL) return 0; //Abort if line is empty

    itemlist[0] = token;
    int len=1;
    while (token != NULL) {
//         std::cout << token << '\n';   // ###################################### DEBUG ####################### 
        token = std::strtok(NULL, DELIMS); 
        if(token ==NULL) break; //Abort if no more tokens
        itemlist[len] = token;
        ++len;
    }
//   for (int n = 1; n < MAX_ITEMS; n++)
//   {
//     
//     // 'NULL' means continue splitting after last successful split
//     char * curr_token = std::strtok(NULL, DELIMS); 
//     if (!curr_token) break;
//     itemlist[n]=curr_token;
// //     itemlist[n] = std::strtok(NULL, DELIMS); 
// //     if (!itemlist[n]) break; // no more tokens
// //     std::cout<<itemlist[n]<<std::endl; // ###################################### DEBUG ####################### 
//     ++len;
//   }
//   std::cout<<"length "<<len<<std::endl; // ###################################### DEBUG ####################### 
  return len;
}






void extract_config(dict &config, std::ifstream &istr)
///Read config file and set up dictionary for initialization of model
{
  while(!istr.eof())
  {
    char buf[MAX_CHARS];        
    istr.getline(buf,MAX_CHARS);

    //Sear for class beginning by $
    if(buf[0]=='$')
    {      
      const char* classes[MAX_ITEMS];
      
      int len = split_line(classes, buf, DELIMITERS1);      
      if(len!=1)
	throw std::runtime_error("ERROR: More than one item as class definition in conifg file.");
      
      std::string clas(classes[0]);
      parmap itemmap; 
      
      while(buf[0]!=0)
      {
	istr.getline(buf,MAX_CHARS);

	const char* items[MAX_ITEMS];
	if(buf[0]==0)break;
	else len = split_line(items, buf, DELIMITERS1);
	
	if(len==2)
	{
	  std::string key(items[0]);
	  double value = atof(items[1]);
	  itemmap[key]=value;
	}
      }      
      //Fill the dictionary with the class and corresponding item map
      config[clas]=itemmap;
    }
  }   
}






conf_reader::conf_reader(const char* filename)
//Implementation of constructor for reader
{
  std::ifstream ifs(filename);
  
  //Check whether specified file is readable. Else raise error and exit.
  if (!ifs.good())
    throw std::runtime_error("ERROR: Cannot open config file. Make sure that a config file is specified on startup.");  
  
  //Set up the dictionary with initial parameters
  extract_config(config, ifs);
  
  ifs.close();
}





//Returns the config dictionarry that has been set up
dict conf_reader::get_config()
{
  return config;
}





  //**************************************//
  //         START FILE READER            //
  //**************************************//



///Fill the string map with the variable and correponding values given in start file
void fill_map(strmap * pmap, std::ifstream * pistream)
{
  while(!pistream->eof())
  {
    char buffer[MAX_CHARS];        
    pistream->getline(buffer,MAX_CHARS);
    const char* words[MAX_ITEMS];
    int len = split_line(words, buffer, DELIMITERS2 );
    //Assert that line is not empty
    if(words[0]!=NULL)
    {
      //Check whether first item is a variable name
      if( words[0][0] =='$')
      {
        std::vector<std::string> items;
        for(int i=1; i < len; ++i)
        {
          items.push_back(std::string(words[i]));
        }
        pmap->operator[](words[0]) = items ;
      }
    }
  }
}



///Get the absolute path to the package root directory
std::string package_root_path()
{
  char buf[1024];
#if MAC_SYS
  char buf_tmp[1024];
  uint32_t size = sizeof(buf_tmp);
  _NSGetExecutablePath(buf_tmp, &size);
  realpath(buf_tmp, buf);
#else
  //Get full path of running C++ executable under LINUX
  int linkcount=readlink("/proc/self/exe", buf, sizeof(buf)-1);  //#####Something is going wrong a this point or after (Uninitialised value was created by a stack allocation) ## FIX
//   std::cout<<"Linksize "<<linkcount<<"\n"; //##################################################### DEBUG
#endif  
  //Split the path at the separator '/'
  std::string buf_str(buf);
  const char* words[MAX_ITEMS];
  std::printf( "Buffered path is %s \n", buf_str.c_str());        //##################################################### DEBUG
  int len = split_line(words, (char *)buf_str.c_str(), SEPARATOR );
  //Strip of the last 2 items of path (exec_dir/exec_name) to get package base path
  std::string path;
  for(int i=0; i<len-2; ++i)
  {
    path.append("/");
    path.append( words[i]);
  }
  return path;
}



///Returns absloute path; 'path' MUST either be absolute or relative to 'base'
std::string abs_path(std::string path, std::string base)
{
  std::string sep(SEPARATOR);
  //If path doesn't begin with '/' prepend base
  if (path[0]!=sep[0])
  {
    base.append("/");
    base.append(path);
    return base;
  //Else return path
  }else{
    return path;
  }
}








///Read input file and extrac settings
settings::settings(const char* startfile)
{
  std::printf( "Reading from input file %s ...\n", startfile);       
  
  std::ifstream input(startfile);
  //Check whether specified file is readable. Else raise error and exit.
  if (!input.good())
    throw std::runtime_error("ERROR: Cannot open settings file. Make sure that file is specified on startup."); 
  strmap inmap;
  //Extract settings strings
  fill_map(&inmap, &input);
  //Get base path to program dir
  std::string basepath = package_root_path();    
  
  
  ///Get VERBOSITY mode
  strmap::iterator it = inmap.find("$VERBOSE");
  if(it == inmap.end()) _verb = false;
  else _verb = atoi(((it->second)[0]).c_str());
  
  
  ///Get MODEL file 
  it= inmap.find("$MODEL");
  //Make sure that a file was specified
  if(it == inmap.end()) throw std::runtime_error("ERROR: No MODEL file  specified.\n\n\t Define in input file as: \'$MODEL = /absolute/path/to/file\'\n\t OR: \'$MODEL = <Mass of Z\'_SSM>\'\n"); 
  else
  {
    try
    {
      _mzssm = atof( ((it->second)[0]).c_str());
      //If file path is given, no conversion possible and value is set to 0
      if(_mzssm==0) throw 1;
      _use_ssm = true;
    }
    catch(int i)
    {
      _mzssm = -1;
      _use_ssm = false;
      _model = abs_path((it->second)[0], basepath);
    }
  } 

  
  ///Get PDF set 
  it= inmap.find("$PDF");
  //Make sure that a file was specified
  if(it == inmap.end())
  {
    _pdfset = basepath;
    _pdfset.append("/mstw/grid/mstw2008lo.00.dat");  // DEFAULT Pdf set
  }
  else  _pdfset = abs_path((it->second)[0], basepath);  
  
  
  ///Get LIMIT directory 
  it = inmap.find("$LIMITS");
  if(it == inmap.end()) throw std::runtime_error("ERROR: No LIMITS directory specified.\n\n\t Define in input file as: \'$LIMITS = /absolute/path/to/dir\'\n"); 
  else  _limdir = abs_path((it->second)[0], basepath);
  //Check whether _limdir exists
  struct stat st;
  if(stat(_limdir.c_str(),&st) != 0)
  {
    std::string error_msg("Limits directory doesn't exist: ");
    error_msg.append(limdir());
    throw std::runtime_error(error_msg.c_str()); 
  }
  
  
  ///Get OUTPUT directory 
  it = inmap.find("$ODIR");
  //Check whether a directory was specified; otherwise defaults to $ZPriMATE/results
  if(it == inmap.end()) 
  {
//     _odir = getenv("HOME");
    _odir = basepath;
    _odir.append("/results");    
    if(_verb) std::printf("INFO: No output directoy specified. Defaults to %s\n\n\t Else define on input as: \'$ODIR = /absolute/path/to/dir\'\n\n", _odir.c_str());
  }
  else  _odir = abs_path((it->second)[0], basepath);   
  //Check whether _odir exists; if not mkdir
  if(stat(_odir.c_str(),&st) != 0)
  {
    int status = mkdir(_odir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); 
    if(status == 0 && _verb)  std::printf("Created output directory %s \n", _odir.c_str());
    else if (status==-1)
    {
      std::printf("ERROR: Could not create directoy %s \n", _odir.c_str());
      throw std::runtime_error("ERROR: Output directory not created\n");
    }
  }
  
  
  ///Get BINNING file 
  it = inmap.find("$BINS");
  //Make sure that a file was specified
  if(it == inmap.end()) throw std::runtime_error("ERROR: No BINNING specified.\n\n\t Define in input file as: \'$BINS = /absolute/path/to/file\'\n\n"); 
  else  _binning = abs_path((it->second)[0], basepath); 
  //Check whether the binning file exists and is readable
  int bin_stat = access(_binning.c_str(), R_OK);
  if (bin_stat < 0) {
    if (errno == ENOENT) {
      std::string error_msg("Binning file doesn't exist: ");
      error_msg.append(_binning);
      throw std::runtime_error(error_msg.c_str()); 
    } else if (errno == EACCES) {
      std::string error_msg("Binning file is not readable: ");
      error_msg.append(_binning);
      throw std::runtime_error(error_msg.c_str());
    } else {
      std::string error_msg("Cannot access binning file: ");
      error_msg.append(_binning);
      throw std::runtime_error(error_msg.c_str());
    }
  }
  
  
  ///Get EFFICIENCIES file 
  it = inmap.find("$EFFICIENCIES");
  //Make sure that a file was specified
  if(it == inmap.end())
  {
    if(_verb) std::printf("INFO: No efficiency file specified.\n\n\t Define in input file as: \'$EFFICIENCIES = /absolute/path/to/file\'\n\n");
    _efficiencies = "";
  }
  else  _efficiencies = abs_path((it->second)[0], basepath);     
     
  
  
  ///Get Collider energy; if not specified defaults to 8TeV 
  it = inmap.find("$EBEAM");
  //Check whether beam energy was specified 
  if(it == inmap.end()) 
  {
    _ebeam= 8000;  // Default value
    if(_verb) std::printf("INFO: No collider energy specified. Defaults to %g GeV.\n\n\t Define in input file as: \'$EBEAM = /absolute/path/to/file\'\n\n", _ebeam);
  }
  else  _ebeam = atof( ((it->second)[0]).c_str());

  
  
  ///Get process code 
  it = inmap.find("$PROC");
  //Check whether process was specified 
  if(it == inmap.end())  throw std::runtime_error("ERROR: No process id specified.\n\n\t Define in input file as: \'$PROC = <id>\'\n\n"); 
  else _proc_id = atoi( ((it->second)[0]).c_str());
  
  
  
  ///Get luminosity
  it = inmap.find("$LUM");
  //Check whether luminosity was specified; if not set to 1 -> get cross section in [fb]
  if(it == inmap.end())  
  {
    if(_verb) std::printf("WARNING: No luminostiy specified. Obtain cross section in [fb].\n\n\t Else define on input as: \'$LUM = <luminosity in fb^-1>\'\n\n");
    _luminosity =1;
  }
  else _luminosity = atof( ((it->second)[0]).c_str());
  
  
  ///Get accuracy for numerical integration; defaults to 1e-3 (for cubature)
  it = inmap.find("$ACC");
  //Check whether luminosity was specified; if not set to 1 -> get cross section in [fb]
  if(it == inmap.end())  
  {
    _acc = 1e-2;
    if(_verb) std::printf("INFO: No accuracy for numerical integration specified. Defaults to %g.\n\n\t Else define on input as: \'$ACC = <accuracy>\'\n\n", _acc);
  }
  else  _acc = atof( ((it->second)[0]).c_str());
  
}




///FUNTIONS to access settings


bool settings::use_ssm()
///Returns false if a model file is specified
{
  return _use_ssm;
}


bool settings::verbose()
///Return verbosity mode
{
  return _verb;
}


int settings::proc_id()
///Returns the process code
{
  return _proc_id;
}


double settings::int_acc()
///Returns accuracy for numerical integration
{
  return _acc;
}


double settings::ebeam()
///Returns the beam energy
{
  return _ebeam;
}


double settings::mzssm()
///Returns the mass of Z'_SSM; -1 means not used
{
  return _mzssm;
}


double settings::luminosity()
///Returns the luminosity in [fb^-1]
{
  return _luminosity;
}


std::string settings::binning()
///Returns the path to the binning file
{
  return _binning;
}


std::string settings::efficiencies()
///Returns the path to the efficiencies file
{
  return _efficiencies;
}


std::string settings::limdir()
///Returns the path to the limits directoy
{
  return _limdir;
}


std::string settings::model()
///Returns the path to the model file
{
  return _model;
}


std::string settings::odir()
///Returns the output directory 
{
  return _odir;
}


std::string settings::pdfset()
///Returns the path to the pdfset
{
  return _pdfset;
}