#include "spectrum_analysis.h"
#include <stdio.h>
#include <fstream>
#include <cstdio>
#include <cstring>






//*******************************************************//
//            Class for scanning spectrum                //
//*******************************************************//


//Constructor
template<class CrossSection>
pheno::SpectrumScanner<CrossSection>::SpectrumScanner(pheno::ZpModel* pmod, CrossSection* phsec)
{
  //Initialize
  _model=pmod;
  _hsec=phsec;
  
  //Create default sampling region
  samplingRegions.clear(); 
  samplingRegions.push_back( new double[3]  );
  double max = 1.5*_model->mzp_(), min=5;
  //Fill array with low, high, step
  samplingRegions.back()[0]=min;        //lower bound
  samplingRegions.back()[1]=max;      //upper bound
  samplingRegions.back()[2]=(max-min)/200;  //sampling step
  //Set default flag
  is_default=true;
}


template<class CrossSection>
void pheno::SpectrumScanner<CrossSection>::reset_sampling()
{
  //Delete all pointers in the container samplingRegions
  for(sampling_scheme::iterator it=samplingRegions.begin(); it!=samplingRegions.end(); ++it)
  {
//     std::printf("Vector is: %g %g %g\n",(*it)[0],(*it)[1],(*it)[2]); //DEBUG
    delete[]  *it;  //Free all array memory
    *it = NULL;     //Make sure that freed memory block is not addressed any more
  }
  samplingRegions.clear(); //Delete all vector elements
}



//Destructor
template<class CrossSection>
pheno::SpectrumScanner<CrossSection>::~SpectrumScanner()
{
  reset_sampling();
}




//Method which adds a new sampling range to the list
template<class CrossSection>
void pheno::SpectrumScanner<CrossSection>::add_interval(double low, double high, double step)
{
  //If is default, clear regions
  if(is_default)
  {
    reset_sampling();
    is_default=false;
  }
  samplingRegions.push_back(new double[3]);
  samplingRegions.back()[0]=low;
  samplingRegions.back()[1]=high;
  samplingRegions.back()[2]=step;
}




//Copy interval array
template<class CrossSection>
void pheno::SpectrumScanner<CrossSection>::add_interval(double* pinterval)
{
  //If is default, clear regions
  if(is_default)
  {
    reset_sampling();
    is_default=false;
  }
  samplingRegions.push_back(new double[3]);
  samplingRegions.back()[0]=pinterval[0];
  samplingRegions.back()[1]=pinterval[1];
  samplingRegions.back()[2]=pinterval[2];
}





//Sample Cross Section in given interval and append data to output file
template<class CrossSection>
void pheno::SpectrumScanner<CrossSection>::sampler(char* outfile, double low, double high, double step)
{
  std::ofstream outf(outfile, std::ofstream::app);
  for(double E=low; E<high; E+=step)
  {
    std::vector<double> res;
    _hsec->crossSections(E, &res);
    char buffer[100];
    std::sprintf(buffer, "%-7g %-15g %-15g %-15g %-15g\n", E, res[0], res[1], res[2], res[3]);
    outf<<buffer;
  }
  outf.close();
}




//Sample the Cross Section over all sampling regions and write data to file
template<class CrossSection>
void pheno::SpectrumScanner<CrossSection>::scan(char* outfile)
{
  //Create file header
  std::ofstream outf(outfile);
  char buffer[100];
  std::sprintf(buffer, "%-7s %-15s %-15s %-15s %-15s\n", "E", "Xtot", "Xint", "Xsig", "Xsm");
  outf<<buffer;
  outf.close();
  
  //Loop over sampling regions
  for(sampling_scheme::iterator it=samplingRegions.begin(); it!=samplingRegions.end(); ++it)
  {
    sampler(outfile, (*it)[0], (*it)[1], (*it)[2]);
  }
  std::printf("Cross section scan written to %s\n", outfile);
}




//Explicit instantiation so that linking works --> ONLY these types can be used
template class pheno::SpectrumScanner<pheno::PartonXSec>;
template class pheno::SpectrumScanner<pheno::HadronXSec>;





//*******************************************************//
//            Function to get binning from file          //
//*******************************************************//



pheno::binning pheno::get_binning(char* binfile)
{
  pheno::binning tmp;
  std::ifstream ifs(binfile);
  const int MAX_CHARS(100);
  char buf[MAX_CHARS];     
  
  //Ignore first line: contains only comments
  ifs.getline(buf,MAX_CHARS);
  
  //Loop over data
  while(!ifs.eof())
  {
    //Read line
    ifs.getline(buf,MAX_CHARS);
    
    char* tokens[2];
    const char* const DELIMITER = "\t";
    
    //Split line into items
    tokens[0] = std::strtok(buf, DELIMITER);
    tokens[1] = std::strtok(NULL, DELIMITER);
//     std::printf("%s\t%s",tokens[0], tokens[1]);
    
    //Insert values into binning vector
    tmp.push_back( std::pair<double,double>(atof(tokens[0]), atof(tokens[1])) );
  }
  
  return tmp;
}

