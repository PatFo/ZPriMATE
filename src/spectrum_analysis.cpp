#include "spectrum_analysis.h"
#include <stdio.h>
#include <fstream>
#include <cstdio>





//*******************************************************//
//            Class for scanning spectrum                //
//*******************************************************//


//Constructor
pheno::SpectrumScanner::SpectrumScanner(pheno::zpmodel* pmod, pheno::HadronXSec* phsec)
{
  //Initialize
  _model=pmod;
  _hsec=phsec;
  
  //Create default sampling region
  samplingRegions.clear(); 
  samplingRegions.push_back( new double[3]  );
  double max = 1.5*_model->mzp_();
  //Fill array with low, high, step
  samplingRegions.back()[0]=0;        //lower bound
  samplingRegions.back()[1]=max;      //upper bound
  samplingRegions.back()[2]=max/200;  //sampling step
  //Set default flag
  is_default=true;
}


void pheno::SpectrumScanner::reset_regions()
{
  //Delete all pointers in the container samplingRegions
  for(std::vector<double*>::iterator it=samplingRegions.begin(); it!=samplingRegions.end(); ++it)
  {
//     std::printf("Vector is: %g %g %g\n",(*it)[0],(*it)[1],(*it)[2]); //DEBUG
    delete[]  *it;  //Free all array memory
    *it = NULL;     //Make sure that freed memory block is not addressed any more
  }
  samplingRegions.clear(); //Delete all vector elements
}



//Destructor
pheno::SpectrumScanner::~SpectrumScanner()
{
  reset_regions();
}




//Method which adds a new sampling range to the list
void pheno::SpectrumScanner::set_interval(double low, double high, double step)
{
  //If is default, clear regions
  if(is_default)
  {
    reset_regions();
    is_default=false;
  }
  samplingRegions.push_back(new double[3]);
  samplingRegions.back()[0]=low;
  samplingRegions.back()[1]=high;
  samplingRegions.back()[2]=step;
}



void pheno::SpectrumScanner::sampler(char* outfile, double low, double high, double step)
{
  std::ofstream outf(outfile, std::ofstream::app);
  for(double E=low; E<high; E+=step)
  {
    std::vector<double> res;
    _hsec->crossSections(E, &res);
    outf<<E<<"\t\t"<<res[0]<<"\t\t"<<res[3]<"\n";
  }
}



/*
//Sampling the cross section
void pheno::SpectrumScanner::scan(char* outfile)
{
  //Define sampling for Z-peak
  float low(5), high(200);
  float step = (high-low)/30;
  sampler(outfile, low, high, step);
  
  //Define sampling between peaks
  double zpspread= _model->wzp_()*10;
  low = high;
  high = _model->mzp_() -zpspread;
  if(low<high) //Make sure lower bound is smaller than upper 
  {
    step =50;
    sampler(outfile, low, high, step);
  }
  low = high;
  high = _model->mzp_() + zpspread;
}*/
