#include "input.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <stdexcept>
#include <cstring>


const int MAX_CHARS(100);
const int MAX_ITEMS(10);
const char* const DELIMITERS1 = " $=";
const char* const DELIMITERS2 = " =";





  //**************************************//
  //         CONFIG READER                //  
  //**************************************//



int split_line(const char** itemlist,  char* line , const char* const DELIMS)
///Split a line into its components separated by the DELIMS
{
  itemlist[0] = std::strtok(line, DELIMS); //Split the string
  std::cout<<itemlist[0]<<std::endl;  // ###################################### DEBUG ####################### 
  
  int len=1;
  for (int n = 1; n < MAX_ITEMS; n++)
  {
    // 'NULL' means continue splitting after last successful split
    itemlist[n] = std::strtok(NULL, DELIMS); 
    if (!itemlist[n]) break; // no more tokens
    std::cout<<itemlist[n]<<std::endl; // ###################################### DEBUG ####################### 
    ++len;
  }
  
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



///Fill the string map with the variable and correponding value given in start file
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
        printf("Loading data\n");

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



settings::settings(const char* startfile, bool verbose)
{
  printf( "This is %s\n", startfile);       //###################################### DEBUG STATEMENT  
  std::ifstream input(startfile);
  
  //Check whether specified file is readable. Else raise error and exit.
  if (!input.good())
    throw std::runtime_error("ERROR: Cannot open settings file. Make sure that file is specified on startup."); 
  //Extract settings strings
  strmap inmap;
  fill_map(&inmap, &input);
  
  
  //Get PDF set 
  strmap::iterator it= inmap.find("$PDF");
  //Make sure that a file was specified
  if(it == inmap.end()) throw std::runtime_error("ERROR: No PDF set specified.\n\n\t Define in input file as: \'$PDF = /absolute/path/to/file\'\n"); 
  else
  {
    _pdfset = (it->second)[0];
    printf("%s\n",_pdfset.c_str());         //###################################### DEBUG STATEMENT   
  }   
  
  //Get BINNING file 
  it = inmap.find("$BINS");
  //Make sure that a file was specified
  if(it == inmap.end()) throw std::runtime_error("ERROR: No BINNING specified.\n\n\t Define in input file as: \'$BINS = /absolute/path/to/file\'\n"); 
  else
  {
    _binning = (it->second)[0];
    printf("%s\n",_binning.c_str());         //###################################### DEBUG STATEMENT  
  }    
  
  
    //Get BINNING file 
  it = inmap.find("$SREGION");
  //Check whether search region was specified 
  if(it == inmap.end()) 
  {
    if(verbose) printf("WARNING: No search region has been specified. Full data set is used.\n");
    _smin = -1; 
    _smax = -1;
    
  }
  else
  {
    _smin = atof( ((it->second)[0]).c_str());
    _smax = atof( ((it->second)[1]).c_str());
    printf("%g\t%g\n",_smin, _smax);   
  }   

  printf("END of Constructor");
}