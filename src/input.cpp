#include "input.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <stdexcept>
#include <cstring>


const int MAX_CHARS(100);
const int MAX_ITEMS(10);
const char* const DELIMITERS = " $=";






  //**************************************//
  //         CONFIG READER                //  
  //**************************************//



int split_line(const char** itemlist,  char* line)
///Split a line into its components separated by the DELIMITERS
{
  itemlist[0] = std::strtok(line, DELIMITERS); //Split the string
//   std::cout<<itemlist[0]<<std::endl;  // ###################################### DEBUG ####################### 
  
  int len=1;
  for (int n = 1; n < MAX_ITEMS; n++)
  {
    // 'NULL' means continue splitting after last successful split
    itemlist[n] = std::strtok(NULL, DELIMITERS); 
    if (!itemlist[n]) break; // no more tokens
//     std::cout<<itemlist[n]<<std::endl; // ###################################### DEBUG ####################### 
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
      
      int len = split_line(classes, buf);      
      if(len!=1)
	throw std::runtime_error("ERROR: More than one item as class definition in conifg file.");
      
      std::string clas(classes[0]);
      parmap itemmap; 
      
      while(buf[0]!=0)
      {
	istr.getline(buf,MAX_CHARS);

	const char* items[MAX_ITEMS];
	if(buf[0]==0)break;
	else len = split_line(items, buf);
	
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




void fill_map(inmap * pmap, std::ifstream * pistream)
{
    char buffer[MAX_CHARS];        
    pistream->getline(buffer,MAX_CHARS);
    
    if(buffer[0]=='$')
    {      
      const char* classes[MAX_ITEMS];
    }
}



settings::settings(const char* startfile)
{
  printf( "This is %s\n", startfile);       //###################################### DEBUG STATEMENT  
  std::ifstream input(startfile);
  
  //Check whether specified file is readable. Else raise error and exit.
  if (!input.good())
    throw std::runtime_error("ERROR: Cannot open settings file. Make sure that file is specified on startup."); 
  
}