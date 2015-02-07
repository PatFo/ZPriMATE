#include "read_config.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <stdexcept>
#include <cstring>


const int MAX_CHARS(100);
const int MAX_ITEMS(10);
const char* const DELIMITER = " $=";



int split_line(const char** itemlist,  char* line)
///Split a line into ites components separated by the DELIMITERS
{
  itemlist[0] = std::strtok(line, DELIMITER); //Split the string
  
  int len=1;
  for (int n = 1; n < MAX_ITEMS; n++)
  {
    // 'NULL' means continue splitting after last successful split
    itemlist[n] = std::strtok(NULL, DELIMITER); 
    if (!itemlist[n]) break; // no more tokens
    ++len;
  }
  
  return len;
}





void extract_config(dict &config, std::ifstream &istr)
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
	throw std::runtime_error("ERROR: More then one item as class definition in conifg file.");
      
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
    throw std::runtime_error("ERROR: Cannot open config file.");  
  std::cout<<"Opened config file "<<filename<<std::endl;
  
  //Set up the dictionary with initial parameters
  extract_config(config, ifs);
  
  ifs.close();
}



dict conf_reader::get_config()
{
  return config;
}
