/*
  Copyright 2018 Andreas Thune.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_DICTREAD_HEADER_INCLUDED
#define OPM_DICTREAD_HEADER_INCLUDED

#endif // OPM_DICTREAD_HEADER_INCLUDED

class DictRead
{
public:
  
  std::map<int,std::string> dict;

  void read_file_and_update(char* innit)
  {
    std::ifstream myfile {innit};
    
    std::string b;
    int a;
    while (myfile>> a>>b)
      {
	//std::cout << a <<" "<<b<<std::endl;
	dict[a]=b;
      }

    //std::cout <<"PEST data: "<< dict[3].data() <<dict[3] <<std::endl;
    
  }
  void write_param()
  {
    for (int i=0; i<dict.size(); ++i)
      {
	std::cout <<i<<": "<<dict[i] <<" ";
      }
    std::cout<<std::endl;
  }

  DictRead ()
  {
    dict[0]  = std::string("1");    //Use trans
    dict[1]  = std::string("0.01"); //tolerance
    dict[2]  = std::string("500");  //Max iter
    dict[3]  = std::string("1.01"); //IMBALANCE_TOL
    dict[4]  = std::string("1");    //Use wells
    dict[5]  = std::string("2");    //Zoltan DEBUG_LEVEL
    dict[6]  = std::string("0");    //Compare par sol to seq sol
    dict[7]  = std::string("0");    //Use node weights 
    dict[8]  = std::string("0");    //Reorder strategy (0=no,1=CMK,2=RCM,3=ghostLast)
    dict[9]  = std::string("0");    //Include ghost rows off-diag adjecency
    dict[10] = std::string("10");   //Number of timing iterations
    dict[11] = std::string("1.0");  //base exp for log weights.
    dict[12] = std::string("/global/D1/homes/andreast/linear_systems/patameters/cpr/fieldCase1.json");
    /*
    dict[7]  = std::string("2");    //DEBUG_LEVEL
    dict[8]  = std::string("10");   //PHG_REFINEMENT_LOOP_LIMIT
    dict[9]  = std::string("None"); //File with edge-weights (d[10]=1)
    dict[10] = std::string("0");    //Use edge-weight file (0=nor,1=read,2=uni,3=log)
    dict[11] = std::string("10.0"); //scaler in trans transform (d[10]=3)
    */
  }
};
