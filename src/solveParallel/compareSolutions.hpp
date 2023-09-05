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

#ifndef OPM_COMPARESOLUTIONS_HEADER_INCLUDED
#define OPM_COMPARESOLUTIONS_HEADER_INCLUDED

#endif // OPM_COMPARESOLUTIONS_HEADER_INCLUDED

template<class Vec, class C>
void compareSolutions(Vec& x, Vec& x_s, C& cc, std::vector<int>& local2global,
		      std::vector<int>& overlapMap)
{
    typedef typename Vec::block_type V;
    int dim = V::dimension;
    
    double l1 = 0.0;
    double lInf = 0.0;
    int rank = cc.rank();

    V L1 = V(0.0);
    V LInf = V(0.0);

    for (int i = 0; i < local2global.size(); ++i)
    {
	int gid = local2global[i];
	if (overlapMap[gid]==0)
	{
	    auto diff = x[i] - x_s[gid];
	    
	    l1 += diff.one_norm();
	    
	    for (int j = 0; j < dim; ++j)
	    {
		L1[j] += std::abs(diff[j]);

		if ( LInf[j] < std::abs(diff[j]) )
		    LInf[j] = std::abs(diff[j]);
	    }

	    double idiff = diff.infinity_norm();
	    if (idiff > lInf)
		lInf = idiff;
	}
    }
    
    double l1_ = cc.sum(l1);
    double lInf_ = cc.max(lInf);

    V L1_ = cc.sum(L1);
    V LInf_ = V(0.0);//cc.max(LInf);
    for (int i = 0; i < dim; ++i)
	LInf_[i] = cc.max(LInf[i]);
    
    //l1_ = l1_/(dim*x.N());
    if (rank==0)
    {
	std::cout << "l1 error: " << l1_ << std::endl;
	std::cout << "lInf error: "<< lInf_ << std::endl;
	std::cout << "Seperated var L1 error: " << L1_ << std::endl;
	std::cout << "Seperated var lInf error: " <<LInf_ <<std::endl;
    }
}
