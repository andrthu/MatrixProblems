template<class Mat>
double infoNNZ(Mat& t)
{
    std::vector<double> vals;
    double max = 0;
    double s = 0;
    for (auto row=t.begin();row!=t.end();++row) {
	auto col = row->begin();
	for (; col!=row->end(); ++col) {
	    auto val = *col;
	    if (val > max)
		max = val;
	    vals.push_back(val);
	    s += val;
	}
    }

    std::cout << max << " "<< vals.size()<<" "<<s/vals.size() <<std::endl;

    int num90 = 0;
    int num50 = 0;
    int num10 = 0;
    int num1 = 0;
    int num01 = 0;
    int num001 = 0;
    int num0 = 0;
    int numA = 0;
    for (int i = 0; i < vals.size(); ++i) {

	auto v = vals[i];
	if (v>0.9*max) {
	    num90++;
	    
	}
	if (v>0.5*max) {
	    num50++;
	    
	}
	if (v>0.1*max) {
	    num10++;
	    
	}
	if (v>0.01*max) {
	    num1++;
	    
	}
	if (v>0.001*max) {
	    num01++;
	    
	}
	if (v>0.0001*max) {
	    num001++;
	}
	if (v == 0)
	    num0++;
	if (v>s/vals.size())
	    numA++;
    }

    std::cout << num90 << " "<< num50 <<" "<< num10 <<" "<< num1 <<" "<< num01<<" "<<num001 <<
	" "<< num0<< " " << numA << std::endl;

    return max;
}

template<class Mat>
std::vector<double> infoNNZ2(Mat& t)
{
    std::vector<double> vals;
    double max = 0;
    double s = 0;
    for (auto row=t.begin();row!=t.end();++row) {
	auto col = row->begin();
	for (; col!=row->end(); ++col) {
	    auto val = *col;
	    if (val > max)
		max = val;
	    vals.push_back(val);
	    s += val;
	}
    }

    std::sort(vals.begin(),vals.end());

    return vals;
}

template<class Mat, class T>
void removeSmallTransNNZ(const Mat& A, Mat& M, T trans, double w)
{
    Dune::MatrixIndexSet op;
    op.resize( A.N(), A.N() );

    for (auto row = A.begin(); row != A.end(); ++row) {

	int d = row.index();
	op.add(d,d);
	auto col = row->begin();
	for (; col != row->end(); ++col) {
	    int nab = col.index();
	    if (trans.exists(d,nab)) {
		if (trans[d][nab] > w)
		    op.add(d,nab);
	    } else {
		op.add(d,nab);
	    }
	}
	
    }

    op.exportIdx(M);

    for (auto row = M.begin(); row != M.end(); ++row) {
	int d = row.index();
	M[d][d] = A[d][d];
	auto col = row->begin();
	for (; col != row->end(); ++col) {
	    int nab = col.index();
	    M[d][nab] = A[d][nab];
	}
    }
}
