//  Created by Aviv Abramovich on 16/10/15.
//  Copyright (c) 2015 Aviv Abramovich. All rights reserved.
/*
 Grid_Multiplication.cpp
 This code generate a random matrix, encrypt them into encrypted matrices grid and multiply them.
 The test also measuring the time it take to each operation
 */

#include <sys/time.h>
#include <assert.h>

#include <vector>
#include <bitset>

#include "float.h"

#include "minlib/keys.h"
#include "minlib/encrypted_number.h"
#include "minlib/zp.h"
#include "minlib/settings.h"

#include "Grid.h"

/* Timers variables and methods*/
time_t time_begin, time_stop;
/*clock_t*/ long long clock_begin, clock_stop;

void resetTimers(string label="") {
    if(label.compare("")!=0)
        cout << label << endl;
    time(&time_begin);
    clock_begin = clock();
}

long long stopTimers(string label="") {
    time(&time_stop);
    clock_stop = clock(); //stop the clocks
    if(label.compare("")!=0)
        cout << "It took " << difftime(time_stop, time_begin) << " seconds and " << clock_stop-clock_begin<< " clock ticks " << label << endl;
    return clock_stop-clock_begin;
}

bool isPrime(long num) {
    if(num < 2)
        return false;
    for(unsigned int i=2; i*i <= num; i++) {
        if(num%i == 0)
            return false;
    }
    return true;
}

long power(long p, long r) {
    if(r ==0)
        return 1;
    long ret = power(p, r/2);
    return ret*ret*(r%2 ? p : 1);
}

void initiate_division (/*vector<vector<int > >&*/ long** division,int rows, int columns)
{
	for(int i=0;i<rows;i++)
	{
		division[i][0]=0;
	}

}

int increase_counter (int counter[],int size,int k)
{
	for(int i=0;i<size;i++)
	{
		if(counter[i]<k-1)
		{
			counter[i]++;
			return 1;
	
		}
		else
			counter[i]=0;
	}
return -1;//The counting is over

}

void  division_by_counter(long** division,int counter[],int counter_size)
{
	for(int i=0;i<counter_size;i++)
	{
		division[counter[i]][0]++;
		division[counter[i]][division[counter[i]][0]]=i;	
	}
}

double cal_distance(/*vector<vector<int > >&*/long**  division,int rows,int columns,PTMatrix* mat)
{
	double sum=0;
	long p,q;
double check =0;
	for(int i=0;i<rows;i++)
	{
		for(int j=1;j<=division[i][0];j++)
		{
			for(int k=j+1;k<=division[i][0];k++)
			{
				p=division[i][j];
				q=division[i][k];
				//check=((double)((mat->get_value(p,p))-2*(mat->get_value(p,q))+(mat->get_value(q,q))))/(double)division[i][0];//  (p^2-2pq+q^2)/|size of group|
				check=((double)((*mat)(p,p))-2*((*mat)(p,q))+((*mat)(q,q)))/(double)division[i][0];//  (p^2-2pq+q^2)/|size of group|
				sum+=check;
			}
		}
	}	

	return sum;
}

Ctxt encrypted_cal_distance(/*vector<vector<int > >&*/long**  division, int rows, int columns, const EncryptedMatrix &mat)
{
	EncryptedNumber sum;
	long p,q;
	EncryptedNumber check;

	ZP zp;

	for(int i=0;i<rows;i++)
	{
		for(int j=1;j<=division[i][0];j++)
		{
			for(int k=j+1;k<=division[i][0];k++)
			{
				p = division[i][j];
				q = division[i][k];
				//  (p^2-2pq+q^2)/|size of group|
				check = (EncryptedNumber(mat(p, p)) - EncryptedNumber(mat(p, q)) - EncryptedNumber(mat(p, q)) + EncryptedNumber(mat(q, q))) * zp.inv(division[i][0]);
				sum += check;
			}
		}
	}	

	return sum.number();
}

/*vector<double> k_means(PTMatrix* m,int k,int field)
{
	if(k<=0)
		return vector<double>();
	
	int division[k][m->getRows()+1];// the groups of vectors. row=group
	int counter[m->getRows()]; //it makes all the k^n option to divide n vectors to k groups
	vector<double> vec;//the vector of answers
	MatSize sqr_m_transpose(m->getColumns(),m->getRows());	
	PTMatrix* m_transpose = new PTMatrix(sqr_m_transpose,field);//transpose to m
	MatSize sqr_mat(m->getRows(),m->getRows());	
	PTMatrix* mat = new PTMatrix(sqr_mat,field);
	
	for(int i=0;i<m->getRows();i++) //build m_transpose
	{
		for(int j=0;j<m->getColumns();j++)
		{
			(*m_transpose)[j][i]=(*m)[i][j];
		}
	}

	for(int i=0;i<m->getRows();i++)//initialize counter
	{
		counter[i]=0;
	}
	initiate_division(division,k,m->getRows()+1);//initialize the 2D array
	*mat=((*m)*(*m_transpose));//calculate all p*q vector multiplications

	do
	{
		division_by_counter(division,counter,m->getRows());//make division by counter
		vec.push_back(cal_distance(division,k,m->getRows()+1,mat));//calculate the distance of specific division
	}
	while(increase_counter(counter,m->getRows(),k)!=-1);//all division options
	

	return vec;

}*/

double k_means(PTMatrix* m,int k,int field)
{
	if(k<=0)
		return -1;
	
	long** division = new long*[k*sizeof(long*)];// the groups of vectors. row=group
	for(int i=0; i<k; i++)
	{
		division[i] = new long[(m->getRows()+1)*sizeof(long)];
	}	
	int counter[m->getRows()]; //it makes all the k^n option to divide n vectors to k groups
	vector<double> vec;//the vector of answers
	MatSize sqr_m_transpose(m->getColumns(),m->getRows());	
	PTMatrix* m_transpose = new PTMatrix(sqr_m_transpose,2/*field*/);//transpose to m
	MatSize sqr_mat(m->getRows(),m->getRows());	
	PTMatrix* mat = new PTMatrix(sqr_mat,field);
	double min=DBL_MAX,temp;

	for(int i=0;i<m->getRows();i++) //build m_transpose
	{
		for(int j=0;j<m->getColumns();j++)
		{
			//cout<<"i="<<i<<" j="<<j<<" value="<<m->get_value(i,j)<<endl;
			//m_transpose->change_value(j,i,m->get_value(i,j));
			(*m_transpose)(j,i)=(*m)(i,j);
		}
	}
	for(int i=0;i<m->getRows();i++)//initialize counter
	{
		counter[i]=0;
	}

	/*std::cout << "Marix    m              ---> " << std::endl;
	m->print();
	std::cout << "Marix    m_transpose    ---> " << std::endl;
	m_transpose->print();*/
	initiate_division(division,k,m->getRows()+1);//initialize the 2D array
	*mat=((*m)*(*m_transpose));//calculate all p*q vector multiplications
	
	/*std::cout << "Marix    mat(result)    ---> " << std::endl;
         mat->print();*/
	
	do
	{

		division_by_counter(division,counter,m->getRows());//make division by counter
		temp=cal_distance(division,k,m->getRows()+1,mat);//calculate the distance of specific division
	
	//cout<<"temp = "<<temp<<endl;
		if(temp<min)
			min=temp;
		initiate_division(division,k,m->getRows()+1);//initialize the 2D array
	}
	while(increase_counter(counter,m->getRows(),k)!=-1);//all division options
	
	for(int i=0; i<k; i++)
	{
		delete(division[i]);
	}	

	delete(division);
	delete(m_transpose);
	delete(mat);

	return min;
}


Ctxt encrypted_k_means(EncryptedMatrix &m, int k, int field)
{
	assert(k > 0);
	
	long** division = new long*[k*sizeof(long*)];// the groups of vectors. row=group
	for(int i=0; i<k; i++)
	{
		division[i] = new long[(m.getRows()+1)*sizeof(long)];
	}	
	int counter[m.getRows()]; //it makes all the k^n option to divide n vectors to k groups
	vector<EncryptedNumber> vec;//the vector of answers
	MatSize sqr_m_transpose(m.getColumns(),m.getRows());	
	EncryptedMatrix* m_transpose = new EncryptedMatrix(sqr_m_transpose,2/*field*/);//transpose to m
	MatSize sqr_mat(m.getRows(),m.getRows());	
	EncryptedMatrix* mat;
	double min=DBL_MAX,temp;


	for(int i=0;i<m.getRows();i++)//initialize counter
	{
		counter[i]=0;
	}

	/*std::cout << "Marix    m              ---> " << std::endl;
	m.print();
	std::cout << "Marix    m_transpose    ---> " << std::endl;
	m_transpose->print();*/
	initiate_division(division,k,m.getRows()+1);//initialize the 2D array
	*mat=(m*(*m_transpose));//calculate all p*q vector multiplications
	
	/*std::cout << "Marix    mat(result)    ---> " << std::endl;
         mat->print();*/
	
	do
	{

		division_by_counter(division,counter,m.getRows());//make division by counter
		vec.push(EncryptedNumber(encrypted_cal_distance(division,k,m.getRows()+1,*mat)));//calculate the distance of specific division
	
	//cout<<"temp = "<<temp<<endl;
	//	if(temp<min)
	//		min=temp;
		initiate_division(division,k,m.getRows()+1);//initialize the 2D array
	}
	while(increase_counter(counter,m.getRows(),k)!=-1);//all division options
	
	for(int i=0; i<k; i++)
	{
		delete(division[i]);
	}	

	delete(division);
	delete(m_transpose);
	delete(mat);

	
	return doMin(vec);
}


int main(int, char **) {
	    long m, r, p,L, c, w, s, d, security, enc1, encMul, recommended;
	    long long EncSec,DecSec, enc, dec, ptMul,k_means_sec,k_means_ticks;
	    char tempChar;
	    bool toEncMult, toPrint;
	    int k;

	    
	    //Scan parameters
	    
	    cout << "Enter HElib's keys paramter. Enter zero for the recommended values" << endl;
	    
	    while(true) {
		cout << "Enter the field of the computations (a prime number): ";
		cin >> p;
		if(isPrime(p))
		    break;
		cout << "Error! p must be a prime number! " << endl;
	    }
	    while(true) {
		recommended = 1;
		cout << "Enter r (recommended " << recommended <<"): ";
		cin >> r;
		if(r == 0)
		    r = recommended;
		if(r > 0)
		    break;
		cout << "Error! r must be a positive number!" << endl;
	    }
	    while(true) {
		recommended = 16;
		cout << "Enter L (recommended " << recommended <<"): ";
		cin >> L;
		if(L == 0)
		    L = recommended;
		if(L > 1)
		    break;
		cout << "Error! L must be a positive number!" << endl;
	    }
	    while(true) {
		recommended = 3;
		cout << "Enter c (recommended " << recommended <<"): ";
		cin >> c;
		if(c == 0)
		    c = recommended;
		if(c > 1)
		    break;
		cout << "Error! c must be a positive number!" << endl;
	    }
	    while(true) {
		recommended = 64;
		cout << "Enter w (recommended " << recommended <<"): ";
		cin >> w;
		if(w == 0)
		    w = recommended;
		if(w > 1)
		    break;
		cout << "Error! w must be a positive number!" << endl;
	    }
	    while(true) {
		recommended = 0;
		cout << "Enter d (recommended " << recommended <<"): ";
		cin >> d;
		if(d >= 0)
		    break;
		cout << "Error! d must be a positive or zero!" << endl;
	    }
	    while(true) {
		recommended = 0;
		cout << "Enter s (recommended " << recommended <<"): ";
		cin >> s;
		if(s >= 0)
		    break;
		cout << "Error! s must be a positive or zero!" << endl;
	    }
	    while(true) {
		recommended = 128;
		cout << "Enter security (recommended " << recommended << "): ";
		cin >> security;
		if(security == 0)
		    security = recommended;
		if(security >= 1)
		    break;
		cout << "Error! security must be a positive number " << endl;
	    }
	    
	    ZZX G;
	    m = FindM(security,L,c,p, d, s, 0);
	    FHEcontext context(m, p, r);
	    // initialize context
	    buildModChain(context, L, c);
	    // modify the context, adding primes to the modulus chain
	    FHESecKey secretKey(context);
	    // construct a secret key structure
	    FHEPubKey& publicKey = secretKey;
	    // an "upcast": FHESecKey is a subclass of FHEPubKey
	    
	    //if(0 == d)
	    G = context.alMod.getFactorsOverZZ()[0];
	    
	    secretKey.GenSecKey(w);
	    // actually generate a secret key with Hamming weight w
	    
	    addSome1DMatrices(secretKey);
	    EncryptedArray ea(context, G);
	    // constuct an Encrypted array object ea that is
	    // associated with the given context and the polynomial G

		Keys::setKeys(&publicKey, &secretKey, &ea, &context);
		ZP::set_global_p(p); assert(r == 1);
		ZP zp;
		Settings<EncryptedNumber>::max_value(100, zp);
	    
	    long nslots = ea.size(), field = power(p,r);
	    cout << "nslots: " << nslots << endl ;
	    cout << "Computations will be modulo " << field << endl;
	    cout << "m: " << m << endl;
	    
	    unsigned int sz1,sz2;
	    while(true) {
		cout << "Enter number of vectors: ";
		cin >> sz1;
		if(sz1 > 1)
		    break;
		cout << "Error! the value must be a positive number!" << endl;
	    }
	    while(true) {
		cout << "Enter the dimension: ";
		cin >> sz2;
		if(sz2 > 1)
		    break;
		cout << "Error! the value must be a positive number!" << endl;
	    }
	    while(true) {
		cout << "Enter k the number of group ";
		cin >> k;
		if(k > 1)
		    break;
		cout << "Error! the value must be a positive number!" << endl;
	    }
	   
	    MatSize atom(nslots,nslots), sqr(sz1,sz2);
	    PTMatrix* mat = new PTMatrix(sqr, 2/*field*/);
	    cout<<"MAT1 SUCCESS"<<endl;

	enc=0;
//	dec=0;
        EncSec=0;
  //      DecSec = 0;
	
	//DecSec=dec/(long)CLOCKS_PER_SEC;
	//cout << "It took " << enc << " clock ticks and "<<EncSec<<" seconds to encrypt the  matrix" <<  endl;
	resetTimers();
	int ans = k_means(mat, k, field);
	k_means_ticks = stopTimers("to do regular k-means");
	k_means_sec = k_means_ticks/(long)CLOCKS_PER_SEC;
	cout << "ANSWER IS " << ans  <<  endl;

	resetTimers();
	EncryptedMatrix enc_mat = mat->encrypt(ea, publicKey);
	enc=stopTimers("to encrypt the matrix");
	EncSec=enc/(long)CLOCKS_PER_SEC;

	resetTimers();
	Ctxt ansENC = encrypted_k_means(enc_mat, k, field);
	k_means_ticks = stopTimers("to do encrypted k-means");
	k_means_sec = k_means_ticks/(long)CLOCKS_PER_SEC;
	vector<long> v = new vector<long>(1);	
	cout << "Decrypting answer .."<<endl;
	cout << "ANSWER IS " << ea.decrypt(ansENC,secretKey,v) <<  endl;

	//cout << "It took " << dec << " clock ticks and " <<DecSec<<" seconds to do k_means" <<  endl;
	//cout << "It took " << dec << " clock ticks and " <<DecSec<<" seconds to decrypt the result vector" <<  endl;
		
		MatSize vecSize(power(k,sz1/2),power(k,sz1/2));
    
   /* PTMatrix vec(vecSize, field);
    cout << "Encrypting vec .." << endl;
    EncryptedMatrix EncVec = vec.encrypt(ea,publicKey);
    cout << "Decrypting vec .."<<endl;
    long decVec = 0;
    resetTimers();
    EncVec.decrypt(ea,secretKey);
    decVec = stopTimers("to decrypt vec");*/



	    
	    cout << "\n\n----------------------------------------Summary------------------------------ " << endl;
	    cout << "p: " << p << ", r: " << r << ", L: " << L << ", c: " << c << ", w: " << w << ", d: " << d << ", s: " << s << ", security: " << 		    security << endl;
	    cout << "nslots: " << nslots << "\nm: " << m << endl;
	cout << "It took " << enc << " clock ticks and "<<EncSec<<" seconds to encrypt the matrix" <<  endl;
	cout << "It took " << k_means_ticks << " clock ticks and " <<k_means_sec<<" seconds to do k_means" <<  endl;
	//cout << "It took " << dec << " clock ticks and " <<DecSec<<" seconds to decrypt the result vector" <<  endl;


	delete(mat);
	    return 0;
}//End main()
