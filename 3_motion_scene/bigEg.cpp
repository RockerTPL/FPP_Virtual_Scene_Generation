//
//  bigEg.cpp
//  MatlabC
//
//  Created by 田沛林 on 2020/5/17.
//  Copyright © 2020 Andy. All rights reserved.
//

#include <iostream>
#include "mex.h"
#include "outil.h"

#define Num 1000
#define NUM 10000000
using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
//    *prhs = { ij, sizeij, ci(int), cj(int), gap(double) }
    double *ijp, *sizeijp, *cip, *cjp, *gapp;
    
    ijp = mxGetPr(prhs[0]);
    sizeijp = mxGetPr(prhs[1]);
    cip = mxGetPr(prhs[2]);
    cjp = mxGetPr(prhs[3]);
    gapp = mxGetPr(prhs[4]);
    
    int sizeij = *sizeijp, ci = *cip, cj = *cjp;
    double gap = *gapp;

    double ij[NUM][2];
    for(int k = 0; k < sizeij; k++)
    {
        ij[k][0] = ijp[k];
        ij[k][1] = ijp[sizeij + k];
    }



    // distance map (ci cj)
    int Dmp[Num][Num];
    for(int i = 0; i < ci; i++)
        for(int j = 0; j < cj; j++)
            Dmp[i][j] = 10000;
    for(int k = 0; k < sizeij; k++)
        Dmp [ int(ij[k][0] - 1) ] [ int(ij[k][1] - 1) ] = 0;


    for(int i = 1; i < ci; i++)
        for(int j = 1; j < cj; j++)
            Dmp[i][j] = min( Dmp[i][j], min(Dmp[i-1][j], Dmp[i][j-1]) + 1 );
    for(int j = cj-2; j >= 0; j--)
        Dmp[ci-1][j] = min( Dmp[ci-1][j], Dmp[ci-1][j+1] + 1 );
    for(int i = ci-2; i >= 0; i--)
        Dmp[i][cj-1] = min( Dmp[i][cj-1], Dmp[i+1][cj-1] + 1 );
    for(int i = ci-2; i >= 0; i--)
        for(int j = cj-2; j >= 0; j--)
            Dmp[i][j] = min( Dmp[i][j], min(Dmp[i+1][j], Dmp[i][j+1]) + 1 );


    for(int i = 0; i < ci; i++)
        for(int j = 0; j < cj; j++)
            if(Dmp[i][j] <= gap)
            {
                ij[sizeij][0] = i+1;
                ij[sizeij][1] = j+1;
                sizeij++;
            }
    
    
    plhs[0] = mxCreateDoubleMatrix(sizeij, 2, mxREAL);
    double* ijres;
    ijres = mxGetPr(plhs[0]);
    for(int k = 0; k < sizeij; k++)
    {
        ijres[k] = ij[k][0];
        ijres[sizeij + k] = ij[k][1];
    }
    
}
