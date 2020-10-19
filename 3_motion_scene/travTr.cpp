//
//  travTr.cpp
//  MatlabC
//
//  Created by ?°æ??? on 2020/5/5.
//  Copyright Â© 2020 Andy. All rights reserved.
//

#pragma clang diagnostic push
#pragma clang diagnostic ignored"-Wshorten-64-to-32"
#include <stdio.h>

#include <iostream>
#include <vector>
#include <math.h>
#include "mex.h"
#include "outil.h"
using namespace std;




void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
//    *prhs = {i, j, List(1 1 nTr), All(3 3 nTr), Cover(ci cj nTr), 
//             temp(2 1), Dir(3 1), crit, nTrp}
    // *plhs = {List(1 1 nTr)}
    
    // é¦???????
    double *ip, *jp, *Listp, *Allp, *Coverp, *tempp, *Dirp, *critp, *nTrp;
    ip = mxGetPr(prhs[0]);
    jp = mxGetPr(prhs[1]);
    Listp = mxGetPr(prhs[2]);
    Allp = mxGetPr(prhs[3]);
    Coverp = mxGetPr(prhs[4]);
    tempp = mxGetPr(prhs[5]);
    Dirp = mxGetPr(prhs[6]);
    critp = mxGetPr(prhs[7]);
    nTrp = mxGetPr(prhs[8]);
    
    int i = *ip, j = *jp, nTr = *nTrp;
    double crit = *critp;
    
    // get element in All and Cover     [idxC2M3D(i, j, k, r, c)]
    int *Alldim, *Coverdim;
    Alldim = (int*)mxGetDimensions(prhs[3]);
    Coverdim = (int*)mxGetDimensions(prhs[4]);
    int Allr = Alldim[0], Allc = Alldim[2];
    int Coverr = Coverdim[0], Coverc = Coverdim[2];
    
    
    
    int k;
    double Edg1[3], Edg2[3], Edg3[3], N[3], Dp, NDir, I[3], AI[3], BI[3], CI[3], cs1[3], cs2[3], cs3[3];
    
    for(int tr = 0; tr < nTr; tr++)
    {
        k = Listp[tr];
        if(k != 0)
        {
            // 3 points of triangle
            double A[3] = {Allp[idxM2C3D(1,1,k,Allr, Allc)],
                           Allp[idxM2C3D(2,1,k,Allr, Allc)],
                           Allp[idxM2C3D(3,1,k,Allr, Allc)]};
            double B[3] = {Allp[idxM2C3D(1,2,k,Allr, Allc)],
                           Allp[idxM2C3D(2,2,k,Allr, Allc)],
                           Allp[idxM2C3D(3,2,k,Allr, Allc)]};
            double C[3] = {Allp[idxM2C3D(1,3,k,Allr, Allc)],
                           Allp[idxM2C3D(2,3,k,Allr, Allc)],
                           Allp[idxM2C3D(3,3,k,Allr, Allc)]};
            
            // get edge of triangle
            vsub(B, A, Edg1);
            vsub(C, B, Edg2);
            vsub(A, C, Edg3);
            // norm vector of surface o triangle
            cross(Edg1, Edg2, N);
            normlz(N);
            // constant of surface equation
            Dp = (-1)*dot(N, A);
            NDir = dot(N, Dirp);
            if(NDir != 0) mul(Dirp, -Dp/NDir, I);  //abs(NDir-0) > abs(crit)
            else{
                I[0] = 0;
                I[1] = 0;
                I[2] = -1;
                Coverp[idxM2C3D(i, j, k, Coverr, Coverc)] = k;
            }
            
            if(I[2] > crit)
            {   
                vsub(I, A, AI);
                vsub(I, B, BI);
                vsub(I, C, CI);
                cross(Edg1, AI, cs1);
                cross(Edg2, BI, cs2);
                cross(Edg3, CI, cs3);
                if(I[2] < tempp[1] &&
                   dot(N, cs1) > crit &&
                   dot(N, cs2) > crit &&
                   dot(N, cs3) > crit)
                {
                    if(tempp[0] != 0)
                        Coverp[idxM2C3D(i, j, (int)tempp[0], Coverr, Coverc)] = tempp[0];
                    tempp[0] = k;
                    tempp[1] = I[2];
                }
            }
        }
    }
    if(tempp[0] != 0)
        Coverp[idxM2C3D(i, j, (int)tempp[0], Coverr, Coverc)] = tempp[0];
}

#pragma clang diagnostic pop
