//
//  getEg.cpp
//  MatlabC
//
//  Created by ?°æ??? on 2020/5/17.
//  Copyright Â© 2020 Andy. All rights reserved.
//

#include <iostream>
#include "mex.h"
#include "outil.h"

#define Num 1000
#define NUM 100000000
using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
//    *prhs = { Alleg(6 nTr), K(3 3), nEg(int), ci(int), cj(int), crit(double), gap(double) }
    double *Alleg, *Kp, *nEgp, *cip, *cjp, *critp, *gapp;
    
    Alleg = mxGetPr(prhs[0]);
    Kp = mxGetPr(prhs[1]);
    
    nEgp = mxGetPr(prhs[2]);
    cip = mxGetPr(prhs[3]);
    cjp = mxGetPr(prhs[4]);
    critp = mxGetPr(prhs[5]);
    gapp = mxGetPr(prhs[6]);
    
    int nEg = *nEgp, ci = *cip, cj = *cjp;
    double crit = *critp, gap = *gapp;
    double K[9];
    for(int i = 0; i < 9; i++) K[i] = Kp[i];


    double pt1x, pt1y, pt1z, pt2x, pt2y, pt2z, ptw,
           pt1[3], pt2[3], ptdir[3],
           pt1uv0[3], pt1uv[2], pt2uv0[3], pt2uv[2],
           *idi0 = new double[NUM], *idj0 = new double[NUM], (*ij)[2] = new double[NUM][2];
    int gp, pt1uvr[2], pt2uvr[2],
        *idi = new int[NUM], *idj = new int[NUM],
        (*id0)[2] = new int[NUM][2], (*id1)[2] = new int[NUM][2];
    int sizeij = 0, sizeid1 = 0;

    for(int e = 0; e < nEg; e++)
    {
        pt1x = Alleg[6*e + 0];
        pt1y = Alleg[6*e + 1];
        pt1z = Alleg[6*e + 2];
        pt2x = Alleg[6*e + 3];
        pt2y = Alleg[6*e + 4];
        pt2z = Alleg[6*e + 5];
        pt1[0] = pt1x; pt1[1] = pt1y; pt1[2] = pt1z;
        pt2[0] = pt2x; pt2[1] = pt2y; pt2[2] = pt2z;


        if(pt1z > 0 && pt2z > 0)
        {
//             if(e==2)
//                 cout << e << " "<< pt1z << " " << pt2z << endl;
            mul(pt1, 1.0/pt1z, pt1);
            dotMV(K, pt1, pt1uv0);
            pt1uv[0] = pt1uv0[0]; pt1uv[1] = pt1uv0[1];
            roundVec(2, pt1uv, pt1uvr);

            mul(pt2, 1.0/pt2z, pt2);
            dotMV(K, pt2, pt2uv0);
            pt2uv[0] = pt2uv0[0]; pt2uv[1] = pt2uv0[1];
            roundVec(2, pt2uv, pt2uvr);

            gp = max(abs(pt1uv[1] - pt2uv[1]), abs(pt1uv[0] - pt2uv[0]));
            // idi??idjï¼??¿åº¦ gpï¼?id0ï¼?gp 2ï¼?
            linspace(pt1uv[1], pt2uv[1], gp, idi0);
            roundVec(gp, idi0, idi);
            linspace(pt1uv[0], pt2uv[0], gp, idj0);
            roundVec(gp, idj0, idj);

            for(int i = 0; i < gp; i++)
            {
                id0[i][0] = idi[i];
                id0[i][1] = idj[i];
            }

            sizeid1 = chooseInt(1, ci, 1, cj, gp, id0, id1);
            for(int i = 0; i < sizeid1; i++)
            {
                ij[sizeij+i][0] = id1[i][0];
                ij[sizeij+i][1] = id1[i][1];
//                 cout << id1[i][0] << " ";
//                 cout << id1[i][1] << endl;
            }
            sizeij = sizeij + sizeid1;
        }


        else if(pt1z > 0 && pt2z <= 0)
        {
//             if(e == 2)
//                 cout << e << " " << pt1z << " " << pt2z << endl;
            vsub(pt1, pt2, ptdir);
            ptw = - pt2z * 1.0 / (pt1z - pt2z);
            mul(ptdir, ptw - crit, ptdir);
            vadd(pt2, ptdir, pt2);


            mul(pt1, 1.0/pt1z, pt1);
            dotMV(K, pt1, pt1uv0);
            pt1uv[0] = pt1uv0[0]; pt1uv[1] = pt1uv0[1];
            roundVec(2, pt1uv, pt1uvr);

            pt2z = pt2[2];
//             if(e==2) cout << pt2z << endl;
            mul(pt2, 1.0/pt2z, pt2);
            dotMV(K, pt2, pt2uv0);
            pt2uv[0] = pt2uv0[0]; pt2uv[1] = pt2uv0[1];
            roundVec(2, pt2uv, pt2uvr);

            gp = max(abs(pt1uv[1] - pt2uv[1]), abs(pt1uv[0] - pt2uv[0]));
            // idi??idjï¼??¿åº¦ gpï¼?id0ï¼?gp 2ï¼?
            linspace(pt1uv[1], pt2uv[1], gp, idi0);
            roundVec(gp, idi0, idi);
            linspace(pt1uv[0], pt2uv[0], gp, idj0);
            roundVec(gp, idj0, idj);

            for(int i = 0; i < gp; i++)
            {
                id0[i][0] = idi[i];
                id0[i][1] = idj[i];
            }
            
            sizeid1 = chooseInt(1, ci, 1, cj, gp, id0, id1);
            
            for(int i = 0; i < sizeid1; i++)
            {
                ij[sizeij+i][0] = id1[i][0];
                ij[sizeij+i][1] = id1[i][1];
            }
            
//             if(e==2)
//                 for(int i = 0; i < sizeid1; i++)
//                     cout << ij[sizeij+i][0] << " " << ij[sizeij+i][1] << " ";
            
            sizeij = sizeij + sizeid1;
        }


        else if(pt1z <= 0 && pt2z > 0)
        {
            vsub(pt2, pt1, ptdir);
            ptw = - pt1z * 1.0 / (pt2z - pt1z);
            mul(ptdir, ptw - crit, ptdir);
            vadd(pt1, ptdir, pt1);

            pt1z = pt1[2];
            mul(pt1, 1.0/pt1z, pt1);
            dotMV(K, pt1, pt1uv0);
            pt1uv[0] = pt1uv0[0]; pt1uv[1] = pt1uv0[1];
            roundVec(2, pt1uv, pt1uvr);

            mul(pt2, 1.0/pt2z, pt2);
            dotMV(K, pt2, pt2uv0);
            pt2uv[0] = pt2uv0[0]; pt2uv[1] = pt2uv0[1];
            roundVec(2, pt2uv, pt2uvr);

            gp = max(abs(pt1uv[1] - pt2uv[1]), abs(pt1uv[0] - pt2uv[0]));
            // idi??idjï¼??¿åº¦ gpï¼?id0ï¼?gp 2ï¼?
            linspace(pt1uv[1], pt2uv[1], gp, idi0);
            roundVec(gp, idi0, idi);
            linspace(pt1uv[0], pt2uv[0], gp, idj0);
            roundVec(gp, idj0, idj);

            for(int i = 0; i < gp; i++)
            {
                id0[i][0] = idi[i];
                id0[i][1] = idj[i];
            }

            sizeid1 = chooseInt(1, ci, 1, cj, gp, id0, id1);
            for(int i = 0; i < sizeid1; i++)
            {
                ij[sizeij+i][0] = id1[i][0];
                ij[sizeij+i][1] = id1[i][1];
            }
            sizeij = sizeij + sizeid1;
        }
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
        if(ijres[k] <1)
            ijres[k] = 1;
        ijres[sizeij + k] = ij[k][1];
        if(ijres[sizeij + k] <1)
            ijres[sizeij + k] = 1;
    }

    
    
    delete[] idi0;
    delete[] idj0;
    delete[] idi;
    delete[] idj;
    delete[] id0;
    delete[] id1;
    delete[] ij;
}
