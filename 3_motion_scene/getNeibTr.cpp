//
//  getNeibTr.cpp
//  MatlabC
//
//  Created by ?°æ??? on 2020/5/5.
//  Copyright Â© 2020 Andy. All rights reserved.
//
#pragma clang diagnostic push
#pragma clang diagnostic ignored"-Wshorten-64-to-32"

#include <iostream>
#include <vector>
#include "mex.h"
#include "outil.h"
using namespace std;




void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // *prhs = {i, j, Map(ci cj), Cover(ci cj nTr), List(1 1 nTr), ci, cj, gap, nTr}
    // *plhs = {List(1 1 nTr)}
    
    double *ip, *jp, *Mapp, *Coverp, *Listp, *cip, *cjp, *gapp, *nTrp;
    ip = mxGetPr(prhs[0]);
    jp = mxGetPr(prhs[1]);
    Mapp = mxGetPr(prhs[2]);
    Coverp = mxGetPr(prhs[3]);
    Listp = mxGetPr(prhs[4]);
    cip = mxGetPr(prhs[5]);
    cjp = mxGetPr(prhs[6]);
    gapp = mxGetPr(prhs[7]);
    nTrp = mxGetPr(prhs[8]);
    
    int Mapm = mxGetM(prhs[2]);
//    int Mapn = mxGetN(prhs[2]);
    int i = *ip, j = *jp, ci = *cip, cj = *cjp, gap = *gapp, nTr = *nTrp;

    
    // get element in Cover   [idxC2M3D(i, j, k, r, c)]
    int *Coverdim;
    Coverdim = (int*)mxGetDimensions(prhs[4]);
    int Coverr = Coverdim[0], Coverc = Coverdim[2];
    
    // add the triangle of last instant
    int curTr = Mapp[idxM2C(i, j, Mapm)];
    if(curTr != 0) Listp[curTr - 1] = curTr;
    
    // add the triangle on same surface into List
    switch (curTr) {
        // Cb
        case 11:
            Listp[10] = 11; Listp[11] = 12; Listp[14] = 15; Listp[20] = 21;
            break;
        case 12:
            Listp[10] = 11; Listp[11] = 12; Listp[16] = 17; Listp[18] = 19;
            break;
        case 13:
            Listp[12] = 13; Listp[13] = 14; Listp[21] = 22; Listp[15] = 16;
            break;
        case 14:
            Listp[12] = 13; Listp[13] = 14; Listp[17] = 18; Listp[19] = 20;
            break;
        case 15:
            Listp[14] = 15; Listp[15] = 16; Listp[10] = 11; Listp[20] = 21;
            break;
        case 16:
            Listp[14] = 15; Listp[15] = 16; Listp[12] = 13; Listp[18] = 19;
            break;
        case 17:
            Listp[16] = 17; Listp[17] = 18; Listp[11] = 12; Listp[21] = 22;
            break;
        case 18:
            Listp[16] = 17; Listp[17] = 18; Listp[13] = 14; Listp[19] = 20;
            break;
        case 19:
            Listp[18] = 19; Listp[19] = 20; Listp[11] = 12; Listp[15] = 16;
            break;
        case 20:
            Listp[18] = 19; Listp[19] = 20; Listp[13] = 14; Listp[17] = 18;
            break;
        case 21:
            Listp[20] = 21; Listp[21] = 22; Listp[10] = 11; Listp[14] = 15;
            break;
        case 22:
            Listp[20] = 21; Listp[21] = 22; Listp[12] = 13; Listp[16] = 17;
            break;
        // Te
        case 23:
            Listp[22] = 23; Listp[23] = 24; Listp[24] = 25; Listp[25] = 26;
            break;
        case 24:
            Listp[22] = 23; Listp[23] = 24; Listp[24] = 25; Listp[25] = 26;
            break;
        case 25:
            Listp[22] = 23; Listp[23] = 24; Listp[24] = 25; Listp[25] = 26;
            break;
        case 26:
            Listp[22] = 23; Listp[23] = 24; Listp[24] = 25; Listp[25] = 26;
            break;
        // Ct
        case 27:
            Listp[26] = 27; Listp[27] = 28; Listp[28] = 29; Listp[32] = 33;
            break;
        case 28:
            Listp[26] = 27; Listp[27] = 28; Listp[31] = 32; Listp[34] = 35;
            break;
        case 29:
            Listp[28] = 29; Listp[29] = 30; Listp[26] = 27; Listp[32] = 33;
            break;
        case 30:
            Listp[28] = 29; Listp[29] = 30; Listp[30] = 31; Listp[35] = 36;
            break;
        case 31:
            Listp[30] = 31; Listp[31] = 32; Listp[32] = 33; Listp[29] = 30;
            break;
        case 32:
            Listp[30] = 31; Listp[31] = 32; Listp[27] = 28; Listp[33] = 34;
            break;
        case 33:
            Listp[32] = 33; Listp[26] = 27; Listp[28] = 29; Listp[30] = 31;
            break;
        case 34:
            Listp[33] = 34; Listp[34] = 35; Listp[35] = 36; Listp[31] = 32;
            break;
        case 35:
            Listp[33] = 34; Listp[34] = 35; Listp[35] = 36; Listp[27] = 28;
            break;
        case 36:
            Listp[33] = 34; Listp[34] = 35; Listp[35] = 36; Listp[29] = 30;
            break;
        default:
            break;
    }
    
    // get neighbor triangles (neighbors and covered by neighbors)
    if(i+gap <= ci)
    {
        if(j + gap <= cj)
        {
            curTr = Mapp[idxM2C(i+gap/2, j+gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;

            curTr = Mapp[idxM2C(i+gap, j+gap/2, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;

            curTr = Mapp[idxM2C(i+gap, j+gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;

            curTr = Mapp[idxM2C(i, j+gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;

            curTr = Mapp[idxM2C(i+gap, j, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            
//            for(int c = 1; c <= nTr; c++)
//            {
//                if(Coverp[idxM2C3D(i+gap/2, j+gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i+gap, j+gap/2, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i+gap, j+gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i, j+gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i+gap, j, c, Coverr, Coverc)] != 0 )
//                    Listp[c-1] = c;
//            }
            
        }
        if(j - gap >= 1)
        {
            curTr = Mapp[idxM2C(i+gap/2, j-gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i+gap, j-gap/2, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i+gap, j-gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i, j-gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i+gap, j, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            
//            for(int c = 1; c <= nTr; c++)
//            {
//                if(Coverp[idxM2C3D(i+gap/2, j-gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i+gap, j-gap/2, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i+gap, j-gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i, j-gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i+gap, j, c, Coverr, Coverc)] != 0 )
//                    Listp[c-1] = c;
//            }
        }
    }
    if(i-gap >= 1)
    {
        if(j + gap <= cj)
        {
            curTr = Mapp[idxM2C(i-gap/2, j+gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i-gap, j+gap/2, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i-gap, j+gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i, j+gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i-gap, j, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            
//            for(int c = 1; c <= nTr; c++)
//            {
//                if(Coverp[idxM2C3D(i-gap/2, j+gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i-gap, j+gap/2, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i-gap, j+gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i, j+gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i-gap, j, c, Coverr, Coverc)] != 0 )
//                    Listp[c-1] = c;
//            }
            
        }
        if(j - gap >= 1)
        {
            curTr = Mapp[idxM2C(i-gap/2, j-gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i-gap, j-gap/2, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i-gap, j-gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i, j-gap, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            curTr = Mapp[idxM2C(i-gap, j, Mapm)];
            if(curTr != 0) Listp[curTr - 1] = curTr;
            
//            for(int c = 1; c <= nTr; c++)
//            {
//                if(Coverp[idxM2C3D(i-gap/2, j-gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i-gap, j-gap/2, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i-gap, j-gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i, j-gap, c, Coverr, Coverc)] != 0 ||
//                   Coverp[idxM2C3D(i-gap, j, c, Coverr, Coverc)] != 0 )
//                    Listp[c-1] = c;
//            }
        }
    }
    
    
    
}

#pragma clang diagnostic pop
