//
//  outil.h
//  MatlabC
//
//  Created by ?°æ??? on 2020/5/5.
//  Copyright Â© 2020 Andy. All rights reserved.
//

#ifndef outil_h
#define outil_h

#include <iostream>
#include <math.h>
using namespace std;

// ç´¢å?è½???
int idxM2C(const int i, const int j, const int m){
    return m * (j-1) + (i-1);
}
int idxC2M(const int i, const int j, const int m){
    return m * j + i;
}
int idxM2C3D(const int i, const int j, const int k, const int row, const int col){
    return row * col * (k-1) + row * (j-1) + (i-1);
}
int idxC2M3D(const int i, const int j, const int k, const int row, const int col){
    return row * col * (k) + row * (j) + (i);
}


// ??ä¹?
double* cross(const double edg1[3], const double edg2[3], double dst[3])
{
    dst[0] = edg1[1] * edg2[2] - edg1[2] * edg2[1];
    dst[1] = edg1[2] * edg2[0] - edg1[0] * edg2[2];
    dst[2] = edg1[0] * edg2[1] - edg1[1] * edg2[0];
    return dst;
}

// ?¹ä?
double dot(const double edg1[3], const double edg2[3])
{
    return edg1[0]*edg2[0] + edg1[1]*edg2[1] + edg1[2]*edg2[2];
}

// ?©é?? * ??????
double* dotMV(double mat[9], double vec[3], double dst[3])
{
    dst[0] = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2];
    dst[1] = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2];
    dst[2] = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2];
    return dst;
}

// ???? * ????
double* mul(double src[3], double l, double dst[3])
{
    dst[0] = src[0] * l;
    dst[1] = src[1] * l;
    dst[2] = src[2] * l;
    return dst;
}

// æ±?æ¨?
double norm(const double src[3])
{
    return sqrt(pow(src[0], 2) + pow(src[1], 2) + pow(src[2], 2));
}

// å½?ä¸?
double* normlz(double dst[3])
{
    double n = sqrt(dst[0]*dst[0] + dst[1]*dst[1] + dst[2]*dst[2]);
    dst[0] /= n;
    dst[1] /= n;
    dst[2] /= n;
    return dst;
}

// ??????
double* vadd(const double edg1[3], const double edg2[3], double dst[3])
{
    dst[0] = edg1[0] + edg2[0];
    dst[1] = edg1[1] + edg2[1];
    dst[2] = edg1[2] + edg2[2];
    return dst;
}

// ??????
double* vsub(const double edg1[3], const double edg2[3], double dst[3])
{
    dst[0] = edg1[0] - edg2[0];
    dst[1] = edg1[1] - edg2[1];
    dst[2] = edg1[2] - edg2[2];
    return dst;
}

// ?????·è?
double* vcp(double dst[3], const double src[3])
{
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    return dst;
}

// ????ç­?å·??°å?? linspace
double* linspace(double begin, double end, int num, double* dst)
{
    
    double gap = (end - begin) / (num - 1);
    double temp = begin;
    dst[0] = begin;
    dst[num-1] = end;
    for(int i = 1; i <= num-2; i++)
    {
        temp += gap;
        dst[i] = temp;
    }
    return dst;
}

// ç­????°ç?(n*2)ä¸?ç¬????¡ä»¶????ç´?
int chooseInt(int min1, int max1, int min2, int max2,
           int size, int src[][2], int dst[][2])
{
    int newsize = 0;
    for(int i = 0; i < size; i++)
    {
        if(src[i][0] >= min1 && src[i][0] <= max1 && src[i][1] >= min2 && src[i][1] <= max2)
        {
            dst[newsize][0] = src[i][0];
            dst[newsize][1] = src[i][1];
            newsize++;
        }
    }
    return newsize;
}

// ?°ç???ç´?????
void roundVec(int size, double* src, int* dst)
{
    for(int i = 0; i < size; i++)
    {
        dst[i] = round(src[i]);
    }
}


#endif /* outil_h */
