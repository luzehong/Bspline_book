//
// Created by luzehong on 22-7-15.
//

#ifndef PARAMETER_H
#define PARAMETER_H
#include <vector>
#include <iostream>
#include "vector_operator.h"
using namespace std;

vector<double> equally_space(vector<vector<double>> &P)
{
    unsigned int n= P.size()-1;
    vector<double> uk(P.size());
    for (int i=0;i<P.size();i++){
       uk[i] = double(i)/n;
    }
    return uk;
}

vector<double> centripetal();
vector<double> chord_length(vector<vector<double>> &P){
}


void SurfaceMeshParams(vector<vector<vector<double>>> &Q,
                       vector<double> &uk, vector<double> &vl)
{/* This is a chord_length parameter method for rectangular distribute points
 * */

    unsigned num = Q.size(), nvm = Q[0].size();
    unsigned n = num-1, m = nvm-1;
    uk[0] = 0;
    uk[n] = 1;
    double total,d;
    vector<double> cds(n+1);
    vector<double> cdv(m+1);

    for(int k=1;k<n;k++) {uk[k] = 0;}

    for(int l=0;l<=m;l++)
    {
        total = 0;
        for(int k=1;k<=n;k++){
            cds[k] = norm(Q[k][l]-Q[k-1][l]);
            total  = total + cds[k];
        }
        if(total == 0.0)
            num = num-1;
        else{
            d = 0.0;
            for(int k=1;k<n;k++)
            {
                d = d + cds[k];
                uk[k] = uk[k] + d/total;
            }
        }
    }
    if(num == 0) {cout<<"error in SurfaceMeshParams() function"<<endl; return;}
    for(int k=1;k<n;k++)  uk[k] = uk[k]/num;
    /*Then the parameters in v-direction can be calculated by the same way!*/


    vl[0] = 0;
    vl[m] = 1;
    for(int l=1;l<m;l++)  {vl[l] = 0;}

    for(int k=0;k<=n;k++){
        total = 0;
        for(int l=1;l<=m;l++){
            cdv[l] = norm(Q[k][l]-Q[k][l-1]);
            total = total + cdv[l];
        }
        if(total == 0.0)  nvm = nvm -1;
        else{
            d = 0.0;
            for(int l=1;l<m;l++)
            {
                d = d+cdv[l];
                vl[l] = vl[l] + d/total;
            }
        }
    }
    if(nvm == 0) {cout<<"error in SurfaceMeshparams() second part"<<endl;return;}
    for(int k=1;k<m;k++) vl[k] = vl[k]/nvm;

}
#endif //PARAMETER_H
