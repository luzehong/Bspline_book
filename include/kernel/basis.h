//
// Created by luzehong on 22-7-15.

#ifndef BASIS_H
#define BASIS_H

#include<vector>
#include "vector_operator.h"
using namespace std;


template <typename T> unsigned int FindSpan(unsigned int n , unsigned int p1, T u, const std::vector<T> &U )
{
    /* Determine the knot span index
     * Input: n, p, u, U
     * Return: the knot span index
     */
    if(u==U[n]) return(n-1);

    unsigned int low = p1, high = n+1;
    unsigned int mid = (low+high)/2;
    while (u< U[mid]||u>= U[mid+1])
    {
        if(u < U[mid])
            high = mid;
        else
            low  = mid;
        mid = (low + high) / 2;
    }
    return(mid);
}

template<typename T> void OneBasisFun(unsigned int p, unsigned int m, const vector<T> &U, int i, T u, T &Nip)
{
    /*Compute the basis function Nip
     * Input: p, m, U, i, u
     * Output: Nip;
     * */
    double N[p+1], saved, Uleft, Uright, temp;

    if((i == 0 && u == U[0]) || (i == m-p-1 && u == U[m])){ /*Special Case*/
        Nip = 1.0; return ;
    }

    if (u< U[i] || u>= U[i+p+1]) { /* Local property*/
        Nip = 0.0; return;
    }

    for (int j=0; j<=p; j++){ /* Initialize zeroth-degree functions*/
        if(u >= U[i+j] && u<U[i+j+1])
            N[j] = 1.0;
        else
            N[j] = 0.0;
    }

    for (int k=1; k<=p; k++){ /*Compute triangular table*/
        if( N[0] == 0.0)
            saved = 0.0;
        else
            saved = ((u-U[i])*N[0])/(U[i+k]-U[i]);
        for( int j=0; j< p-k+1; j++)
        {
            Uleft = U[i+j+1];
            Uright = U[i+j+k+1];
            if (N[j+1] == 0.0){
                N[j] = saved;
                saved = 0.0;
            }
            else{
                temp = N[j+1]/(Uright-Uleft);
                N[j] = saved + (Uright-u)*temp;
                saved = (u-Uleft)*temp;
            }
        }
    }
    Nip = N[0];


}

template <typename T> void BasisFuns(unsigned int i, double u, unsigned int p, const vector<T> &U, vector<T> &N)
{
    /* Compute the non-vanishing basis functions
     * Input: i, u, p,U;
     * Output: N
     */
    double left[p+1], right[p+1], saved, temp;

    N[0] = 1.0;
    for (int j=1; j<=p; j++){
        left[j] = u - U[i+1-j];
        right[j] = U[i+j] - u;
        saved = 0.0;
        for (int r=0; r<j; r++){
            temp = N[r]/(right[r+1] + left[j-r]);
            N[r] = saved + right[r+1] * temp;
            saved = left[j-r] * temp;
        }
        N[j] = saved;
    }

}

#endif
