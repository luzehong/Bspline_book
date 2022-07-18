//
// Created by luzehong on 22-7-15.
//

#ifndef EVALUATE_H
#define EVALUATE_H

#include <vector>
#include "basis.h"
#include "vector_operator.h"

using namespace std;

template <typename T>
void CurvePoint(unsigned int n, unsigned int p, const vector<T> &U, const vector<vector<T>> &P,double u, vector<T> &C)
{
    /*compute curve point
     * Input: n,p,U,P,u
     * Output: C
     */
    int span = FindSpan(n,p,u,U);
    vector<double> N(p+1);
    BasisFuns(span,u,p,U,N);
    for(int i=0;i<P[0].size();i++) {
        C[i] = 0;
    }
    for(int i=0;i<=p;i++)
        C = C+N[i]*P[span-p+i]; //We  redefined the vector operator *, +.
}

#endif //EVALUATE_H
