//
// Created by luzehong on 22-7-17.
//

#ifndef BSPLINESURFACE_H
#define BSPLINESURFACE_H
#include<vector>
#include<iostream>
#include "../kernel/parameter.h"
#include "bsplinexx.h"

using namespace std;
class Bsplinesurface{
private:
    unsigned int degree_u; /*B-spline polynomial degree*/
    unsigned int degree_v;
    unsigned int num_u; /*Nuber of control points*/
    unsigned int num_v; /*Nuber of control points*/
    unsigned int knotnumber_u;
    unsigned int knotnumber_v;
    unsigned int dimension; /*The Bspline's dimension: 2-d or 3-d;*/
    vector<vector<vector<double>>> ControlPointsMatrix;
    vector<double> KnotVector_u;
    vector<double> KnotVector_v;

public:
    void GlobalSurfaceInterp(vector<vector<vector<double>>> &InterPoints, unsigned int pu = 3, unsigned int pv = 3);
    void GlobalSurfaceInterp(vector<vector<vector<double>>> &InterPoints, vector<vector<double>> &Parameters,
                             unsigned  int pu=3, unsigned int pv =3);
    vector<vector<vector<double>>> evalute(vector<vector<vector<double>>> &u_v);

    vector<double> KnotVector_initial(unsigned int nu, unsigned int p=3);
    static vector<double> KnotVector_initial(vector<double> &uk, unsigned int p = 3);
    static vector<vector<vector<double>>> ControlPointsMatrix_initial(unsigned  int nu, unsigned int nv, unsigned int di=3);

    void showKnotVector();
    void showControlPointSMatrix();

};

void Bsplinesurface::GlobalSurfaceInterp(vector<vector<vector<double>>> &InterPoints, unsigned int pu, unsigned int pv)
{
    degree_u = pu;
    degree_v = pv;
    num_u    = InterPoints.size();
    num_v    = InterPoints[0].size();
    knotnumber_u = num_u + degree_u +1;
    knotnumber_v = num_v + degree_v +1;
    dimension= InterPoints[0][0].size();
    ControlPointsMatrix = ControlPointsMatrix_initial(num_u,num_v,dimension);
    vector<double> uk(num_u),vk(num_v);


    SurfaceMeshParams(InterPoints,uk,vk);
    KnotVector_u = KnotVector_initial(uk);
    KnotVector_v = KnotVector_initial(vk);


    /* Then we will calculate controlpoints in each direction*/
    vector<vector<double>> Q(num_u,vector<double>(dimension));
    vector<vector<double>> Rrow(num_v,vector<double>(dimension));
    vector<vector<vector<double>>> R(num_u,vector<vector<double>>(num_v,vector<double>(dimension)));
    vector<vector<double>> Controlu(num_u,vector<double>(dimension));
    vector<vector<double>> Controlv(num_v,vector<double>(dimension));


    Bspline UdirectionSpline;
    Bspline VdirectionSpline;


    for(int l=0; l<num_v;l++){
        for(int i=0;i<num_u;i++)
            Q[i] = InterPoints[i][l];
        UdirectionSpline.GlobalCurveInterp(Q);
        Controlu = UdirectionSpline.getControlPoints();
        for(int i=0;i<num_u;i++) R[i][l] = Controlu[i];
    }

    for(int i=0;i<num_u;i++){
        for(int j=0;j<num_v;j++) Rrow[j] = R[i][j];
        VdirectionSpline.GlobalCurveInterp(Rrow);
        Controlv = VdirectionSpline.getControlPoints();
        for(int j=0;j<num_v;j++)  ControlPointsMatrix[i][j] = Controlv[j];
    }

}

void Bsplinesurface::GlobalSurfaceInterp(vector<vector<vector<double>>> &InterPoints, vector<vector<double>> &Parameters,
                                         unsigned int pu, unsigned int pv){

}

vector<vector<vector<double>>>
Bsplinesurface::ControlPointsMatrix_initial(unsigned int nu, unsigned int nv, unsigned int di) {
    vector<vector<vector<double>>> M(nu,vector<vector<double>>(nv,vector<double>(di)));
    return M;
}

void Bsplinesurface::showControlPointSMatrix() {
    for(auto & i : ControlPointsMatrix) {
        for (int j = 0; j < ControlPointsMatrix[0].size(); j++) {
            cout<<"(";
            for (int k = 0; k < ControlPointsMatrix[0][0].size(); k++)
                cout << i[j][k] << " ";
            cout<<")";
        }
        cout<<endl;
    }

}

vector<double> Bsplinesurface::KnotVector_initial(unsigned int nu, unsigned int p) {
    /* If not define the KnotVector, it will be defined automatically in [0,1] */
    unsigned int n = nu-1; /*这里为了与The NURBS Book的书一致，那本书上n+1才是控制点的个数，m+1才是节点适量的长度，其中m=n+p+1*/
    unsigned int m = p + n +1;
    vector<double> Knot(m+1);
    for (int i=0;i<=p;i++){
        Knot[i] = 0;
    }
    for (unsigned int i= m-p;i<=m;i++){
        Knot[i] = 1;
    }
    for (unsigned int i=1;i<=n-p;i++){
        Knot[i+p] = double(i)/(n-p+1);
    }
    return Knot;
}

vector<double> Bsplinesurface::KnotVector_initial(vector<double> &uk, unsigned int p) {
    unsigned int n = uk.size()-1;
    unsigned int m = n + p +1;
    vector<double> Knot(m+1);
    double sum;
    for (int i=0;i<=p;i++){
        Knot[i] = 0;
    }
    for (unsigned int i= m-p;i<=m;i++){
        Knot[i] = 1;
    }
    for (unsigned int j=1;j<=n-p;j++){
        sum = 0;
        for (unsigned int i=j ;i<=j+p-1;i++){
            sum = sum + uk[i];
        }
        Knot[j+p] = sum/p;
    }
    return Knot;
}

vector<vector<vector<double>>> Bsplinesurface::evalute(vector<vector<vector<double>>> &u_v) {
   /*Evalute the points in B-spline surface;
    * Input: parameter matrxi u_v: vector<vector<double>>
    * Output:Points matri: vector<vector<vector<double>>>*/
   vector<vector<vector<double>>> Points(u_v.size(),vector<vector<double>>(u_v[0].size(),vector<double>(dimension)));
   unsigned int uspan,vspan,uind,vind;
   unsigned int  p = degree_u,q=degree_v;
   vector<vector<double>> temp(p+1,vector<double>(dimension));
   vector<double> S(dimension);

   vector<double> Nu(p+1),Nv(q+1);
   for(int i=0;i<u_v.size();i++){
       for(int j=0;j<u_v[0].size();j++){

           uspan = FindSpan(num_u,p,u_v[i][j][0],KnotVector_u);
           BasisFuns(uspan,u_v[i][j][0],p,KnotVector_u,Nu);
           vspan = FindSpan(num_v,q,u_v[i][j][1],KnotVector_v);
           BasisFuns(vspan,u_v[i][j][1],q,KnotVector_v,Nv);
           uind = uspan - p;

           for(int k=0;k<dimension;k++){ S[k] = 0.0;}

           for(int l=0;l<=q;l++)
           {
               for(int k=0;k<dimension;k++){temp[l][k] = 0.0;} //temp = 0.0;

               vind = vspan - q;
               for(int k=0;k<=p;k++)
               {
                   temp[l] = temp[l] + Nu[k] * ControlPointsMatrix[uind + k][vind+l];
               }
               S = S + Nv[l]*temp[l];
           }
           Points[i][j]= S;
       }
   }
    return Points;
}

void Bsplinesurface::showKnotVector() {
    cout<<"KnotVector_u:"<<endl;
    for(double i : KnotVector_u)
        cout<<i<<endl;
    cout<<"KnotVector_v:"<<endl;
    for(double i : KnotVector_v)
        cout<<i<<endl;

}


#endif //BSPLINESURFACE_H
