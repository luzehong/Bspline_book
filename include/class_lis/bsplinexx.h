//
// Created by luzehong on 22-7-15.
//
//我们希望将B样条的所有内容都放到B样条模板中去，然后模板中对应的函数与The NURBS Book保持一致;
#ifndef BSPLINEXX_H
#define BSPLINEXX_H
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "../kernel/parameter.h"

using namespace std;
class Bspline{
private:
    unsigned int degree; /*B-spline polynomial degree*/
    unsigned int num; /*Nuber of control points*/
    unsigned int knotnumber;
    unsigned int dimension; /*The Bspline's dimension: 2-d or 3-d;*/
    vector<vector<double>> ControlPoints;
    vector<double> KnotVector;
public:
    static vector<double> KnotVector_initial( unsigned int nu, unsigned int p = 3);
    static vector<double> KnotVector_initial(vector<double> &uk, unsigned int p = 3);
    static vector<vector<double>> ControlPoints_initial(unsigned int nu,unsigned int di);
    void initial( int p, int nu,  int d, const vector<vector<double>> &point, const vector<double> &kn);
    void initial(int p, int nu, int d, const vector<vector<double>> &point);
    void GlobalCurveInterp(vector<vector<double>> &Q, unsigned int p = 3);
    void GlobalCurveInterp(vector<vector<double>> &Q,vector<double> &uk, unsigned int p = 3);
    vector<vector<double>> evaluate(vector<double> &uk); /* Calculate the points corresponding to the parameter vectors uk */

    unsigned int getnum() const{return num;};
    unsigned int getdegree() const{return degree;};
    unsigned int getknotnumber() const{return knotnumber;};
    unsigned int getdimension() const{return dimension;};
    vector<vector<double>> getControlPoints(){return ControlPoints;};
    vector<double> getKnotVector(){return KnotVector;};
    double getControlPointsbyindex( unsigned int i, unsigned int j){return ControlPoints[i][j];};
    double getKnotVectorbyindex(unsigned int i){return KnotVector[i];};
    void ShowContorlPoints();

};

vector<double> Bspline::KnotVector_initial(unsigned int nu, unsigned int p) {

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

vector<double> Bspline::KnotVector_initial(vector<double> &uk, unsigned int p) {
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

vector<vector<double>> Bspline::ControlPoints_initial(unsigned int nu, unsigned int di) {
    vector<vector<double>> p(nu,vector<double>(di));
    return p;
}
void Bspline::initial( int p, int nu,  int d, const vector<vector<double>> &point, const vector<double> &kn) {
    degree = p;
    num    = nu;
    knotnumber = nu+p+1;
    dimension = d;
    ControlPoints = point;
    KnotVector = kn;
}
void Bspline::initial(int p, int nu, int d, const vector<vector<double>> &point) {
    degree = p;
    num    = nu;
    knotnumber = nu+p+1; /* m+1 = n+1+p+1 as m = n+p+1 in The NURBS Book;*/
    dimension  = d;
    ControlPoints = point;
    KnotVector = KnotVector_initial(nu,p);

}



void Bspline::GlobalCurveInterp(vector<vector<double>> &Q, unsigned int p){

    /*Initial the B-spline parameter by the interpolated points Q and the polynomial degree p;*/
    degree = p;
    num = Q.size();
    knotnumber = num + p +1;
    dimension = Q[0].size();
    ControlPoints = ControlPoints_initial(num,dimension);
    vector<double> uk(Q.size());

    uk = equally_space(Q);
    KnotVector = KnotVector_initial(uk,p);

    vector<double> N(p+1);
    Eigen::MatrixXd A(num,num);
    Eigen::VectorX<double> x(num),b(num);

    A.setZero();
    x.setZero();
    b.setZero();

    /*Calculate the i-th row of the linear equation matrix*/
    unsigned int n = num-1, m = n+p+1, span;
    for(int i=0; i<=n; i++){
        span = FindSpan(num,p,uk[i],KnotVector);
        BasisFuns(span,uk[i],p,KnotVector,N);
        for(int j=0; j<p+1;j++)
            A(i, span - p + j) = N[j];
    }
    /*Solve the linear equation*/
    Eigen::MatrixXd InA = A.inverse();
    for (int j=0; j<dimension; j++) {
        for (int i = 0; i < num; i++)
        {
            b(i) = Q[i][j];
        }

        x = InA*b;
        for (int k = 0; k < num; k++)
        {
            ControlPoints[k][j] = x(k);
        }
    }

}

void Bspline::GlobalCurveInterp(vector<vector<double>> &Q,vector<double> &uk, unsigned int p){
    degree = p;
    num = Q.size();
    knotnumber = num + p +1;
    dimension = Q[0].size();
    ControlPoints = ControlPoints_initial(num,dimension);

    KnotVector = KnotVector_initial(uk,p);
    vector<double> N(p+1);
    Eigen::MatrixXd A(num,num);
    Eigen::VectorX<double> x(num),b(num);

    A.setZero();
    x.setZero();
    b.setZero();

    /*Calculate the i-th row of the linear equation matrix*/
    unsigned int n = num-1, m = n+p+1, span;
    for(int i=0; i<=n; i++){
        span = FindSpan(num,p,uk[i],KnotVector);
        BasisFuns(span,uk[i],p,KnotVector,N);
        for(int j=0; j<p+1;j++) A(i, span - p + j) = N[j];
    }


    /*Solve the linear equation*/
    Eigen::MatrixXd InA = A.inverse();
    for (int j=0; j<dimension; j++) {
        for (int i = 0; i < num; i++)
        {
            b(i) = Q[i][j];
        }

        x = InA*b;
        for (int k = 0; k < num; k++)
        {
            ControlPoints[k][j] = x(k);
        }
    }

}



void Bspline::ShowContorlPoints() {
    cout<<"ControlPoints:"<<endl;
    for(auto & ControlPoint : ControlPoints)
    {
        for(auto &ControlPoin : ControlPoint)
            cout<<ControlPoin<<" ";
        cout<<endl;}
}

vector<vector<double>> Bspline::evaluate(vector<double> &uk) {
    vector<vector<double>> Points(uk.size(),vector<double>(dimension));
    vector<double> C(dimension);
    for(int i=0;i<uk.size();i++){
        CurvePoint(num,degree,KnotVector,ControlPoints,uk[i],C);
        for(int j=0;j<dimension;j++) Points[i][j]= C[j];
    }
    return Points;
}

#endif //BSPLINEXX_H
