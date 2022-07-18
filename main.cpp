#include <iostream>
#include <fstream>
#include <vector>

#include "include/nurbs.h"
#include <Eigen/Dense>

int main() {
    std::cout << "Hello, World!" << std::endl;
    /*
    //We will construct a simple B-spline function and test the above three functions
    int p = 2;
    double u;
    vector<double> U = {0,0,0,1,2,3,4,5,6,6,6};
    //测试一下能不能输出离散点;

    vector<vector<double>> P = {{0,0},{1,1},{2,4},{3,9},{4,16},{5,25},{6,36},{7,49},{8,64}};
    vector<double> point(2);

    for(int i=0;i<P.size();i++) {
        for (int j = 0; j < P[0].size(); j++)
            cout << P[i][j] << " ";
        cout<<endl;
    }
    Bspline s1;
    s1.GlobalCurveInterp(P);
    s1.ShowContorlPoints();

    //Let‘s using the constructed b-spline curve, and output the curve points
    vector<double> uk(101);
    for(int i=0;i<=100;i++)
        uk[i] = double(i)/100;

    vector<vector<double>> Points(101,vector<double>(2));
    Points = s1.calculate(uk);


    ofstream ofs;
    ofs.open("./01.txt",ios::out);
    for (auto & Point : Points) {
        ofs<<Point[0]<<" "<<Point[1]<<endl;
    }
    ofs.close();
    */
    //Let's construct a B-spline surface!


    vector<vector<vector<double>>> Q(5,(vector<vector<double>>(5,vector<double>(3))));

    for(int i=0;i<5;i++) {
        for (int j = 0; j < 5; j++){
            Q[i][j][0] = i-2;
            Q[i][j][1] = j-2;
            Q[i][j][2] = (i-2)*(i-2)+(j-2)*(j-2);
        }
    }
    Bsplinesurface S1;
    S1.GlobalSurfaceInterp(Q);
  //  cout<<"showControlPointsMatrix:"<<endl;
  //  S1.showControlPointSMatrix();
  //  S1.showKnotVector();

    //Evaluate the points and test the function;
    vector<vector<vector<double>>> P(21,vector<vector<double>>(21,vector<double>(3)));
    vector<vector<vector<double>>> para(21,vector<vector<double>>(21,vector<double>(2)));

    for(int i=0;i<21;i++) {
        for (int j = 0; j < 21; j++) {
            para[i][j][0] = i * 0.05;
            para[i][j][1] = j * 0.05;
        }

    }


    cout<<"main Checked!"<<endl;
    P = S1.evalute(para);

    for(int i=0;i<P.size();i++) {
        for (int j = 0; j < P[0].size(); j++) {
            cout << "(" << P[i][j][0] << "," << P[i][j][1] << "," << P[i][j][2] << ") ";
        }
        cout<<endl;
    }

    // Output points;
    ofstream ofs;
    ofs.open("./01.txt",ios::out);
    for (auto & P1 : P) {
        for(auto & P2 : P1){
            for(auto & P3 : P2) {
                ofs << P3 << " ";
            }
            ofs<<endl;
        }

    }
    ofs.close();



}
