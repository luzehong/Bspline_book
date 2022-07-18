//
// Created by luzehong on 22-7-15.
//

#ifndef VECTOR_OPERATOR_H
#define VECTOR_OPERATOR_H
#include <vector>
#include <cmath>
using namespace std;

vector<double> operator+(const vector<double> &v1, const vector<double> &v2){
    vector<double> v3(v1.size());
    for(int i=0;i<v1.size();i++){
        v3[i] = v1[i]+ v2[i];
    }
    return v3;
}

vector<double> operator-(const vector<double> &v1, const vector<double> &v2){
    vector<double> v3(v1.size());
    for(int i=0;i<v1.size();i++){
        v3[i] = v1[i]- v2[i];
    }
    return v3;
}

vector<double> operator*(const double &d, const vector<double> &v1){
    vector<double> v3(v1.size());
    for(int i=0;i<v1.size();i++){
        v3[i] = v1[i]*d;
    }
    return v3;
}
double norm(const vector<double> &V){
    /*This function is used to calculate the norm of a vector*/
    double result = 0;
    for(double i : V) //这是一种新的遍历数组的方法。
        result = result + i*i;
    result = sqrt(result);
    return result;
}


#endif //VECTOR_OPERATOR_H
