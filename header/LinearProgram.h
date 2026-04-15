#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>

using namespace std;

class LinearProgram {
public:
    void Input();
    void PrintCoeffMat(bool phase1 = true) const;
    bool IsUnitVector(int col) const;
    void Pivot(int pivot_row, int pivot_col, bool phase1 = true);
    void Simplex(bool phase1 = true);
    virtual void Solve();
    virtual void PrintSolution();
private:
    int n, m;
    vector <vector <double>> A;
    vector <double> target;
    double eps = 1e-9;
};