#pragma once
#include <vector>

using namespace std;

enum class SolutionState : short {
    NoFeasibleSol = 0,
    UnboundedSol = 1,
    HasSolution = 2
};

class LinearProgram {
public:
    void Input();
    void PrintCoeffMat(bool phase1 = true) const;
    bool IsUnitVector(int col) const;
    void Pivot(int pivot_row, int pivot_col, bool phase1 = true);
    bool Simplex(bool phase1 = true);
    SolutionState TwoPhaseMethod();
    virtual void Solve();
    virtual void PrintSolution();
private:
    int n, m;
    vector <vector <double>> A;
    vector <double> target;
    double eps = 1e-9;
};
