#pragma once
#include "LinearProgram.hpp"

class IntLinearProgram : public LinearProgram {
public:
    bool IsIntSolution();
    bool DualSimplex();
    void PrintSolution() override;
    void Solve() override;
};
