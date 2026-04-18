#include "IntLinearProgram.hpp"
#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <cstdint>

bool IntLinearProgram::IsIntSolution() {
    for (int i = 0; i < n; i++) if (IsUnitVector(i)){
        double val = 0;
        for (int j = 0; j < m; j++) {
            if (A[j][i] == 1) {
                val = A[j][n];
                break;
            }
        }
        if (abs(val - round(val)) > eps) return false; 
    }
    return true;
}

void IntLinearProgram::PrintSolution() {
    cout << "Optimal value: " << A[m][n] << '\n';
    cout << "Optimal solution:\n";
    for (int i = 0; i < n; i++) {
        int64_t val = 0;
        if (IsUnitVector(i)) for (int j = 0; j < m; j++) {
            if (abs(A[j][i] - 1) < eps) {
                val = A[j][n];
                break;
            }
        }
        cout << "x" << i << " = " << val << '\n';
    }
}

bool IntLinearProgram::DualSimplex() {
    int COLWIDTH = 15;
    for (int j = 0; j < n; j++) {
        cerr << setw(COLWIDTH) << right << (j < n ? "x" + to_string(j) : "y" + to_string(j - n));
        cerr << " |";
    }
    cerr << setw(COLWIDTH) << "RHS";
    cerr << "\n" << string((n + 2) * COLWIDTH + 8, '-') << '\n';
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n; j++) {
            cerr << setw(COLWIDTH) << right << A[i][j];
            if (j < n) {
                cerr << " |";
            }
        }
        cerr << '\n';
    }
    cerr << "\n\n";
    while (true) {
        int pivot_row = -1;
        double min_rhs = -eps;
        for (int i = 0; i < m; i++) {
            if (A[i][n] < min_rhs) {
                min_rhs = A[i][n];
                pivot_row = i;
            }
        }
        if (pivot_row == -1) break;
        int pivot_col = -1;
        double best_ratio = 0;
        for (int j = 0; j < n; j++) {
            if (A[pivot_row][j] < -eps) {
                double ratio = abs(A[m][j] / A[pivot_row][j]);
                if (pivot_col == -1 || best_ratio > ratio) {
                    best_ratio = ratio;
                    pivot_col = j;
                }
            }
        }
        if (pivot_col == -1) return false;
        double val = A[pivot_row][pivot_col];
        for (int j = 0; j <= n; j++) {
            A[pivot_row][j] /= val;
            if (abs(A[pivot_row][j]) < eps) A[pivot_row][j] = 0;
        }
        for (int i = 0; i <= m; i++) if (i != pivot_row) {
            double f = A[i][pivot_col];
            if (abs(f) > eps) for (int j = 0; j <= n; j++) {
                A[i][j] -= f * A[pivot_row][j];
                if (abs(A[i][j]) < eps) {
                    A[i][j] = 0;
                }
            }
        }
        for (int j = 0; j < n; j++) {
            cerr << setw(COLWIDTH) << right << (j < n ? "x" + to_string(j) : "y" + to_string(j - n));
            cerr << " |";
        }
        cerr << setw(COLWIDTH) << "RHS";
        cerr << "\n" << string((n + 2) * COLWIDTH + 8, '-') << '\n';
        for (int i = 0; i <= m; i++) {
            for (int j = 0; j <= n; j++) {
                cerr << setw(COLWIDTH) << right << A[i][j];
                if (j < n) {
                    cerr << " |";
                }
            }
            cerr << '\n';
        }
        cerr << "\n\n";
    }
    return true;
}

void IntLinearProgram::Solve() {
    int COLWIDTH = 15;
    auto print_term = [&](double val, int var_idx) {
        if (abs(val) < eps) {
            cerr << string(COLWIDTH, ' ');
        } else {
            ostringstream oss;
            if (val > 0) {
                oss << "+";
            }
            oss << fixed << setprecision(2) << val << " * x" << var_idx;
            cerr << left << setw(COLWIDTH) << oss.str();
        }
    };

    SolutionState state = TwoPhaseMethod();
    if (state == SolutionState::NoFeasibleSol) return void(cout << "No feasible solution!!!");
    if (IsIntSolution()) return PrintSolution();
    int ite = 100;
    int real_n = n;
    cerr << "\nStarting Gomory Cuts...\n\n";
    while ((not IsIntSolution()) and (ite--)) {
        pair <int, int> chose = make_pair(-1, -1);
        double best_f = -1;
        for (int i = 0; i < n; i++) if (IsUnitVector(i)) {
            for (int j = 0; j < m; j++) if (A[j][i] == 1) {
                double f = A[j][n] - floor(A[j][n]);
                if (abs(f) > eps && abs(f) < 1 - eps) {
                    if (f > best_f) {
                        best_f = f;
                        chose = make_pair(j, i);
                    }
                }
            }
        }
        if (best_f == -1) break;
        for (int i = 0; i <= m; i++) {
            A[i].push_back(0);
            swap(A[i][n], A[i][n + 1]);
        }
        n++;
        A.push_back(vector <double> (n + 1, 0));
        swap(A[m], A[m + 1]);
        for (int i = 0; i <= n; i++) {
            if (i == n - 1) A[m][i] = 1;
            else if (i == n) A[m][i] = -best_f;
            else {
                double f = A[chose.first][i] - floor(A[chose.first][i]);
                if (f < eps || f > 1 - eps) f = 0;
                A[m][i] = -f;
            }
            if (abs(A[m][i]) < eps) A[m][i] = 0;
        }
        m++;
        cerr << "Gomory Cut: New constraint: ";
        for (int j = 0; j < n; j++) {
            print_term(A[m - 1][j], j);
        }
        cerr << " = " << fixed << setprecision(2) << A[m - 1][n + m] << "\n\n";
        if (!DualSimplex()) return void(cout << "Cant find integer solution!!!");
    }
    if (not IsIntSolution()) {
        cout << "Cant find integer solution!!!";
    } else {
        PrintSolution();
    }
}
