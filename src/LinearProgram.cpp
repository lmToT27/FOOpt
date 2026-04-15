#include "LinearProgram.h"

using namespace std;

void LinearProgram::Input() {
    try {
        cin >> n >> m;
        if (cin.fail() || n <= 0 || m <= 0) {
            throw invalid_argument("Invalid input. Please enter positive integers.");
        }

        target.resize(n);
        for (double &x : target) {
            cin >> x;
            if (cin.fail()) {
                throw invalid_argument("Invalid input. Please enter a real number.");
            }
        }

        A.resize(m + 1, vector <double> (n + m + 1));
        for (int i = 0; i < m; i++) {
            A[i][n + i] = 1.0;
            for (int j = 0; j <= n; j++) {
                double x; cin >> x;
                if (cin.fail()) {
                    throw invalid_argument("Invalid input. Please enter a real number.");
                }
                try {
                    if (j != n) {
                        A[i][j] = x;
                        A[m][j] -= x;
                    } else {
                        if (x < 0) {
                            throw invalid_argument("RHS must be a non-negative real number.");
                        }
                        A[i][j + m] = x;
                        A[m][j + m] -= x;
                    }
                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << '\n';
                    exit(0);
                }
            }
        }
        cerr << '\n';
    } catch (const invalid_argument& e) {
        cerr << "Error: " << e.what() << endl;
    }

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

    cerr << "Maximize:\n";
    for (int i = 0; i < n; i++) {
        print_term(target[i], i);
    }

    cerr << "\nConstraints:\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            print_term(A[i][j], j);
        }
        // Giữ nguyên vế phải của phương trình
        cerr << " = " << fixed << setprecision(2) << A[i][n + m] << '\n';
    }
    cerr << "\n\n";
    exit(0);
}
void LinearProgram::PrintCoeffMat(bool phase1) const {
    int COLWIDTH = 15;
    for (int j = 0; j < n + (phase1 ? m : 0); j++) {
        cerr << setw(COLWIDTH) << right << (j < n ? "x" + to_string(j) : "y" + to_string(j - n));
        if (j < n + (phase1 ? m : 0)) {
            cerr << " |";
        }
    }
    cerr << setw(COLWIDTH) << "RHS";
    cerr << "\n" << string((n + (phase1 ? m : 0) + 2) * COLWIDTH + 8, '-') << '\n';
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n + (phase1 ? m : 0); j++) {
            cerr << setw(COLWIDTH) << right << A[i][j];
            if (j < n + (phase1 ? m : 0)) {
                cerr << " |";
            }
        }
        cerr << '\n';
    }
}
void LinearProgram::Pivot(int pivot_row, int pivot_col, bool phase1) {
    double f = A[pivot_row][pivot_col];
    for (int j = 0; j <= n + (phase1 ? m : 0); j++) {
        A[pivot_row][j] /= f;
        if (abs(A[pivot_row][j]) < eps) A[pivot_row][j] = 0;
    }
    for (int i = 0; i <= m; i++) {
        if (i != pivot_row) {
            double g = A[i][pivot_col];
            for (int j = 0; j <= n + (phase1 ? m : 0); j++) {
                A[i][j] -= g * A[pivot_row][j];
                if (abs(A[i][j]) < eps) A[i][j] = 0;
            }
        }
    }
}
bool LinearProgram::IsUnitVector(int col) const {
    bool ok = false;
    for (int i = 0; i <= m; i++) {
        if (A[i][col] == 1) {
            if (ok) return false;
            ok = true;
        } else if (A[i][col] != 0) {
            return false;
        }
    }
    return ok;
}
void LinearProgram::Simplex(bool phase1) {
    while (true) {
        int pivot_col = -1;
        for (int i = 0; i < n + (phase1 ? m : 0); i++) if (A[m][i] < 0) {
            if (pivot_col == -1 || A[m][i] < A[m][pivot_col]) {
                pivot_col = i;
            }
        }
        if (pivot_col == -1) break;
        int pivot_row = -1;
        double min_ratio = 1e18;
        for (int i = 0; i < m; i++) if (A[i][pivot_col] > 0) {
            double ratio = A[i][n + (phase1 ? m : 0)] / A[i][pivot_col];
            if (ratio < min_ratio) {
                min_ratio = ratio;
                pivot_row = i;
            }
        }
        if (pivot_row == -1) {
            cerr << "Unbounded solution!!!" << endl;
            exit(0);
        }
        Pivot(pivot_row, pivot_col, phase1);
        PrintCoeffMat(phase1);
        cerr << "\n\n";
    }
}
void LinearProgram::Solve() {
    cerr << "PHASE 1:\n\n";
    PrintCoeffMat();
    cerr << "\n\n";
    Simplex();
    if (A[m][n + m] < 0) {
        cerr << "No feasible solution" << '\n';
        exit(0);
    }
    PrintCoeffMat();
    while (true) {
        int pivot_col = -1;
        for (int i = n; i < n + m; i++) if (IsUnitVector(i)) {
            pivot_col = i;
            break;
        }
        if (pivot_col == -1) break;
        int pivot_row = -1;
        for (int i = 0; i < m; i++) if (A[i][pivot_col] == 1) {
            pivot_row = i;
            break;
        }
        // assert(pivot_row != -1);
        int tmp = -1;
        for (int i = 0; i < n; i++) {
            if (A[pivot_row][i] != 0) {
                tmp = i;
                break;
            }
        }
        if (tmp != -1) {
            Pivot(pivot_row, tmp);
        } else {
            // delete row pivot_row and column pivot_col
            vector <vector <double>> newA(m, vector <double> (n + m));
            for (int i = 0, r = 0; i <= m; i++) {
                if (i == pivot_row) continue;
                for (int j = 0, c = 0; j <= n + m; j++) {
                    if (j == pivot_col) continue;
                    newA[r][c] = A[i][j];
                    c++;
                }
                r++;
            }
            A = newA;
            n--;
            m--;
        }
    }
    cerr << "Done Phase 1!!!\n\nPHASE 2:\n\n";
    vector <int> basis;
    for (int i = 0; i < n; i++) {
        if (IsUnitVector(i)) {
            basis.push_back(i);
        }
    }
    // delete y variables:
    vector <vector <double>> newA(m + 1, vector <double> (n + 1));
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j < n; j++) newA[i][j] = A[i][j];
        newA[i][n] = A[i][n + m];
    }
    A = newA;
    for (int i = 0; i < n; i++) {
        A[m][i] = -target[i];
        if (abs(A[m][i]) < eps) A[m][i] = 0;
    }
    PrintCoeffMat(false);
    cerr << "\n\n";
    for (int &j : basis) {
        for (int i = 0; i < m; i++) {
            if (A[i][j] == 1) {
                for (int k = 0; k <= n; k++) {
                    A[m][k] += A[i][k] * target[j];
                }
                break;
            }
        }
    }
    PrintCoeffMat(false);
    cerr << "\n\n";
    Simplex(false);
}
void LinearProgram::PrintSolution() {
    cerr << "Optimal value: " << A[m][n] << '\n';
    cerr << "Optimal solution:\n";
    for (int i = 0; i < n; i++) {
        double val = 0;
        for (int j = 0; j < m; j++) {
            if (A[j][i] == 1 && IsUnitVector(i)) {
                val = A[j][n];
                break;
            }
        }
        cerr << "x" << i << " = " << val << '\n';
    }
}
/*
    Input co dang:
    n m
    c1 c2 c3 ... cn
    a11 a12 a13 ... a1n b1
    a21 a22 a23 ... a2n b2
    ...
    am1 am2 am3 ... amn bm 
    Trong do:
    - n: so bien
    - m: so rang buoc
    - aij: he so cua bien xj trong rang buoc i
    - bi: phan tu ben phai cua rang buoc i
    - ci: he so cua bien xi trong ham muc tieu
    - Dieu kien: can tim max ham muc tieu, cac rang buoc co dang a1x1 + a2x2 + ... + anxn = bi >= 0, xi >= 0
*/