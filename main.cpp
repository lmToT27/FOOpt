#include "IntLinearProgram.hpp"
#include <iostream>
using namespace std;

int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0); cout.tie(0);
    if (fopen("TASK.INP", "r")) {
        freopen("TASK.INP", "r", stdin);
        freopen("TASK.OUT", "w", stdout);
    }
    IntLinearProgram ilp;
    ilp.Input();
    ilp.Solve();
    return 0;
}
