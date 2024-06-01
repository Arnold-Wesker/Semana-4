// Primero genere las matrices A y B con los datos de entrada. Use el algoritmo de Strassen. Este divide las matrices en cuadrantes y 
// realiza 7 multiplicaciones recursivas. Para implementarlo use funciones para sumar y restar matrices. Estas me permitieron calcular 
// las 7 multiplicaciones eficientemente. Además, use funciones para dividir a la matriz en cuadrantes y unir los cuadrantes en una 
// matriz. Estas son necesarias para la recursividad del algoritmo. La función de división me permitió dividir las matrices en 
// submatrices más pequeñas, mientras que la función de unión me permitió combinar las submatrices resultantes en la matriz final. 
// Finalmente, calculé el vector V a partir de la matriz C.
#include <vector>
#include <iostream>

using namespace std;

vector<vector<uint64_t>> add(vector<vector<uint64_t>> A, vector<vector<uint64_t>> B) {
    int n = A.size();
    vector<vector<uint64_t>> C(n, vector<uint64_t>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

vector<vector<uint64_t>> sub(vector<vector<uint64_t>> A, vector<vector<uint64_t>> B) {
    int n = A.size();
    vector<vector<uint64_t>> C(n, vector<uint64_t>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

void split(vector<vector<uint64_t>>& cd, vector<vector<uint64_t>>& gc, size_t iB, size_t jB) {
    for(size_t i1 = 0, i2 = iB; i1 < gc.size(); i1++, i2++)
        for(size_t j1 = 0, j2 = jB; j1 < gc.size(); j1++, j2++)
            gc[i1][j1] = cd[i2][j2];
}

void join(vector<vector<uint64_t>>& gc, vector<vector<uint64_t>>& cd, size_t iB, size_t jB) {
    for(size_t i1 = 0, i2 = iB; i1 < gc.size(); i1++, i2++)
        for(size_t j1 = 0, j2 = jB; j1 < gc.size(); j1++, j2++)
            cd[i2][j2] = gc[i1][j1];
}

vector<vector<uint64_t>> strassen(vector<vector<uint64_t>> A, vector<vector<uint64_t>> B) {
    int n = A.size();

    if (n == 1) {
        vector<vector<uint64_t>> C(1, vector<uint64_t>(1));
        C[0][0] = A[0][0] * B[0][0];
        return C;
    }

    int newSize = n / 2;
    vector<vector<uint64_t>> a11(newSize, vector<uint64_t>(newSize)), a12(newSize, vector<uint64_t>(newSize)), a21(newSize, vector<uint64_t>(newSize)), a22(newSize, vector<uint64_t>(newSize));
    vector<vector<uint64_t>> b11(newSize, vector<uint64_t>(newSize)), b12(newSize, vector<uint64_t>(newSize)), b21(newSize, vector<uint64_t>(newSize)), b22(newSize, vector<uint64_t>(newSize));

    split(A, a11, 0 , 0); split(A, a12, 0 , newSize); split(A, a21, newSize, 0); split(A, a22, newSize, newSize);
    split(B, b11, 0 , 0); split(B, b12, 0 , newSize); split(B, b21, newSize, 0); split(B, b22, newSize, newSize);

    vector<vector<uint64_t>> p1 = strassen(add(a11, a22), add(b11, b22));
    vector<vector<uint64_t>> p2 = strassen(add(a21, a22), b11);
    vector<vector<uint64_t>> p3 = strassen(a11, sub(b12, b22));
    vector<vector<uint64_t>> p4 = strassen(a22, sub(b21, b11));
    vector<vector<uint64_t>> p5 = strassen(add(a11, a12), b22);
    vector<vector<uint64_t>> p6 = strassen(sub(a21, a11), add(b11, b12));
    vector<vector<uint64_t>> p7 = strassen(sub(a12, a22), add(b21, b22));

    vector<vector<uint64_t>> c11 = add(sub(add(p1, p4), p5), p7);
    vector<vector<uint64_t>> c12 = add(p3, p5);
    vector<vector<uint64_t>> c21 = add(p2, p4);
    vector<vector<uint64_t>> c22 = add(sub(add(p1, p3), p2), p6);

    vector<vector<uint64_t>> C(n, vector<uint64_t>(n));
    join(c11, C, 0 , 0); join(c12, C, 0 , newSize); join(c21, C, newSize, 0); join(c22, C, newSize, newSize);

    return C;
}

int main() {
    uint32_t n, i, j, d1, p1, r1, m1, d2, p2, r2, m2;
    cin >> n >> p1 >> d1 >> r1 >> m1 >> p2 >> d2 >> r2 >> m2;
    vector< vector<uint64_t> > A(n, vector<uint64_t>(n)), B(n, vector<uint64_t>(n));
    vector< vector<uint64_t> > C(n, vector<uint64_t>(n, 0));
    vector<uint64_t> V(n, 0);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            d1 = d1 * p1 + r1;
            d2 = d2 * p2 + r2;
            A[i][j] = d1 >> (32 - m1);
            B[i][j] = d2 >> (32 - m2); 
        }
    }

    C = strassen(A, B);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            V[i] ^= C[i][j];
        }
        cout << V[i] << " ";
    }

    return 0;
}
