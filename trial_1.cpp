#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

void vector_assign(vector<double>& D_0, vector<double>& D_p1, vector<double>& D_p2, 
                   vector<double>& D_n1, vector<double>& D_n2, int grids, int n1, vector<double>& B_source, double length) {

    double cndt, heat, length;
    double delta_z, delta_y;
    double ae, aw, an, as, ap;
    cout << "Enter the value of q :";
    cin >> heat;
    cout << "insert conductivity :";
    cin >> cndt;
    delta_z = length / double(grids);
    delta_y = delta_z;
    ae = cndt * delta_z / delta_y;
    aw = cndt * delta_z / delta_y;
    an = cndt * delta_y / delta_z;
    as = cndt * delta_y / delta_z;

    for ( int i =0; i < grids; ++i) {
        for(int j =0; j < grids; ++j) {
                    if (i == 0 && j==0) {
                        D_p1[i*grids+j] = -1 * ae;
                        D_p2[i*grids+j] = -1 * an;
                        D_n1[i*grids+j] = 0.0;
                        D_n2[i*grids+j] = 0.0;
                        D_0[i*grids+j] = ae + an;
                        B_source[i*grids+j] = -4.0 * heat * delta_z * delta_y/ length + heat * delta_z + heat * delta_y;
                    }
                    else if (i == 0 && j == grids - 1) {
                        D_p1[i*grids+j] = 0.0;
                        D_p2[i*grids+j] = -1 * an;
                        D_n1[i*grids+j] = -1 * aw;
                        D_n2[i*grids+j] = 0.0;
                        D_0[i*grids+j] = aw + an;
                        B_source[i*grids+j] = -4.0 * heat * delta_z * delta_y/ length + heat * delta_z + heat * delta_y;
                   }
                    else if (i == grids - 1 && j == 0) {
                        D_p1[i*grids+j] = -1 * ae;
                        D_p2[i*grids+j] = 0.0;
                        D_n1[i*grids+j] = 0.0;
                        D_n2[i*grids+j] = -1 * as;
                        D_0[i*grids+j] = ae + as;
                        B_source[i*grids+j] = -4.0 * heat * delta_z * delta_y/ length + heat * delta_z + heat * delta_y;
                    }
                    else if (i == grids - 1 && j == grids - 1) {
                        D_p1[i*grids+j] = 0.0;
                        D_p2[i*grids+j] = 0.0;
                        D_n1[i*grids+j] = -1 * aw;
                        D_n2[i*grids+j] = -1 * as;
                        D_0[i*grids+j] = aw + as;
                        B_source[i*grids+j] = -4.0 * heat * delta_z * delta_y/ length + heat * delta_y + heat * delta_z;
                    }
                    else if ( i == 0 && j!=0 && j != grids - 1) {
                        D_p1[i*grids+j] = -1 * ae;
                        D_p2[i*grids+j] = -1 * an;
                        D_n1[i*grids+j] = -1 * aw;
                        D_n2[i*grids+j] = 0.0;
                        D_0[i*grids+j] = ae + an + aw;
                        B_source[i*grids+j] = -4.0 * heat * delta_z * delta_y/ length + heat * delta_y;
                    }
                    else if ( j == 0 && i!=0 && i != grids - 1) {
                        D_p1[i*grids+j] = -1 * ae;
                        D_p2[i*grids+j] = -1 * an;
                        D_n1[i*grids+j] = 0.0;
                        D_n2[i*grids+j] = -1 * as;
                        D_0[i*grids+j] = ae + an + as;
                        B_source[i*grids+j] = -4.0 * heat * delta_z * delta_y/ length + heat * delta_z;
                    }
                    else if ( j == grids - 1 && i!=0 && i!= grids - 1) {
                        D_p1[i*grids+j] = 0.0;
                        D_p2[i*grids+j] = -1 * an;
                        D_n1[i*grids+j] = -1 * aw;
                        D_n2[i*grids+j] = -1 * as;
                        D_0[i*grids+j] = aw + an + as;
                        B_source[i*grids+j] = -4.0 * heat * delta_z * delta_y/ length + heat * delta_z;
                    }
                    else if ( i == grids - 1 && j!=0 && j != grids - 1) {
                        D_p1[i*grids+j] = -1 * ae;
                        D_p2[i*grids+j] = 0.0;
                        D_n1[i*grids+j] = -1 * aw;
                        D_n2[i*grids+j] = -1 * as;
                        D_0[i*grids+j] = ae + aw + as;
                        B_source[i*grids+j] = -4.0 * heat * delta_z * delta_y/ length + heat * delta_y;
                    }
                    else {
                        D_p1[i*grids+j] = -1 * ae;
                        D_p2[i*grids+j] = -1 * an;
                        D_n1[i*grids+j] = -1 * aw;
                        D_n2[i*grids+j] = -1 * as;
                        D_0[i*grids+j] = ae + an + aw + as;
                        B_source[i*grids+j] = -4.0 * heat * delta_z * delta_y/ length;
                    }
        }
    }

}

void vector_assemble(double** temperature, const vector<double>& D_0, 
                    const vector<double>& D_p1, const vector<double>& D_p2, 
                    const vector<double>& D_n1, const vector<double>& D_n2, 
                    int grids, int n1) {
                        int k=0, l=0, m=0, n=0, o=0;

                        for (int i = 0; i < n1; ++i){
                            for (int j = 0; j < n1; ++j) {
                                if ( i == j ) {
                                    temperature[i][j] = D_0[k];
                                    k++;
                                }
                                else if ( j == i + 1){
                                    temperature[i][j] = D_p1[l];
                                    l++;
                                }
                                else if ( j == i - 1){
                                    temperature[i][j] = D_n1[m+1];
                                    m++;
                                }
                                else if ( j == i + grids){
                                    temperature[i][j] = D_p2[n];
                                    n++;
                                }
                                else if ( j == i - grids){
                                    temperature[i][j] = D_n2[o+grids];
                                    o++;
                                }
                                else
                                    temperature[i][j] = 0.0;

                            }
                        }
}

void gauss_seidel(double** temperature, const std::vector<double>& B_source, int n1, vector<double>& x_result) {

    int max_iter = 10000, iter;
    double tolerance = 1e-6, error = 0.0;
    for (iter = 0; iter < max_iter; ++iter) {
        std::vector<double> x_old = x_result;

        for (int i = 0; i < n1; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n1; ++j) {
                if (j != i)
                    sum += temperature[i][j] * x_result[j];
            }
            x_result[i] = (B_source[i] - sum) / temperature[i][i];
        }

        error = 0.0;
        for (int i = 0; i < n1; ++i)
            error += std::abs(x_result[i] - x_old[i]);

        if (error < tolerance) {
            break;
        }
    }
    cout << "error is "<<error <<" at iteration "<<iter<< endl;
}

void result_print(vector<double>& x_result, int grids, double length, std::string& filename)) {   // vector<double>& x_result
    
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    double delta_x = length / double(grids);

    for (int i = 0; i < grids; ++i) {
        for (int j = 0; j < grids; ++j) {
            int idx = i * grids + j;
            file << i * delta_x << " " << j * delta_x << " " << x_result[idx] << "\n";
        }
        file << "\n";  
    }

    file.close();

}

int main() {
    int grids, n1, length;
    cout << "Enter the grid \nFor example \n10 for For 10 x 10 \n20 for 20 x 20 \n";
    cin >> grids;
    n1 = grids * grids;

    cout << "Enter the value of d :";
    cin >> length;

    double** temperature = new double*[n1];
    for (int i = 0; i < n1; i++) 
        temperature[i] = new double[n1];
    for ( int i =0; i < n1; ++i)
        for(int j =0; j < n1; ++j)
            temperature[i][j] = 0.0;
    vector<double> D_0(n1, 0.0);
    vector<double> D_p1(n1, 0.0);
    vector<double> D_p2(n1, 0.0);
    vector<double> D_n1(n1, 0.0);
    vector<double> D_n2(n1, 0.0);
    vector<double> B_source(n1, 0.0);
    vector<double> x_result(n1, 0.0);

    for ( int i =0; i < n1; ++i) {
        D_0[i] = 0.0;
        D_p1[i] = 0.0;
        D_p2[i] = 0.0;
        D_n1[i] = 0.0;
        D_n2[i] = 0.0;
    }

    vector_assign(D_0, D_p1, D_p2, D_n1, D_n2, grids, n1, B_source, length);
    vector_assemble(temperature, D_0, D_p1, D_p2, D_n1, D_n2, grids, n1);
    gauss_seidel(temperature, B_source, n1, x_result);
    result_print(x_result, grids, length, "temperature.dat");

    for (int i = 0; i < n1; i++) 
        delete[] temperature[i];
    return 0;
}

    

