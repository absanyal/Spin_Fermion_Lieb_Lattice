#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <complex>
#include "Matrix.h"
#include "diag.h"
#include "parameters.h"
#include "rngesus.hpp"
#include <iomanip>
#include <cassert>
#define PI acos(-1.0)

using namespace std;

parameters prm;
xorshift64 rng;

pair<int, int> kinv(int);
double lorentzian(double, double);
// cd phasex(int), phasey(int);

const double kB = 1.38 * pow(10, -23);
cd imagi = cd(0, 1);

typedef vector<cd> Mat_1_cd;
typedef vector<Mat_1_cd> Mat_2_cd;
typedef vector<Mat_2_cd> Mat_3_cd;

double mu;
double global_present_mu;
double temperature;
double measure_counter;
double theta1, theta2;

Matrix<cd> H, H_cluster, H_k;
Matrix<cd> H_up_k, H_dn_k; // these matrices used only if Spin is good quantum number.
Matrix<cd> T_intra, T_inter_px, T_inter_py;
Matrix<cd> theta;
Matrix<cd> phi;
Matrix<cd> eigs_k_all;
vector<double> eigs_, eigs_cluster_, eigs_k, eigs_up_k, eigs_dn_k;
Matrix<cd> SiSj;
Matrix<cd> Sq;

int ns;
int accepted = 0;
int total_change = 0;

Matrix<cd> sigma_x, sigma_y, sigma_z;

Mat_3_cd eigenstate_k;

pair<int, int> kinv(int M)
{
    int x = (int)M / (3 * prm.Ly);
    int y = M % (3 * prm.Ly);
    return std::make_pair(x, y);
}

int k(int x, int y)
{
    return y + prm.Ly * x;
}

// pair<int, int> kinv_cluster(int M)
// {
//     int x = (int)M / prm.c_Ly;
//     int y = M % prm.Ly;
//     return std::make_pair(x, y);
// }

// int k_cluster(int x, int y)
// {
//     return y + prm.c_Ly * x;
// }

double getmu()
{
    double target_N, n1;
    double mu1, mu2, mutemp, muout;
    bool converged;
    int nstates = int(2.0 * 3.0 * prm.Lx * prm.Ly);
    target_N = prm.filling * 2.0 * 3.0 * prm.Lx * prm.Ly;
    mu1 = eigs_[0];
    mu2 = eigs_[nstates - 1];
    mutemp = (eigs_[0] + eigs_[nstates - 1]) / 2.0;
    // mutemp = global_present_mu;
    for (int i = 0; i < 40000; i++)
    {
        n1 = 0.0;
        for (int j = 0; j < nstates; j++)
        {
            n1 +=
                    double(1.0 / (exp((eigs_[j] - mutemp) *
                                      (1.0 / temperature)) +
                                  1.0));
        }
        if (abs(target_N - n1) < double(0.00001))
        {
            converged = true;
            break;
        }
        else
        {
            if (n1 < target_N)
            {
                mu1 = mutemp;
                mutemp = 0.5 * (mutemp + mu2);
            }
            else
            {
                mu2 = mutemp;
                mutemp = 0.5 * (mutemp + mu1);
            }
        }
    }
    if (!converged)
    {
        cout << "mu not converged, stopping at N= " << n1 << endl;
    }
    global_present_mu = mutemp;
    return mutemp;
}

// double getmu_cluster()
// {
//     double target_N, n1;
//     double mu1, mu2, mutemp, muout;
//     bool converged;
//     int nstates = int(2.0 * prm.c_Lx * prm.c_Ly);
//     target_N = prm.filling * 2.0 * prm.c_Lx * prm.c_Ly;
//     mu1 = eigs_cluster_[0];
//     mu2 = eigs_cluster_[nstates - 1];
//     mutemp = (eigs_cluster_[0] + eigs_cluster_[nstates - 1]) / 2.0;
//     // mutemp = global_present_mu;
//     for (int i = 0; i < 40000; i++)
//     {
//         n1 = 0.0;
//         for (int j = 0; j < nstates; j++)
//         {
//             n1 +=
//                 double(1.0 / (exp((eigs_cluster_[j] - mutemp) *
//                                   (1.0 / temperature)) +
//                               1.0));
//         }
//         if (abs(target_N - n1) < double(0.00001))
//         {
//             converged = true;
//             break;
//         }
//         else
//         {
//             if (n1 < target_N)
//             {
//                 mu1 = mutemp;
//                 mutemp = 0.5 * (mutemp + mu2);
//             }
//             else
//             {
//                 mu2 = mutemp;
//                 mutemp = 0.5 * (mutemp + mu1);
//             }
//         }
//     }
//     if (!converged)
//     {
//         cout << "mu not converged, stopping at N= " << n1 << endl;
//     }
//     global_present_mu = mutemp;
//     return mutemp;
// }

double filter(double x)
{
    if (x > 2.0 * M_PI)
    {
        x -= 2.0 * M_PI;
    }
    else if (x < 0)
    {
        x += 2.0 * M_PI;
    }
    return x;
}

double getQuantumEnergy()
{
    double q_energy = 0;
    double this_mu = getmu();
    for (int i = 0; i < 2.0 * 3.0 * prm.Lx * prm.Ly; i++)
    {
        q_energy +=
                eigs_[i] * 1.0 /
                (1.0 + exp((eigs_[i] - this_mu) * (1.0 / temperature)));
    }
    return q_energy;
}

void makeTB()
{
    int i_new, ix_new, iy_new;
    // ix = 0;
    // iy = 0;
    int col_index, row_index;
    H.resize(2 * 3 * prm.Lx * prm.Ly, 2 * 3 * prm.Lx * prm.Ly);
    for (int ix = 0; ix < prm.Lx; ix++)
    {
        for (int iy = 0; iy < prm.Ly; iy++)
        {
            for (int alpha = 0; alpha < 3; alpha++)
            {
                for (int spin = 0; spin < 2; spin++)
                {
                    int i = iy + prm.Ly * ix;

                    //intra unit cell
                    col_index = i + (alpha * prm.Lx * prm.Ly) +
                            spin * 3 * prm.Lx * prm.Ly;

                    for (int beta = 0; beta < 3; beta++)
                    {
                        row_index = i + (beta * prm.Lx * prm.Ly) +
                                spin * 3 * prm.Lx * prm.Ly;
                        H(row_index, col_index) += T_intra(beta, alpha);
                        if (alpha == beta)
                        {
                            H(row_index, col_index) +=
                                    (2.0 * spin - 1.0) * prm.B;
                        }
                    }

                    //inter +x
                    ix_new = (ix + 1) % prm.Lx;
                    iy_new = iy;
                    i_new = iy_new + prm.Ly * ix_new;
                    col_index = i + (alpha * prm.Lx * prm.Ly) +
                            spin * 3 * prm.Lx * prm.Ly;
                    for (int beta = 0; beta < 3; beta++)
                    {
                        row_index = i_new + (beta * prm.Lx * prm.Ly) +
                                spin * 3 * prm.Lx * prm.Ly;
                        H(row_index, col_index) += T_inter_px(beta, alpha);
                        H(col_index, row_index) +=
                                conj(T_inter_px(beta, alpha));
                    }

                    //inter +y
                    ix_new = ix;
                    iy_new = (iy + 1) % prm.Ly;
                    i_new = iy_new + prm.Ly * ix_new;
                    col_index = i + (alpha * prm.Lx * prm.Ly) +
                            spin * 3 * prm.Lx * prm.Ly;
                    for (int beta = 0; beta < 3; beta++)
                    {
                        row_index = i_new + (beta * prm.Lx * prm.Ly) +
                                spin * 3 * prm.Lx * prm.Ly;
                        H(row_index, col_index) += T_inter_py(beta, alpha);
                        H(col_index, row_index) +=
                                conj(T_inter_py(beta, alpha));
                    }
                }
            }
        }
    }

    // Rashba SOC (STrictly for Lieb lattice)
    vector<double> bond_vector;
    bond_vector.resize(3);
    Matrix<cd> rashba_mat; //CONVENTION OF X,Y is Fig-1(a) of Physics Letters A 381 (2017) 944-948
    rashba_mat.resize(2, 2);
    for (int ix = 0; ix < prm.Lx; ix++)
    {
        for (int iy = 0; iy < prm.Ly; iy++)
        {
            for (int alpha = 0; alpha < 3; alpha++)
            {

                int i = iy + prm.Ly * ix;

                //intra unit cell
                for (int beta = 0; beta < 3; beta++)
                {
                    if (alpha == 0 && beta == 1)
                    {
                        bond_vector[0] = 1.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }
                    else if (alpha == 1 && beta == 0)
                    {
                        bond_vector[0] = -1.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }
                    else if (alpha == 0 && beta == 2)
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = -1.0;
                        bond_vector[2] = 0.0;
                    }
                    else if (alpha == 2 && beta == 0)
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = 1.0;
                        bond_vector[2] = 0.0;
                    }
                    else
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }

                    // sigma_x d_y - d_x sigma_y
                    for (int p = 0; p < 2; p++)
                    {
                        for (int q = 0; q < 2; q++)
                        {
                            rashba_mat(p, q) =
                                    sigma_x(p, q) * bond_vector[1] -
                                    sigma_y(p, q) * bond_vector[0];
                        }
                    }

                    for (int spin_i = 0; spin_i < 2; spin_i++)
                    {
                        for (int spin_j = 0; spin_j < 2; spin_j++)
                        {

                            col_index = i + (alpha * prm.Lx * prm.Ly) +
                                    spin_i * 3 * prm.Lx * prm.Ly;

                            row_index = i + (beta * prm.Lx * prm.Ly) +
                                    spin_j * 3 * prm.Lx * prm.Ly;
                            H(row_index, col_index) += imagi * prm.lambda_R * rashba_mat(spin_j, spin_i);
                            //cout<<col_index<<"  "<<row_index<<"  "<< imagi*prm.lambda_R * rashba_mat(spin_j, spin_i)<<endl;
                        }
                    }
                }

                // inter +X
                ix_new = (ix + 1) % prm.Lx;
                iy_new = iy;
                i_new = iy_new + prm.Ly * ix_new;
                for (int beta = 0; beta < 3; beta++)
                {

                    if (alpha == 1 && beta == 0)
                    {
                        bond_vector[0] = 1.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }
                    else
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }
                    for (int spin_i = 0; spin_i < 2; spin_i++)
                    {
                        for (int spin_j = 0; spin_j < 2; spin_j++)
                        {

                            // sigma_x d_y - d_x sigma_y
                            for (int p = 0; p < 2; p++)
                            {
                                for (int q = 0; q < 2; q++)
                                {
                                    rashba_mat(p, q) =
                                            sigma_x(p, q) * bond_vector[1] -
                                            sigma_y(p, q) * bond_vector[0];
                                }
                            }

                            col_index = i + (alpha * prm.Lx * prm.Ly) +
                                    spin_i * 3 * prm.Lx * prm.Ly;

                            row_index = i_new + (beta * prm.Lx * prm.Ly) +
                                    spin_j * 3 * prm.Lx * prm.Ly;
                            H(row_index, col_index) +=
                                    imagi *
                                    prm.lambda_R * rashba_mat(spin_j, spin_i);
                            H(col_index, row_index) +=
                                    conj(imagi *
                                         prm.lambda_R * rashba_mat(spin_j, spin_i));
                        }
                    }
                }

                // inter +Y
                ix_new = ix;
                iy_new = (iy + 1) % prm.Ly;
                i_new = iy_new + prm.Ly * ix_new;
                for (int beta = 0; beta < 3; beta++)
                {

                    if (alpha == 2 && beta == 0)
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = -1.0;
                        bond_vector[2] = 0.0;
                    }
                    else
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }
                    for (int spin_i = 0; spin_i < 2; spin_i++)
                    {
                        for (int spin_j = 0; spin_j < 2; spin_j++)
                        {

                            // sigma_x d_y - d_x sigma_y
                            for (int p = 0; p < 2; p++)
                            {
                                for (int q = 0; q < 2; q++)
                                {
                                    rashba_mat(p, q) =
                                            sigma_x(p, q) * bond_vector[1] -
                                            sigma_y(p, q) * bond_vector[0];
                                }
                            }

                            col_index = i + (alpha * prm.Lx * prm.Ly) +
                                    spin_i * 3 * prm.Lx * prm.Ly;

                            row_index = i_new + (beta * prm.Lx * prm.Ly) +
                                    spin_j * 3 * prm.Lx * prm.Ly;
                            H(row_index, col_index) +=
                                    imagi *
                                    prm.lambda_R * rashba_mat(spin_j, spin_i);
                            H(col_index, row_index) +=
                                    conj(imagi *
                                         prm.lambda_R * rashba_mat(spin_j, spin_i));
                            //cout<<"here 1"<<imagi*prm.lambda_R * rashba_mat(spin_j, spin_i)<<endl;
                        }
                    }
                }
            }
        }
    }
}

void initialize_T()
{
    //intra
    T_intra.resize(3, 3);
    T_intra(0, 1) = prm.t;
    T_intra(0, 2) = prm.t;
    T_intra(1, 0) = conj(T_intra(0, 1));
    T_intra(2, 0) = conj(T_intra(0, 2));

    //inter +x
    T_inter_px.resize(3, 3);
    T_inter_px(0, 1) = prm.t;

    //inter +y
    T_inter_py.resize(3, 3);
    T_inter_py(0, 2) = prm.t;
}

// void makeTB_cluster()
// {
//     H_cluster.resize(2 * prm.c_Lx * prm.c_Ly, 2 * prm.c_Lx * prm.c_Ly);
//     int ns_c;
//     ns_c = prm.c_Lx * prm.c_Ly;
//     for (int i = 0; i < prm.c_Lx; i++)
//     {
//         for (int j = 0; j < prm.c_Ly; j++)
//         {
//             // cout << "TB step " << i << "," << j << endl;
//             int pxy = k_cluster(i, j);

//             int pxp1y = k_cluster((i + 1) % prm.c_Lx, j);
//             int pxm1y = k_cluster((i - 1 + prm.c_Lx) % prm.c_Lx, j);

//             int pxyp1 = k_cluster(i, (j + 1) % prm.c_Ly);
//             int pxym1 = k_cluster(i, (j - 1 + prm.c_Ly) % prm.c_Ly);

//             int pxp1ym1 = k_cluster((i + 1) % prm.c_Lx,
//                                     (j - 1 + prm.c_Ly) % prm.c_Ly);
//             int pxm1yp1 = k_cluster((i - 1 + prm.c_Lx) % prm.c_Lx,
//                                     (j + 1) % prm.c_Ly);

//             H_cluster(pxp1y, pxy) = cd(prm.tx, 0.0);
//             H_cluster(pxm1y, pxy) = cd(prm.tx, 0.0);

//             H_cluster(pxyp1, pxy) = cd(prm.ty, 0.0);
//             H_cluster(pxym1, pxy) = cd(prm.ty, 0.0);

//             H_cluster(pxp1ym1, pxy) = cd(prm.td, 0.0);
//             H_cluster(pxm1yp1, pxy) = cd(prm.td, 0.0);

//             H_cluster(pxp1y + ns_c, pxy + ns_c) = cd(prm.tx, 0.0);
//             H_cluster(pxm1y + ns_c, pxy + ns_c) = cd(prm.tx, 0.0);

//             H_cluster(pxyp1 + ns_c, pxy + ns_c) = cd(prm.ty, 0.0);
//             H_cluster(pxym1 + ns_c, pxy + ns_c) = cd(prm.ty, 0.0);

//             H_cluster(pxp1ym1 + ns_c, pxy + ns_c) = cd(prm.td, 0.0);
//             H_cluster(pxm1yp1 + ns_c, pxy + ns_c) = cd(prm.td, 0.0);
//         }
//     }
// // }

void makeHund()
{
    // Hund's terms
    int index, i;
    int ns_temp = 3.0 * prm.Lx * prm.Ly;

    for (int ix = 0; ix < prm.Lx; ix++)
    {
        for (int iy = 0; iy < prm.Ly; iy++)
        {
            for (int alpha = 0; alpha < 3; alpha++)
            {

                i = ix + prm.Lx * iy;
                index = i + (alpha * prm.Lx * prm.Ly);

                H(index, index) +=
                        0.5 * prm.JH * cos(theta(alpha, index));
                H(index + ns_temp, index + ns_temp) +=
                        -0.5 * prm.JH * cos(theta(alpha, index));
                H(index, index + ns_temp) += 0.5 *
                        prm.JH * sin(theta(0, index)) *
                        exp(-imagi * phi(0, index));
                H(index + ns_temp, index) += 0.5 * prm.JH *
                        sin(theta(alpha, index)) *
                        exp(imagi * phi(0, index));
            }
        }
    }
}

// void makeHund_cluster(int center_site)
// {
//     // theta.print();
//     int i_x, i_y, i_new, i_new_x, i_new_y, cx, cy;

//     int ns_c;
//     ns_c = prm.c_Lx * prm.c_Ly;

//     cx = kinv(center_site).first;
//     cy = kinv(center_site).second;

//     for (int i = 0; i < prm.c_Lx * prm.c_Ly; i++)
//     {
//         i_x = kinv_cluster(i).first;
//         i_y = kinv_cluster(i).second;

//         i_new_x = cx - (prm.c_Lx / 2) + i_x;
//         i_new_y = cy - (prm.c_Ly / 2) + i_y;

//         i_new_x = (i_new_x + prm.Lx) % prm.Lx;
//         i_new_y = (i_new_y + prm.Ly) % prm.Ly;

//         i_new = k(i_new_x, i_new_y);
//         // cout << i_new << " ";
//         // if ((i+1) % prm.c_Ly == 0)
//         // {
//         //     cout << endl;
//         // }
//         // cout << i << " ";
//         // if ((i+1) % prm.c_Ly == 0)
//         // {
//         //     cout << endl;
//         // }
//         // assert(i_new != 0);

//         H_cluster(i, i) = 0.5 * prm.JH * cos(theta(0, i_new));
//         H_cluster(i + ns_c, i + ns_c) = -0.5 * prm.JH * cos(theta(0, i_new));
//         H_cluster(i, i + ns_c) =
//             0.5 * prm.JH * sin(theta(0, i_new)) * exp(-imagi * phi(0, i_new));
//     }
// }

double perturb()
{
    double p;
    p = (rng.random() * 2.0 * prm.window - prm.window) * 2.0 * M_PI;
    return p;
}

double lnP(int i)
{
    double sum = 0.0;
    double mu = getmu();
    for (int lambda = 0; lambda < eigs_.size(); lambda++)
    {
        sum += log((1.0 + exp(-(1 / (temperature)) * (eigs_[lambda] - mu))));
        // cout << sum << endl;
    }
    return sum;
}

// double lnP_cluster(int i)
// {
//     double sum = 0.0;
//     double mu = getmu_cluster();
//     // cout << mu << endl;
//     for (int lambda = 0; lambda < eigs_cluster_.size(); lambda++)
//     {
//         sum +=
//             log((1.0 + exp(-(1 / (temperature)) * (eigs_cluster_[lambda] -
//                                                    mu))));
//         // cout << sum << endl;
//     }
//     return sum;
// }

double lorentzian(double x, double x0)
{
    return (1.0 / 3.1415926535) *
            0.5 * prm.G / (pow(x - x0, 2) + pow(0.5 * prm.G, 2));
}

void make_TB_k(double kx, double ky)
{
    double Tx, Ty, TI, Rx, Ry, sx, sy, cx, cy, Splus, Sminus;
    H_k.resize(2 * 3, 2 * 3);
    sx = sin(kx / 2.0);
    sy = sin(ky / 2.0);
    cx = cos(kx / 2.0);
    cy = cos(ky / 2.0);
    Tx = 2.0 * prm.t * cx;
    Ty = 2.0 * prm.t * cy;
    Rx = 2.0 * prm.lambda_R * sx;
    Ry = 2.0 * prm.lambda_R * sy;
    TI = 4.0 * prm.lambda * sx * sy;
    Splus = (prm.Delta_c + prm.Delta_s);
    Sminus = (prm.Delta_c - prm.Delta_s);

    // Rx = 0;
    // Ry = 0;
    // TI = 0;

    H_k(0, 0) = prm.B + Splus;
    H_k(0, 1) = Tx;
    H_k(0, 2) = Ty;
    H_k(0, 3) = 0;
    H_k(0, 4) = -imagi * Rx;
    H_k(0, 5) = -Ry;

    H_k(1, 0) = Tx;
    H_k(1, 1) = prm.B - Splus;
    H_k(1, 2) = -imagi * TI;
    H_k(1, 3) = -imagi * Rx;
    H_k(1, 4) = 0;
    H_k(1, 5) = 0;

    H_k(2, 0) = Ty;
    H_k(2, 1) = imagi * TI;
    H_k(2, 2) = prm.B - Splus;
    H_k(2, 3) = -Ry;
    H_k(2, 4) = 0;
    H_k(2, 5) = 0;

    H_k(3, 0) = 0;
    H_k(3, 1) = imagi * Rx;
    H_k(3, 2) = -Ry;
    H_k(3, 3) = -prm.B + Sminus;
    H_k(3, 4) = Tx;
    H_k(3, 5) = Ty;

    H_k(4, 0) = imagi * Rx;
    H_k(4, 1) = 0;
    H_k(4, 2) = 0;
    H_k(4, 3) = Tx;
    H_k(4, 4) = -prm.B - Sminus;
    H_k(4, 5) = imagi * TI;

    H_k(5, 0) = -Ry;
    H_k(5, 1) = 0;
    H_k(5, 2) = 0;
    H_k(5, 3) = Ty;
    H_k(5, 4) = -imagi * TI;
    H_k(5, 5) = -prm.B - Sminus;
}

void make_TB_k_spin_reoslved(double kx, double ky)
{
    double Tx, Ty, TI, sx, sy, cx, cy, Splus, Sminus;
    H_up_k.resize(3, 3);
    H_dn_k.resize(3, 3);
    sx = sin(kx / 2.0);
    sy = sin(ky / 2.0);
    cx = cos(kx / 2.0);
    cy = cos(ky / 2.0);
    Tx = 2.0 * prm.t * cx;
    Ty = 2.0 * prm.t * cy;
    TI = 4.0 * prm.lambda * sx * sy;
    Splus = (prm.Delta_c + prm.Delta_s);
    Sminus = (prm.Delta_c - prm.Delta_s);

    // Rx = 0;
    // Ry = 0;
    // TI = 0;

    H_up_k(0, 0) = prm.B + Splus;
    H_up_k(0, 1) = Tx;
    H_up_k(0, 2) = Ty;


    H_up_k(1, 0) = Tx;
    H_up_k(1, 1) = prm.B - Splus;
    H_up_k(1, 2) = -imagi * TI;


    H_up_k(2, 0) = Ty;
    H_up_k(2, 1) = imagi * TI;
    H_up_k(2, 2) = prm.B - Splus;


    H_dn_k(0, 0) = -prm.B + Sminus;
    H_dn_k(0, 1) = Tx;
    H_dn_k(0, 2) = Ty;

    H_dn_k(1, 0) = Tx;
    H_dn_k(1, 1) = -prm.B - Sminus;
    H_dn_k(1, 2) = imagi * TI;

    H_dn_k(2, 0) = Ty;
    H_dn_k(2, 1) = -imagi * TI;
    H_dn_k(2, 2) = -prm.B - Sminus;
}

void initialize_sigma()
{
    sigma_x.resize(2, 2);
    sigma_y.resize(2, 2);
    sigma_z.resize(2, 2);

    // X
    sigma_x(0, 0) = 0.0;
    sigma_x(0, 1) = 1.0;
    sigma_x(1, 0) = 1.0;
    sigma_x(1, 1) = 0.0;

    // y
    sigma_y(0, 0) = 0.0;
    sigma_y(0, 1) = -imagi;
    sigma_y(1, 0) = imagi;
    sigma_y(1, 1) = 0.0;

    // Z
    sigma_z(0, 0) = 1.0;
    sigma_z(0, 1) = 0.0;
    sigma_z(1, 0) = 0.0;
    sigma_z(1, 1) = -1.0;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Enter name of input file";
    }
    string inputfile = argv[1];
    prm.load(inputfile);

    ns = 3.0 * prm.Lx * prm.Ly;

    rng.set_seed(prm.seed);

    initialize_T();
    initialize_sigma();

    temperature = prm.T;
    double totalQE;
    // cout << "Stopping T = " << prm.T_stop << endl;
    double epsilon = pow(10, -6);
    double theta_old, theta_new, P_old, P_new, P_ratio, r;
    int pos;

    //////////////////////TEST STUFF HERE//////////////////////////////
    bool Spin_resolved=true;
    if(Spin_resolved){
        assert(prm.lambda_R==0.0);
    }
    string outfname = "data/dos.dat";
    ofstream dosfile;
    dosfile.open(outfname);
    if (prm.k_only == 0)
    {
        theta.resize(3, ns);
        phi.resize(3, ns);
        SiSj.resize(ns, ns);
        // Initially fill theta with random values
        for (int pos = 0; pos < ns; pos++)
        {
            for (int orbital = 0; orbital < 3; orbital++)
            {

                theta(orbital, pos) = rng.random() * 2.0 * M_PI;
                phi(orbital, pos) = 0.0;
            }
        }
        cout << "Populated initial random angles" << endl;

        cout << "Populated random theta array." << endl;
        cout << "Real space calculations" << endl;
        makeTB();
        //makeHund();
        cout << H.n_row() << ", " << H.n_col() << endl;
        //        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx"<<endl;
        //        H.print();
        //        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXx"<<endl;

        Diagonalize('N', H, eigs_);

        cout << eigs_.size() << endl;

        // int counter = 0;
        double full_rho = 0;
        for (int xi = 0; xi <= prm.numw; xi++)
        {
            double x = prm.minw + xi * ((prm.maxw - prm.minw) / prm.numw);
            // cout << x;

            for (int i = 0; i < eigs_.size(); i++)
            {
                full_rho += (1.0 / eigs_.size()) * lorentzian(x, eigs_[i]);
            }
            dosfile << x << "\t" << full_rho << endl;
            full_rho = 0;
        }
    }
    dosfile.close();

    // double ky = M_PI;

    ofstream bandfile;
    bandfile.open("data/bands.txt");
    ofstream band_ky_pi;
    band_ky_pi.open("data/band_ky_pi");

    eigs_k_all.resize(prm.Lx * prm.Ly, 2 * 3);

    eigenstate_k.resize(6);
    for (int band = 0; band < 6; band++)
    {
        eigenstate_k[band].resize(prm.Lx * prm.Ly);
        for (int momenta = 0; momenta < prm.Lx * prm.Ly; momenta++)
        {
            eigenstate_k[band][momenta].resize(6);
        }
    }

    Matrix<cd> F_mat; //F1, F2, F3, F4, F5;
    F_mat.resize(6, prm.Lx * prm.Ly);

    Matrix<cd> F_mat_orgnl; //F1, F2, F3, F4, F5;
    F_mat_orgnl.resize(6, prm.Lx * prm.Ly);

    cout << "plotting bands" << endl;

    int row_counter = 0;
    for (int nx = 0; nx < prm.Lx; nx++)
    {
        for (int ny = 0; ny < prm.Ly; ny++)
        {
            double kx = (2 * M_PI * nx / (1.0 * prm.Lx));
            double ky = (2 * M_PI * ny / (1.0 * prm.Ly));
            double n = ny + prm.Ly * nx;


            if(!Spin_resolved){
            make_TB_k(kx, ky);
            Diagonalize('V', H_k, eigs_k);
            for (int band = 0; band < 6; band++)
            {
                for (int component = 0; component < 6; component++)
                {
                    eigenstate_k[band][n][component] = H_k(component, band);
                }
            }
            }
            else{
                eigs_k.resize(6);
                make_TB_k_spin_reoslved(kx, ky);
                Diagonalize('V', H_up_k, eigs_up_k);
                Diagonalize('V', H_dn_k, eigs_dn_k);
                for (int band = 0; band < 3; band++)
                {
                    for (int component = 0; component < 3; component++)
                    {
                        eigenstate_k[2*band][n][component] = H_up_k(component, band);
                        eigenstate_k[2*band][n][component+3] = 0.0;
                        eigenstate_k[(2*band)+1][n][component+3] = H_dn_k(component, band);
                        eigenstate_k[(2*band)+1][n][component] =0.0;

                    }
                    eigs_k[2*band]=eigs_up_k[band];
                    eigs_k[(2*band) + 1]=eigs_dn_k[band];
                }

            }





            cd val;
            val = 0;
            for (int band = 0; band < 6; band++)
            {
                for (int component = 0; component < 6; component++)
                {
                    val += conj(eigenstate_k[band][n][component]) * eigenstate_k[band][n][component];
                }
            }
            val = sqrt(abs(val));

            for (int band = 0; band < 6; band++)
            {
                for (int component = 0; component < 6; component++)
                {
                    eigenstate_k[band][n][component] = (1.0 / val) * eigenstate_k[band][n][component];
                }
            }

            // eigenstate_k_0.print();

            bandfile << kx / M_PI << "\t" << ky / M_PI;
            if (ky / M_PI == prm.fixed_ky)
            {
                band_ky_pi << kx / M_PI;
            }
            for (int i = 0; i < eigs_k.size(); i++)
            {
                // cout << "\t" << eigs_k[i];
                bandfile << "\t" << eigs_k[i];
                eigs_k_all(row_counter, i) = real(eigs_k[i]);
                if (ky / M_PI == prm.fixed_ky)
                {
                    band_ky_pi << "\t" << eigs_k[i];
                }
            }
            row_counter += 1;
            // cout << endl;
            bandfile << endl;
        }
        band_ky_pi << endl;
        bandfile << endl;
    }

    // Fixing global phase
    for (int nx = 0; nx < prm.Lx; nx++)
    {
        for (int ny = 0; ny < prm.Ly; ny++)
        {
            int n = ny + prm.Ly * nx;
            for (int band = 0; band < 6; band++)
            {
                cd phase = 0.0;
                // for (int comp = 0; comp < 6; comp++)
                // {
                //     phase += eigenstate_k[band][n][comp];
                // }
                // phase = phase * (1.0 / (abs(phase)));
                // cout << n << " " << band << " " << phase << endl;
                phase = eigenstate_k[5][35][0];
                phase = phase * (1.0 / abs(phase));
                //cout << n << " " << band << " " << phase << endl;
                for (int comp = 0; comp < 6; comp++)
                {
                    // eigenstate_k[band][n][comp] =
                    //     conj(phase) *
                    //     eigenstate_k[band][n][comp];
                    //                    eigenstate_k[band][n][comp] =
                    //                        conj(phase) *
                    //                        eigenstate_k[band][n][comp];
                }
            }
        }
    }

    // eigs_k_all.print();
    cd Ux_k, Uy_k, Ux_kpy, Uy_kpx;
    vector<cd> F_bands;
    F_bands.resize(6);
    vector<cd> Chern_num;
    Chern_num.resize(6);

    cd nkp1_nkp1p2, nkp1_nkp1, nk_nkp2, nk_nk, nkp2_nkp1p2, nkp2_nkp2, nk_nkp1;
    vector<cd> F_bands_orgnl;
    F_bands_orgnl.resize(6);
    vector<cd> Chern_num_orgnl;
    Chern_num_orgnl.resize(6);
    for (int band = 0; band < 6; band++)
    {
        string file_Fk="Fk_band"+to_string(band)+".txt";
        ofstream fl_Fk_out(file_Fk.c_str());
        fl_Fk_out<<"#nx  ny  tilde_F(nx,ny).real()  tilde_F(nx,ny).imag()  ArgofLog.real()  ArgofLog.imag()"<<endl;
        fl_Fk_out<<"#Extra momentum point for pm3d corners2color c1"<<endl;

        string file_Fk_orgnl="Fk_original_band"+to_string(band)+".txt";
        ofstream fl_Fk_orgnl_out(file_Fk_orgnl.c_str());
        fl_Fk_orgnl_out<<"#nx  ny  F(nx,ny).real()*(2pi/Lx)*(2pi/Ly)  F(nx,ny).imag()*(2pi/Lx)*(2pi/Ly)"<<endl;

        F_bands[band] = 0.0;
        F_bands_orgnl[band] = 0.0;
        for (int nx = 0; nx < prm.Lx; nx++)
        {
            for (int ny = 0; ny < prm.Ly; ny++)
            {
                int n = ny + prm.Ly * nx;
                int n_left, n_right, nx_left, ny_left, nx_right, ny_right;

                //U1_k
                Ux_k = 0;
                n_left = n;
                nx_right = (nx + 1) % prm.Lx;
                ny_right = ny;
                n_right = ny_right + prm.Ly * nx_right;
                for (int comp = 0; comp < 6; comp++)
                {
                    Ux_k +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }
                Ux_k = Ux_k * (1.0 / abs(Ux_k));

                //U2_kpx
                Uy_kpx = 0;
                nx_left = (nx + 1) % prm.Lx;
                ny_left = ny;
                n_left = ny_left + prm.Ly * nx_left;
                nx_right = nx_left;
                ny_right = (ny_left + 1) % prm.Ly;
                n_right = ny_right + prm.Ly * nx_right;
                for (int comp = 0; comp < 6; comp++)
                {
                    Uy_kpx +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }
                Uy_kpx = Uy_kpx * (1.0 / abs(Uy_kpx));

                //U1_kpy
                Ux_kpy = 0;
                nx_left = nx;
                ny_left = (ny + 1) % prm.Ly;
                n_left = ny_left + prm.Ly * nx_left;
                nx_right = (nx_left + 1) % prm.Lx;
                ny_right = ny_left;
                n_right = ny_right + prm.Ly * nx_right;
                for (int comp = 0; comp < 6; comp++)
                {
                    Ux_kpy +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }
                Ux_kpy = Ux_kpy * (1.0 / abs(Ux_kpy));

                //U2_k
                Uy_k = 0;
                nx_left = nx;
                ny_left = ny;
                n_left = ny_left + prm.Ly * nx_left;
                nx_right = nx_left;
                ny_right = (ny_left + 1) % prm.Ly;
                n_right = ny_right + prm.Ly * nx_right;
                for (int comp = 0; comp < 6; comp++)
                {
                    Uy_k +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }
                Uy_k = Uy_k * (1.0 / abs(Uy_k));

                // Calculating tilde F12
                F_mat(band, n) = log(Ux_k *
                                     Uy_kpx *
                                     conj(Ux_kpy) * conj(Uy_k));

                F_bands[band] += F_mat(band, n);




                //Original Fk
                //nkp1_nkp1p2
                nkp1_nkp1p2=0.0;
                nx_left = (nx + 1)%prm.Lx;
                ny_left = ny;
                n_left = ny_left + prm.Ly * nx_left;
                nx_right = (nx + 1)%prm.Lx;
                ny_right = (ny + 1)%prm.Ly;
                n_right = ny_right + prm.Ly * nx_right;
                for (int comp = 0; comp < 6; comp++)
                {
                    nkp1_nkp1p2 +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }

                //nkp1_nkp1
                nkp1_nkp1=0.0;
                nx_left = (nx + 1)%prm.Lx;
                ny_left = ny;
                n_left = ny_left + prm.Ly * nx_left;
                n_right = n_left;
                for (int comp = 0; comp < 6; comp++)
                {
                    nkp1_nkp1 +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }

                // nk_nkp2
                nk_nkp2=0.0;
                nx_left = nx;
                ny_left = ny;
                n_left = ny_left + prm.Ly * nx_left;
                nx_right = nx;
                ny_right = (ny + 1)%prm.Ly;
                n_right = ny_right + prm.Ly * nx_right;
                for (int comp = 0; comp < 6; comp++)
                {
                    nk_nkp2 +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }

                // nk_nk
                nk_nk=0.0;
                n_left = n;
                n_right = n;
                for (int comp = 0; comp < 6; comp++)
                {
                    nk_nk +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }

                // nkp2_nkp1p2
                nkp2_nkp1p2=0.0;
                nx_left = nx;
                ny_left = (ny+1)%prm.Ly;
                n_left = ny_left + prm.Ly * nx_left;
                nx_right = (nx + 1)%prm.Lx;
                ny_right = (ny + 1)%prm.Ly;
                n_right = ny_right + prm.Ly * nx_right;
                for (int comp = 0; comp < 6; comp++)
                {
                    nkp2_nkp1p2 +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }



                // nkp2_nkp2
                nkp2_nkp2=0.0;
                nx_left = nx;
                ny_left = (ny+1)%prm.Ly;
                n_left = ny_left + prm.Ly * nx_left;
                nx_right = nx;
                ny_right = (ny + 1)%prm.Ly;
                n_right = ny_right + prm.Ly * nx_right;
                for (int comp = 0; comp < 6; comp++)
                {
                    nkp2_nkp2 +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }


                // nk_nkp1
                nk_nkp1=0.0;
                nx_left = nx;
                ny_left = ny;
                n_left = ny_left + prm.Ly * nx_left;
                nx_right = (nx + 1)%prm.Lx;
                ny_right = ny;
                n_right = ny_right + prm.Ly * nx_right;
                for (int comp = 0; comp < 6; comp++)
                {
                    nk_nkp1 +=
                            conj(eigenstate_k[band][n_left][comp]) *
                            eigenstate_k[band][n_right][comp];
                }

                //Original F12 calculation
                F_mat_orgnl(band, n) = (   (nkp1_nkp1p2 - nkp1_nkp1)*((1.0*prm.Ly)/(2.0*PI))
                                         + (-nk_nkp2 + nk_nk)*((1.0*prm.Ly)/(2.0*PI))
                                         + (-nkp2_nkp1p2 + nkp2_nkp2)*((1.0*prm.Lx)/(2.0*PI))
                                         + (nk_nkp1 - nk_nk)*((1.0*prm.Lx)/(2.0*PI))
                                        )*((2.0*M_PI)/(1.0*prm.Ly))*((2.0*M_PI)/(1.0*prm.Lx));

                F_bands_orgnl[band] += F_mat_orgnl(band, n);




                fl_Fk_out<<nx<<"  "<<ny<<"  "<<F_mat(band, n).real()<<"  "<<F_mat(band, n).imag()<<
                           "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
                           "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;


                fl_Fk_orgnl_out<<nx<<"  "<<ny<<"  "<<F_mat_orgnl(band, n).real()<<"  "<<F_mat_orgnl(band, n).imag()<<endl;

                if(ny==prm.Ly-1){//For pm3d corners2color c1
                    fl_Fk_out<<nx<<"  "<<ny<<"  "<<F_mat(band, n).real()<<"  "<<F_mat(band, n).imag()<<
                               "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
                               "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;
                }
            }
            if(nx==prm.Lx-1){//For pm3d corners2color c1
                fl_Fk_out<<endl;
                for(int ny_=0;ny_<prm.Ly;ny_++){
                    int n_ = ny_ + prm.Ly * nx;
                    fl_Fk_out<<nx<<"  "<<ny_<<"  "<<F_mat(band, n_).real()<<"  "<<F_mat(band, n_).imag()<<
                               "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
                               "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;
                }
            }
            fl_Fk_out<<endl;
            fl_Fk_orgnl_out<<endl;
        }


        Chern_num[band] = (-1.0 * imagi / (2 * M_PI)) * F_bands[band];
        Chern_num_orgnl[band] = (-1.0 * imagi / (2 * M_PI)) * F_bands_orgnl[band];
        fl_Fk_out<<"#Chern no*2pi*Iota= "<<F_bands[band].real()<<"  "<<F_bands[band].imag()<<endl;
        fl_Fk_orgnl_out<<"#Sum_{nx,ny}F(nx,ny)*(2pi/Lx)*(2pi/Ly) = "<<F_bands_orgnl[band].real()<<"  "<<F_bands_orgnl[band].imag()<<endl;
        fl_Fk_orgnl_out<<"# Lx and Ly are number of unit cells in x and y direction respectively, each unit cell has three atoms"<<endl;
        cout << "tilde Chern number [" << band << "] = " << Chern_num[band].real() << " " << Chern_num[band].imag() << endl;
        cout << "Chern number [" << band << "] = " << Chern_num_orgnl[band].real() << " " << Chern_num_orgnl[band].imag() << endl;

    }

    bandfile.close();
    band_ky_pi.close();
    cout << "Bands completed" << endl;

    outfname = "data/band_path.dat";
    ofstream bandpath;
    bandpath.open(outfname);

    // Plot along paths
    int path_point_counter = 0;

    // Gamma(0,0) to X(pi,0)
    for (int nx = 0; nx < prm.Lx / 2; nx++)
    {
        double kx = (2.0 * M_PI * nx / (1.0 * prm.Lx));
        double ky = 0.0;
        make_TB_k(kx, ky);
        // H_k.print();
        // cin >> tempnumber;
        Diagonalize('N', H_k, eigs_k);
        bandpath << path_point_counter << "\t"
                 << kx / M_PI << "\t" << ky / M_PI;
        for (int i = 0; i < eigs_k.size(); i++)
        {
            // cout << "\t" << eigs_k[i];
            bandpath << "\t" << eigs_k[i];
        }
        bandpath << endl;
        path_point_counter += 1;
    }

    // X(pi,0) to M(pi, pi)
    for (int ny = 0; ny < prm.Ly / 2; ny++)
    {
        double kx = M_PI;
        double ky = (2.0 * M_PI * ny / (1.0 * prm.Ly));
        make_TB_k(kx, ky);
        // H_k.print();
        // cin >> tempnumber;
        Diagonalize('N', H_k, eigs_k);
        bandpath << path_point_counter << "\t"
                 << kx / M_PI << "\t" << ky / M_PI;
        for (int i = 0; i < eigs_k.size(); i++)
        {
            // cout << "\t" << eigs_k[i];
            bandpath << "\t" << eigs_k[i];
        }
        bandpath << endl;
        path_point_counter += 1;
    }

    // M(pi, pi) to Gamma(0,0)
    for (int ny = prm.Ly / 2; ny >= 0; ny--)
    {
        double ky = (2.0 * M_PI * ny / (1.0 * prm.Ly));
        double kx = ky;
        make_TB_k(kx, ky);
        // H_k.print();
        // cin >> tempnumber;
        Diagonalize('N', H_k, eigs_k);
        bandpath << path_point_counter << "\t"
                 << kx / M_PI << "\t" << ky / M_PI;
        for (int i = 0; i < eigs_k.size(); i++)
        {
            // cout << "\t" << eigs_k[i];
            bandpath << "\t" << eigs_k[i];
        }
        bandpath << endl;
        path_point_counter += 1;
    }

    cout << "Plotted along paths" << endl;
    bandpath.close();

    outfname = "data/dos_k.dat";
    dosfile.open(outfname);

    // int counter = 0;
    for (int xi = 0; xi <= prm.numw; xi++)
    {
        double x = prm.minw + xi * ((prm.maxw - prm.minw) / prm.numw);
        dosfile << x;

        for (int ec = 0; ec < 6; ec++)
        {
            double rho_ec = 0;
            for (int rc = 0; rc < prm.Lx * prm.Ly; rc++)
            {
                rho_ec += (1.0 / (6 * prm.Lx * prm.Ly)) * lorentzian(x, real(eigs_k_all(rc, ec)));
            }
            dosfile << "\t" << rho_ec;
        }
        dosfile << endl;
    }

    // dosfile.close();

    // dosfile << x << '\t' << fullrho << endl;

    dosfile.close();

    cout << "DoS plotted" << endl;

    //______________________________________________________________________
    // cout << "Exiting without running MC loops" << endl;

    return 0;
    ///////////////////////////////////////////////////////////////////

    // TB Part
    if (prm.TCA == 0)
    {
        makeTB();

        cout << "Populated TB matrix." << endl;

        measure_counter = 0;

        // ofstream tvsE;
        // tvsE.open("data/T_vs_E.txt");

        cout << "Performing ED" << endl;
        while (temperature >= prm.T_stop - epsilon)
        {
            totalQE = 0;
            for (int t = 0; t < prm.sweeps; t++)
            {
                total_change = 0;
                accepted = 0;
                // cout << "Sweep number: " << t + 1 << endl;
                for (int i = 0; i < prm.Lx; i++)
                {
                    for (int j = 0; j < prm.Ly; j++)
                    {
                        for (int orbital = 0; orbital < 3; orbital++)
                        {
                            // cout << "Sweep: " << t + 1
                            //      << " (" << i + 1 << ", " << j + 1 << ")\n";
                            // double theta_old, theta_new, P_old, P_new, P_ratio, r;
                            // int pos;
                            pos = k(i, j);

                            theta_old = real(theta(orbital, pos));
                            makeTB();
                            makeHund();
                            Diagonalize('V', H, eigs_);
                            P_old = lnP(pos);

                            theta_new = filter(theta_old + perturb());
                            theta(orbital, pos) = theta_new;
                            makeTB();
                            makeHund();
                            Diagonalize('V', H, eigs_);
                            P_new = lnP(pos);

                            r = rng.random();
                            P_ratio = exp(P_new - P_old);

                            if (r < P_ratio)
                            {
                                // cout << "Accepted" << endl;
                                accepted += 1;
                                theta(orbital, pos) = theta_new;
                            }
                            else
                            {
                                // cout << "Rejected" << endl;
                                theta(orbital, pos) = theta_old;
                            }
                            total_change += 1;
                        }
                    }
                }

                // Quantum energy measurement only for ED
                double presentQE = getQuantumEnergy();
                totalQE += presentQE;
                cout << "T: " << temperature << ", SN: " << t + 1 << ", AR: "
                     << (accepted * 1.0) / (total_change * 1.0) << ", En: "
                     << presentQE << ", mu: " << getmu() << endl;

                // Do measurements if thermalized
                if (abs(double(temperature - prm.T_stop)) < pow(10, -4))
                {
                    if ((t + 1) >= prm.measure_after &&
                            (t - int(prm.measure_after) + 1) %
                            int(prm.measure_every) ==
                            0)
                    {
                        // Measurements here
                        // measure_counter += 1.0;
                        // cout << "Measuring Si.Sj" << endl;
                        // for (int i = 0; i < ns; i++)
                        // {
                        //     for (int j = 0; j < ns; j++)
                        //     {
                        //         theta1 = real(theta(0, i));
                        //         theta2 = real(theta(0, j));
                        //         SiSj(i, j) += cos(theta1 - theta2);
                        //     }
                        // }
                    }
                }
            }
            // tvsE << temperature << "\t" << totalQE / (1.0 * prm.sweeps) << endl;

            if (temperature > 2.0)
            {
                temperature -= 0.5;
            }
            if (temperature > 0.5 && temperature <= 2.0)
            {
                temperature -= 0.1;
            }
            if (temperature <= 0.5)
            {
                temperature -= 0.01;
            }
        }
    }
    // else
    // {
    // cout << "Performing TCA" << endl;
    // measure_counter = 0;
    // while (temperature >= prm.T_stop - epsilon)
    // {
    //     // cout << "T loop" << endl;
    //     for (int t = 0; t < prm.sweeps; t++)
    //     {
    //         // cout << "t loop" << endl;
    //         total_change = 0;
    //         accepted = 0;
    //         // cout << "Sweep number: " << t + 1 << endl;
    //         for (int i = 0; i < prm.Lx; i++)
    //         {
    //             // cout << "i loop" << endl;
    //             for (int j = 0; j < prm.Ly; j++)
    //             {
    //                 pos = k(i, j);

    //                 theta_old = real(theta(0, pos));
    //                 // cout << "cluster TB 1 made" << endl;
    //                 makeTB_cluster();
    //                 // H_cluster.print();
    //                 // cout << "cluster Hund 1 made" << endl;
    //                 makeHund_cluster(pos);
    //                 // cout << theta(0, 0) << endl;
    //                 // H_cluster.print();
    //                 // Matrix<cd> H_cluster_copy;
    //                 // H_cluster_copy = H_cluster;
    //                 Diagonalize('N', H_cluster, eigs_cluster_);
    //                 // for (int p=0; p < eigs_cluster_.size(); p++){
    //                 //     cout << eigs_cluster_[p] << " ";
    //                 // }
    //                 // cout << endl;
    //                 P_old = lnP_cluster(pos);
    //                 // cout << "P_old = " << P_old << endl;

    //                 theta_new = filter(theta_old + perturb());
    //                 theta(0, pos) = theta_new;
    //                 makeTB_cluster();
    //                 // cout << "cluster TB 2 made" << endl;
    //                 // cout << "theta_old= " << theta_old << " theta_new= "
    //                 //      << theta_new << endl;
    //                 makeHund_cluster(pos);
    //                 // cout << theta(0, 0) << endl;
    //                 // H_cluster.print();
    //                 // for (int p; p < 2 * prm.c_Lx * prm.c_Ly; p++){

    //                 // }
    //                 // cout << "cluster Hund 2 made" << endl;
    //                 Diagonalize('N', H_cluster, eigs_cluster_);
    //                 // for (int p=0; p < eigs_cluster_.size(); p++){
    //                 //     cout << eigs_cluster_[p] << " ";
    //                 // }
    //                 // cout << endl;
    //                 P_new = lnP_cluster(pos);
    //                 // cout << "P_new = " << P_new << endl;
    //                 // assert(false);

    //                 r = rng.random();
    //                 P_ratio = exp(P_new - P_old);

    //                 if (r < P_ratio)
    //                 {
    //                     // cout << "Accepted" << endl;
    //                     accepted += 1;
    //                     theta(0, pos) = theta_new;
    //                 }
    //                 else
    //                 {
    //                     // cout << "Rejected" << endl;
    //                     theta(0, pos) = theta_old;
    //                 }
    //                 total_change += 1;
    //             }
    //         }

    //         cout << "T: " << temperature << ", SN: " << t + 1 << ", AR: "
    //              << (accepted * 1.0) / (total_change * 1.0) << endl;

    //         // Do measurements if thermalized
    //         if (abs(double(temperature - prm.T_stop)) < pow(10, -4))
    //         {
    //             if ((t + 1) >= prm.measure_after &&
    //                 (t - int(prm.measure_after) + 1) %
    //                         int(prm.measure_every) ==
    //                     0)
    //             {
    //                 // Measurements here
    //                 measure_counter += 1.0;
    //                 cout << "Measuring Si.Sj" << endl;
    //                 for (int i = 0; i < ns; i++)
    //                 {
    //                     for (int j = 0; j < ns; j++)
    //                     {
    //                         theta1 = real(theta(0, i));
    //                         theta2 = real(theta(0, j));
    //                         SiSj(i, j) += cos(theta1 - theta2);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     if (temperature > 2.0)
    //     {
    //         temperature -= 0.5;
    //     }
    //     if (temperature > 0.5 && temperature <= 2.0)
    //     {
    //         temperature -= 0.1;
    //     }
    //     if (temperature <= 0.5)
    //     {
    //         temperature -= 0.01;
    //     }
    // }
    // }

    // Printing the final angles
    cout << ":::::::::::Angles:::::::::::::" << endl;
    for (int i = 0; i < prm.Lx; i++)
    {
        for (int j = 0; j < prm.Ly; j++)
        {
            double pos;
            pos = k(i, j);
            cout << fixed << setprecision(3) << setfill('0');
            cout << real(theta(0, pos)) / (2.0 * 3.1415926535) << "\t";
        }
        cout << endl;
    }
    // tvsE.close();

    ofstream Srplot;
    Srplot.open("data/Sr_vs_r.txt");
    // Average out by diving SiSj
    cout << "Times measured = " << measure_counter << endl;
    for (int i = 0; i < ns; i++)
    {
        for (int j = 0; j < ns; j++)
        {
            // theta1 = real(theta(0, i));
            // theta2 = real(theta(0, j));
            SiSj(i, j) = SiSj(i, j) / measure_counter;
            Srplot << i << "\t" << j << "\t"
                   << real(SiSj(i, j)) << "\t" << imag(SiSj(i, j)) << endl;
        }
        Srplot << endl;
    }
    Srplot.close();

    Sq.resize(prm.Lx, prm.Ly);

    // cout << ":::::::::::S(q):::::::::::" << endl;

    double qx, qy, ix, iy, jx, jy;

    ofstream sqplot;
    sqplot.open("data/Sq.txt");

    for (int qnx = 0; qnx < prm.Lx; qnx++)
    {
        for (int qny = 0; qny < prm.Ly; qny++)
        {
            for (int i = 0; i < ns; i++)
            {
                for (int j = 0; j < ns; j++)
                {
                    qx = 2.0 * 3.1415926535 * qnx / double(prm.Lx);
                    qy = 2.0 * 3.1415926535 * qny / double(prm.Ly);
                    ix = kinv(i).first;
                    iy = kinv(i).second;
                    jx = kinv(j).first;
                    jy = kinv(j).second;
                    // cout << ix << " "
                    //      << iy << " " << jx << " " << jy << endl;
                    Sq(qnx, qny) +=
                            (1.0 / double(ns)) *
                            exp(imagi * qx * (ix - jx)) *
                            exp(imagi * qy * (iy - jy)) * SiSj(i, j);
                }
            }
            // cout << qx << "\t" << qy << "\t" << real(Sq(qnx, qny))
            //      << "\t" << imag(Sq(qnx, qny)) << endl;
            sqplot << qx << "\t" << qy << "\t" << real(Sq(qnx, qny))
                   << "\t" << imag(Sq(qnx, qny)) << endl;
        }
        // cout << endl;
        sqplot << endl;
    }

    sqplot.close();
    return 0;
}

// Diagonalize('V', H, eigs_);
// cout << eigs_.size() << endl;
// cout << "Finished diagonalization." << endl;
// ofstream testdos;
// testdos.open("testdos.txt");
// for (int i = 0; i < eigs_.size(); i++)
// {
//     testdos << eigs_[i] << endl;
// }
// testdos.close();
// cout << getmu() << endl;
