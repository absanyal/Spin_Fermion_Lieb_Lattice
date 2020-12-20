#ifndef DIAG_H
#define DIAG_H
#include <vector>
#include <complex>
#include "Matrix.h"

// using namespace std;
using std::complex;
using std::cout;
using std::endl;
using std::pair;
using std::vector;

extern "C" void zheev_(char *, char *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, double *, int *);

void Diagonalize(char option, Matrix<cd> &Ham_, vector<double> &eigs_)
{

    char jobz = option;
    char uplo = 'U'; //WHY ONLY 'L' WORKS?
    int n = Ham_.n_row();
    int lda = Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3 * n - 2);
    int info;
    int lwork = -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(), eigs_.end(), 0);
    // query:
    zheev_(&jobz, &uplo, &n, &(Ham_(0, 0)), &lda, &(eigs_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz, &uplo, &n, &(Ham_(0, 0)), &lda, &(eigs_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    if (info != 0)
    {
        std::cerr << "info=" << info << "\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // cout << "" << endl;
    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}
}
// //*******************************************************************************//
#endif