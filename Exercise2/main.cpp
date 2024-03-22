#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


double errore_rel_PALU(MatrixXd& A, VectorXd& b, VectorXd& xSOL, VectorXd& x){
    x = A.lu().solve(b);
    double err_rel  = (x-xSOL).norm()/xSOL.norm(); //la funzione norm calcola la norma euclidea
    return err_rel;
}

double errore_rel_QR(MatrixXd &A, VectorXd &b, VectorXd &xSOL, VectorXd& x){
    x = A.householderQr().solve(b);
    double err_rel  = (x-xSOL).norm()/xSOL.norm(); //la funzione norm calcola la norma euclidea
    return err_rel;
}


int main()
{
    unsigned int n = 2;
    //inizializzo i tre sistemi lineari
    MatrixXd A1(n,n);
    A1 << 5.547001962252291e-01,  -3.770900990025203e-02,
        8.320502943378437e-01, -9.992887623566787e-01;
    VectorXd b1(n,1);
    b1 << -5.169911863249772e-01,
        1.672384680188350e-01;

    MatrixXd A2(n,n);
    A2 << 5.547001962252291e-01, -5.540607316466765e-01,
        8.320502943378437e-01, -8.324762492991313e-01;

    VectorXd b2(n,1);
    b2 << -6.394645785530173e-04,
        4.259549612877223e-04;

    MatrixXd A3(n,n);
    A3 << 5.547001962252291e-01, -5.547001955851905e-01,
        8.320502943378437e-01, -8.320502947645361e-01;
    VectorXd b3(n,1);
    b3 << -6.400391328043042e-10,
        4.266924591433963e-10;

    VectorXd xSOL =- VectorXd::Ones(n,1);

    VectorXd x1P = VectorXd::Zero(2,1);
    VectorXd x1QR = VectorXd::Zero(2,1);
    double err_rel_1_PALU = errore_rel_PALU(A1, b1, xSOL, x1P);
    double err_rel_1_QR = errore_rel_QR(A1, b1, xSOL, x1QR);
    cout << "L'errore relativo associato alla risoluzione con la fattorizzazione PA = LU del sistema 1 viene: " << err_rel_1_PALU << endl;
    cout << "L'errore relativo associato alla risoluzione con la fattorizzazione A = QR del sistema 1 viene: " << err_rel_1_QR << endl;

    VectorXd x2P = VectorXd::Zero(2,1);
    VectorXd x2QR = VectorXd::Zero(2,1);
    double err_rel_2_PALU = errore_rel_PALU(A2, b2, xSOL, x2P);
    double err_rel_2_QR = errore_rel_QR(A2, b2, xSOL, x2QR);
    cout << "L'errore relativo associato alla risoluzione con la fattorizzazione PA = LU del sistema 2 viene: " << err_rel_2_PALU << endl;
    cout << "L'errore relativo associato alla risoluzione con la fattorizzazione A = QR del sistema 2 viene: " << err_rel_2_QR << endl;

    VectorXd x3P = VectorXd::Zero(2,1);
    VectorXd x3QR = VectorXd::Zero(2,1);
    double err_rel_3_PALU = errore_rel_PALU(A3, b3, xSOL, x3P);
    double err_rel_3_QR = errore_rel_QR(A3, b3, xSOL, x3QR);
    cout << "L'errore relativo associato alla risoluzione con la fattorizzazione PA = LU del sistema 3 viene: " << err_rel_3_PALU << endl;
    cout << "L'errore relativo associato alla risoluzione con la fattorizzazione A = QR del sistema 3 viene: " << err_rel_3_QR << endl;


  return 0;
}
