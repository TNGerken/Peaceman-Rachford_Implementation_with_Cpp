#include <fstream>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <string>
using namespace Eigen;

const double nu=1.0;
const double dt=0.0005;
const int Mx=20;
const int My=20;
const double N=1;

const double dx=1.0 / Mx;
const double dy=1.0 / My;

const double r_x=nu*dt/(dx*dx);
const double r_y=nu*dt/(dy*dy);

const int K = static_cast<int>(N/dt);

MatrixXd A_x=MatrixXd::Zero(Mx-1,Mx-1);
MatrixXd B_x=MatrixXd::Zero(Mx-1,1);
MatrixXd A_y=MatrixXd::Zero(My-1,My-1);
MatrixXd B_y=MatrixXd::Zero(My-1,1);

int K2 =static_cast<int>(2*N/dt+1);
std::vector<MatrixXd> u(K2,MatrixXd::Zero(My+1,Mx+1));


double phi_0(double x, double y){
    return sin(M_PI*x)*sin(2*M_PI*y);
}

Eigen::VectorXd THOM(const Eigen::MatrixXd& A, int m, const Eigen::MatrixXd& phi){
    using namespace Eigen;
    VectorXd alpha =VectorXd::Zero(m);
    VectorXd s= VectorXd::Zero(m);
    VectorXd a(m-1), b(m), c(m-1);

    // Extract diagonals properly
    for (int i = 0; i < m-1; ++i) {
        a(i)= A(i+1,i); // Lower diagonal
        c(i)= A(i,i+1); // Upper diagonal
    }

    for (int i =0; i <m; ++i){
        b(i) =A(i,i); // Main diagonal
    }

    // Flatten phi to get d
    VectorXd d= Eigen::Map<const VectorXd>(phi.data(),phi.size());

    VectorXd phi_2= VectorXd::Zero(m);

    // Forward sweep
    alpha(0)= b(0);
    s(0)= d(0);

    for (int i=1; i< m; ++i){
        alpha(i)=b(i)-a(i-1)*c(i-1)/alpha(i-1);
        s(i) = d(i) - a(i-1)*s(i-1)/alpha(i-1);
    }

    // Backward substitution
    phi_2(m-1)=s(m-1)/alpha(m-1);

    for (int i=m-2;i>=0; --i) {
        phi_2(i)=(s(i)-c(i)*phi_2(i+1))/alpha(i);
    }
    return phi_2;
}

//Exporting to excel
void write_csv_slice(const std::vector<std::vector<std::vector<double>>>& u, double t, double dt, const std::string& filename) {
int timestep = static_cast<int>(t / dt);

std::ofstream file(filename);
if (!file) {
    std::cerr << "Could not open file: " << filename << "\n";
    return;
}

// Optional header: x0, x1, x2, ...
for (int j=0; j < u[0].size(); ++j) {
    file << "x" << j;
    if (j < u[0].size() - 1) file << ",";
}
file << "\n";

for (int i=0; i<u.size(); ++i) {
    for (int j=0; j<u[i].size(); ++j) {
        file << std::setprecision(10) << u[i][j][timestep];
        if (j < u[i].size() - 1) file << ",";
    }
    file << "\n";
}

file.close();
std::cout << "Wrote: " << filename << "\n";
}

int main() {
    //Defining A_x and A_y
    for (int i=0; i<Mx-1;++i) {
        for (int j=0; j<Mx-1; ++j){
            if (i==j){
                A_x(i,j)=r_x+1;
                A_y(i,j) = r_y+1;
            }
            else if (i==j-1){
                A_x(i,j)=-r_x/2;
                A_y(i,j)=-r_y/2;
            }
            else if (i==j+1){
                A_x(i,j)=-r_x/2;
                A_y(i,j)=-r_y/2;
            }
        }
    }

    // Initial conditions
    for (int i=1; i<My; ++i){
        for (int j=1; j<Mx; ++j){
            //u[0](i,j)=phi_0(dx*j,dy*i);
            u[0](i,j) = phi_0(dx*i,dy*j);
        }
    }

    for (int n=1; n<= static_cast<int>(N/dt); ++n){
        // X-direction sweep
        for (int i=1; i<My; ++i){
            for (int j=1; j<Mx; ++j){
                B_x(j-1, 0)=(r_y/2.0)*u[2*n-2](j,i+1)+(1-r_y)*u[2*n-2](j,i)+(r_y/2.0)*u[2*n-2](j,i-1);
            }
            Eigen::VectorXd result = THOM(A_x, Mx-1, B_x);
            for (int j=1; j<Mx; ++j){
                u[2*n-1](i,j)=result(j-1);
            }
        }

        //Enforcing boundary conditions
        for (int i=0; i<=My; ++i){
            u[2*n-1](i,0)=0;// Left boundary
            u[2*n-1](i,Mx)=0;// Right boundary
        }
        for (int j=0; j<=Mx; ++j){
            u[2*n-1](0,j)=0;// Bottom boundary
            u[2*n-1](My,j)=0;// Top boundary
        }

        //Y-direction sweep
        for (int j=1; j<Mx; ++j){
            for (int i=1; i<My; ++i){
                B_y(i - 1, 0) = (r_x / 2.0) * u[2 * n - 1](j+1,i) + (1 - r_x) * u[2 * n - 1](j,i) + (r_x / 2.0) * u[2 * n - 1](j-1,i);
            }
            Eigen::VectorXd result = THOM(A_y, My-1, B_y);
            for (int i=1; i<My; ++i){
                u[2*n](i,j)=result(i-1);
            }
        }

        //Enforcing boundary conditions
        for (int i=0; i<=My; ++i){
            u[2*n](i,0)=0;// Left boundary
            u[2*n](i,Mx)=0;// Right boundary
        }
        for (int j=0; j<=Mx; ++j){
            u[2*n](0,j)=0;// Bottom boundary
            u[2*n](My,j)=0;// Top boundary
        }
    }

    // End of numerical portion

    //Exporting to Excel
    int steps = static_cast<int>(2*N/dt+1);

    std::vector<std::vector<std::vector<double>>> u_export(My+1, std::vector<std::vector<double>>(Mx+1, std::vector<double>(steps)));
    for (int t=0; t<steps; ++t) {
        for (int i=0;i<=My;++i) {
            for (int j = 0;j<= Mx;++j) {
                u_export[i][j][t]=u[t](i,j);
            }
        }
    }

    write_csv_slice(u_export, 2*0.06, dt, "T_0.06.csv");
    write_csv_slice(u_export, 2*0.1, dt, "T_0.1.csv");
    write_csv_slice(u_export, 2*0.9, dt, "T_0.9.csv");

    return 0;
}
