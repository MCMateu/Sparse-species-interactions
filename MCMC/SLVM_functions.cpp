#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <random>
#include <tuple>
#include <vector>
#include <omp.h>


using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// FUNCTION TO INTEGRATE THE SLVM USING EULER-MARUYAMA
MatrixXd SLVM(
    const MatrixXd& M,
    const MatrixXd& b,
    double tau,
    int nSim,
    const VectorXd& timeList,
    double abundanceStart,
    int S
) {
    // Initialize random number generator
    random_device rd;
    mt19937 gen(rd()+omp_get_thread_num());
    normal_distribution<> normalDist(0, 1);

    // Initialize initial state
    VectorXd x0 = VectorXd::Constant(S, abundanceStart);

    // Compute time step
    double dt = (timeList(timeList.size() - 1) - timeList(0)) / (timeList.size() - 1);
    double sqrt_dt = sqrt(dt);

    // Initialize abundance matrix
    MatrixXd abundanceMat=MatrixXd::Zero(timeList.size(), S * nSim);
    abundanceMat.setZero();

    VectorXd noise;

    //#pragma omp parallel for
    for (int sim = 0; sim < nSim; ++sim) {
        abundanceMat.row(0).segment(sim * S, S) = x0;
        
        // Time loop
        for (int t = 1; t < timeList.size(); ++t) {

            // Dynamics
            VectorXd x = abundanceMat.row(t - 1).segment(sim * S, S);
            VectorXd z = generateGaussianNoise(0, 1, S);
            noise = b * z;
            VectorXd xt = x + (x * dt) / tau + (dt / tau) * ( x.cwiseProduct(M * x) ) + sqrt_dt * noise.cwiseProduct(x);
            abundanceMat.row(t).segment(sim * S, S) = xt;
          
        }
    }

    int T = abundanceMat.rows();

    int skip_transitory=2/dt; //number of time-points to skip

    int N=0.05/dt;

    int transformedRows = ( (T - skip_transitory) * nSim ) / N;

    MatrixXd transformedTable=MatrixXd::Zero(transformedRows, S);
    
    int transformedRow = 0;

    //#pragma omp parallel for
    for (int sim = 0; sim < nSim; ++sim) {
        
        for (int t = skip_transitory-1; t < T-N; t+=N) {

            //cout<<"t"<<t<<" transformedRow "<<transformedRow<<"/"<<transformedRows<<endl;

            transformedTable.row(transformedRow) = abundanceMat.row(t).segment(sim * S, S);

            transformedRow++;

        }
    }

    return transformedTable;

}



//Function to produce random Lotka-Volterra Matrix
MatrixXd InteractionMatrix(int S, double C, double muInt, double sigmaInt, const vector<double>& K) {
    /*
    Interaction matrix with K in the diagonal and Gaussian off-diagonal fluctuations;

    INPUT
    - S: species number;
    - C: connectance;
    - muInt: Gaussian mean of the off-diagonal elements;
    - sigmaInt: Gaussian variance of the off-diagonal elements;
    - K: vector of maximum carrying capacities 

    OUTPUT
    - M: SxS matrix where the i,j element provides the interaction weight between species i and j
    */

    // Create a random number generator
    random_device rd;
    mt19937 gen(rd()+omp_get_thread_num());
    normal_distribution<double> dist(muInt, sigmaInt);
    uniform_real_distribution<double> connectance_dist(0.0, 1.0);

    // Create the interaction matrix A
    MatrixXd A = MatrixXd::Zero(S, S);
    for (int i = 0; i < S; i++) {
        for (int j = 0; j < S; j++) {
            A(i, j) = dist(gen);
        }
    }

    // Apply connectance and diagonal adjustments to A
    for (int i = 0; i < S; i++) {
        for (int j = 0; j < S; j++) {
            double random_value = connectance_dist(gen);
            if (random_value > C) {
                A(i, j) = 0.0;
            }
        }
        A(i, i) = 0.0;
    }

    // Create the interaction matrix M
    MatrixXd M = MatrixXd::Zero(S, S);
    for (int i = 0; i < S; i++) {
        for (int j = 0; j < S; j++) {
            M(i, j) = A(i, j) - (i == j ? 1.0 / K[i] : 0.0);
        }
    }

    return M;
}


// Function to create the correlation noise matrix
MatrixXd EnvironmentalFilterMatrix(int S) {
    // Generate random orthogonal matrix
    MatrixXd s = MatrixXd::Random(S, S);
    Eigen::HouseholderQR<MatrixXd> qr(s);
    MatrixXd orthogonalMatrix = qr.householderQ();

    // Generate random diagonal matrix
    random_device rd;
    mt19937 gen(rd()+omp_get_thread_num());
    uniform_real_distribution<> dis(0.0, 1.0);
    VectorXd diagonal = VectorXd::NullaryExpr(S, [&]() { return dis(gen); });
    MatrixXd diagonalMatrix = diagonal.asDiagonal();

    // Compute correlation noise matrix
    MatrixXd D = orthogonalMatrix * diagonalMatrix * orthogonalMatrix.transpose();

    return D;
}

MatrixXd EquivalentEnvFiltMatrix(const MatrixXd& D) {

    /*
    Function to generate n-dimensional Gaussian Random numbers with
    correlation matrix D
    
    The result is a matrix, b, with the coefficients for a linear 
    combination of non-correlated noise that produce the same effect
    of the correlated n-dimensional gaussian
    */


    int S = D.rows();
    MatrixXd b = MatrixXd::Zero(S, S);

    for (int j = 0; j < S; j++) {
        for (int i = 0; i < j + 1; i++) {
            double z = D(j, i);
            for (int k = 0; k < i; k++) {
                z = z - b(j, k) * b(i, k);
            }
            if (i < j) {
                b(j, i) = z / b(i, i);
            } else {
                b(j, j) = sqrt(z);
            }
        }
    }

    return b;
}

/*
tuple<double, double, double,double> mat_statistics(const MatrixXd& M){

    const VectorXd interactionList = Eigen::Map<const VectorXd>(M.data(), M.size());

    int nonzeroInteractions = (interactionList.array() != 0).count();
    int totalInteractions = interactionList.size();

    double C = (1.0*nonzeroInteractions) / totalInteractions;
    double C_eff = (interactionList.array() >= (interactionList.array().maxCoeff()*0.01)).count() / (totalInteractions*1.0);

    double mu = interactionList.mean();

    VectorXd center = (interactionList.array() - mu);
    double sigma = center.norm() / sqrt(interactionList.size());

    return make_tuple(C,C_eff,mu,sigma);

}
*/

tuple<double, double, double,double> mat_statistics(const MatrixXd& M, const vector<pair<int, int>>& upTrList, const vector<pair<int, int>>& lowTrList){

    VectorXd interactionList(M.cols()*(M.cols()-1));

    //#pragma omp parallel for
    for (int i = 0; i < upTrList.size(); ++i) {
        int row = upTrList[i].first;
        int col = upTrList[i].second;
        interactionList(i) = M(row, col);
    }

    int k=upTrList.size();

   //#pragma omp parallel for
   for (int i = 0; i < lowTrList.size(); ++i) {
        int row = lowTrList[i].first;
        int col = lowTrList[i].second;
        interactionList(k) = M(row, col);
        k++;
    }

    double C;
    double C_eff;
    double sigma;
    double mu;

    //CONECTIVITY

    int nonzeroInteractions = (interactionList.array() != 0).count();
    int totalInteractions = interactionList.size();

    cout<<totalInteractions<<"   "<<(1.0*nonzeroInteractions)<<endl;

    if(nonzeroInteractions !=0){

	    C = (1.0*nonzeroInteractions) / totalInteractions;
	    
	    //EFFECTIVE CONECTANCE
            VectorXd AbsinteractionList = interactionList.cwiseAbs();

            double C_eff = (AbsinteractionList.array() >= (AbsinteractionList.array().maxCoeff()*0.01)).count() / (totalInteractions*1.0);

	    //MEAN AND STANDARD DEVIATION

	    VectorXd nonZeroList (nonzeroInteractions);
	    k=0;
	    for(int i=0;i<totalInteractions;i++){

            if(interactionList(i)!=0){ 
                
                nonZeroList(k)=interactionList(i);
                k++;
		    
		    }

	    }

	    mu = nonZeroList.mean();

	    VectorXd center = (nonZeroList.array() - mu);
	    sigma = center.norm() / sqrt(nonzeroInteractions);
     }
     else{
     	mu=0;
     	sigma=0;
     	C=0;
     	C_eff=0;
     }

    return make_tuple(C,C_eff,mu,sigma);

}

