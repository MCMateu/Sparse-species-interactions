//Libraries

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <tuple>
#include <filesystem>
#include <Eigen/Dense>
#include "functions.cpp" 
#include "SLVM_functions.cpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
namespace fs = std::filesystem;

//*****************************//
//****  GLOBAL CONSTANTS   ****//
//*****************************//

//BIOME DATA
const string occupancy = "0.95";
const string biomeFold = "/home/jose/small-perturbation/EmpiricalData&PCA/PearsMat";
const string biomeName = "Lake";
const string biomeDataDir = biomeFold + "/" + biomeName + "/" + "o_" + occupancy + "/" + "PearsMat.csv";

//BINS FOR ABUNDANCE CORRELATION DISTRIBUTION
const double binCorrMin = -1.0;
const double binCorrMax = 1.0;
const int binCorrNum = 100;
const VectorXd binCorrList = linspace(binCorrMin, binCorrMax, binCorrNum);

//EMPIRICAL ABUNDANCE CORRELATION DISTRIBUTION
tuple<VectorXd , double, MatrixXd> result = EmpiricalCorrelationDistribution(biomeDataDir, binCorrList,true);

//READ EACH ELEMENT OF THE TUPLE "result"
VectorXd  freqList = get<0>(result);
int S = get<1>(result);
MatrixXd correlationMatrix = get<2>(result);


//MEAN ABUNDANCE DISTRIBUTION (MAD) 
const double binMadMin = -2.5;
const double binMadMax = 2.5;
const double binMadNum = 30;
const VectorXd  binMadList = linspace(binMadMin,binMadMax,binMadNum);
const VectorXd  freqMadIdeal = gaussianDistribution(binMadList,0,1);


//********************************************//
//*******    SIMULATION PARAMETERS     *******//
//********************************************//

const double tStart = 0;
const double tEnd = 7;
const double dt = 0.05;
const VectorXd  timeList = arrange(tStart, tEnd, dt);

//Initial Guess
const double tau = 0.1;
const double muNoise = 0;
const double sigmaNoise = 0.1;

const double muInt = 0;
const double sigmaInt = 0;

const int nSim = 5;
const double C=0;
const double abundanceStart=0.1;

int main (void){


    cout<<"Species number: "<<S<<endl;

    //********************************************//
    //*******           FOLDERS            *******//
    //********************************************//

    //Main seed folder
    string OutFold = "/home/jose/MEGA/CPP_MCMC/diag_seeds";
    if (!fs::is_directory(OutFold)){

        fs::create_directory(OutFold);

    }

    //biome
    string biomeOutFold = OutFold +"/" + biomeName; 
    if (!fs::is_directory(biomeOutFold)){

        fs::create_directory(biomeOutFold);

    }

    //occupancy
    string occupancyOutFold = biomeOutFold  + "/" + "o_"+ occupancy;
    if (!fs::is_directory(occupancyOutFold)){

        fs::create_directory(occupancyOutFold);

    }


    //********************************************//
    //*******        DIAGONAL SEEDS        *******//
    //********************************************//

    //Diagonal seed list
    VectorXd muLogNormal_List(1);
             muLogNormal_List << 0.1;
         
    VectorXd sigmaLogNormal_List(1); 
             sigmaLogNormal_List << 0.5;

    vector<double> K(S);

    for(int i=0; i < muLogNormal_List.size(); i++){

        for(int j=0; j < sigmaLogNormal_List.size() ; j++){

            double muLN = muLogNormal_List(i);
            double sigmaLN = sigmaLogNormal_List(j);

            string expFold = occupancyOutFold  + "/" + "muLN_"+  as_string(muLN)+"_sigmaLN_"+ as_string(sigmaLN);
            if (!fs::is_directory(expFold)){

                fs::create_directory(expFold);

            }

            // Set up random number generator
            random_device rd;
            mt19937 generator(rd()+omp_get_thread_num());
            lognormal_distribution<double> ln_distribution(muLN, sigmaLN);

            
            for(int k=0; k < 10; k++){
                
                for(int i=0;i<S;i++)
                {
                     K[i] = ln_distribution(generator);
                }


                //Correlation Matrix
                //MatrixXd W=EnvironmentalFilterMatrix(S);

                //Equivalent correlation Matrix
                //MatrixXd b=EquivalentEnvFiltMatrix(W); 

                VectorXd sigmaNoiseVector(S);
                sigmaNoiseVector.fill(sigmaNoise);
                MatrixXd b = sigmaNoiseVector.asDiagonal();

                writeCSV(b,"/home/jose/dif.csv");

                // Variable declarations
                double leadingEigenvalue = 0.0;
                int extinctionsNum = 0;
                int anomalies = 0;
                int infinity = 0;
                MatrixXd M;

                do {

                    // Interaction Network
                    M = InteractionMatrix(S, C, muInt, sigmaInt, K);

                    // Eigenvalue
                    leadingEigenvalue = LeadingEigenvalue(M);

                    // Extinctions
                    extinctionsNum = ExtinctionsNumber(M);

                    // Dynamics
                    MatrixXd abundanceMat = SLVM(M, b, tau, nSim, timeList, abundanceStart, S);

                    // Anomalies
                    //Eigen::VectorXd flattenedVector = Eigen::Map<Eigen::VectorXd>(matrix.data(), matrix.size());
                    VectorXd abundancesFlat = Eigen::Map<VectorXd>(abundanceMat.data(), abundanceMat.size());
                    anomalies = (abundancesFlat.array().isNaN()).count();

                    // Infinity
                    infinity = (abundancesFlat.array().isInf()).count();

                    // Print the values
                    cout << infinity << " " << anomalies << " " << leadingEigenvalue << " " << extinctionsNum << endl;

                } while (leadingEigenvalue >= 0 || extinctionsNum > 0 || anomalies > 0 || infinity > 0);

                //Save matrix
                writeCSV(M, expFold+"/"+"diag_"+as_string(k)+".csv");
                
            }

        }

    }

    return 0;

}

