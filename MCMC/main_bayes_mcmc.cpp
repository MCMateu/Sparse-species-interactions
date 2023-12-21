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

using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
namespace fs = std::filesystem;

//*****************************//
//****  GLOBAL CONSTANTS   ****//
//*****************************//

//BIOME DATA
const string occupancy = "0.5";
const string biomeFold = "/home/jose/small-perturbation/EmpiricalData&PCA/PearsMat";
const string biomeName = "Seawater";
const string biomeDataDir = biomeFold + "/" + biomeName + "/" + "o_" + occupancy + "/" + "PearsMat.csv";

//BINS FOR ABUNDANCE CORRELATION DISTRIBUTION
const double binCorrMin = -1.0;
const double binCorrMax = 1.0;
const int binCorrNum = 100;
const VectorXd binCorrList = linspace(binCorrMin, binCorrMax, binCorrNum);

//EMPIRICAL ABUNDANCE CORRELATION DISTRIBUTION
tuple<VectorXd , double, MatrixXd> result = EmpiricalCorrelationDistribution(biomeDataDir, binCorrList,true);

//READ EACH ELEMENT OF THE TUPLE "result"
VectorXd  biomeFreqList = get<0>(result);
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

double muInt = 0;
double sigmaInt = 0;

const int nSim = 5;
double C=0;
double C_eff;
const double abundanceStart=0.1;

//********************************************//
//*******        MCMC PARAMETERS       *******//
//********************************************//

const int numIterations = 100000;
const double sigmaCorr = 2;
const double sigmaMad = 0.3;
const double sigmaTaylor = 0.1;
const int numItSave = 25;
const double epsilonScale = 0.03;

vector<pair<int, int>> upTrList = triu_indices(S, S, 1); //upper triangle indexes
vector<pair<int, int>> lowTrList = tril_indices(S, S, 1); //lower triangle indexes

int main (void){


    // RANDOM NUMBER GENERATORS

        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<double> uniform01_R(0.0, 1.0);

    /************************************************************************/

    /***********       READ SEED DATA AND CREAT OUTPUT FOLDER     ***********/   

    /************************************************************************/


    //READ DIAGONAL PARAMETERS FOR 

    //Choose diagonal seed
    int diag=0;

    //Choose seed num

    int seedNum=1;

    string SeedFold = "/home/jose/MEGA/CPP_MCMC/diag_seeds/"+biomeName+"/o_"+occupancy;

    //List of ensembles of diagonals
    Eigen::Matrix<string, Eigen::Dynamic, 1> diag_parametersFolder = listDir(SeedFold); 

    string folder_name = diag_parametersFolder(diag);

    VectorXd diagParams = get_nums(folder_name);
    
    double muLN=diagParams[0]; 
    double sigmaLN=diagParams[1];

    // CREATE OUTPUT FOLDERS

    // main
    const string mainOutFold = "/home/jose/MEGA/CPP_MCMC/data/";

    // biome
    const string biomeOutFold = mainOutFold + biomeName;

    if (!fs::is_directory(biomeOutFold)){
        fs::create_directory(biomeOutFold);
    }

    // occupancy
    const string occupancyOutFold = biomeOutFold + "/" + "o_" + occupancy;
    if (!fs::is_directory(occupancyOutFold)){
        fs::create_directory(occupancyOutFold);
    }

    // parameters
    const string expOutFold = occupancyOutFold + "/seed_" + as_string(seedNum) + "_binsCorr_" + as_string(binCorrNum) +
                            "_binsMAD_" + as_string(binMadNum) + "_nSim_" + as_string(nSim) + "_sigmaCorr_" +
                            as_string(sigmaCorr) + "_sigmaMad_" + as_string(sigmaMad) + "_sigmaTaylor_"+as_string(sigmaTaylor)+"_sigmaNoise_" +
                            as_string(sigmaNoise);

    if (!fs::is_directory(expOutFold)){
        fs::create_directory(expOutFold);
    }

    
    //const string expOutFold="/home/jose/single-chain";

    /**************************************/

    /*   FIRST RUN OR READ LAST MATRIX    */

    /**************************************/

    //write first run if it the is the first time you run the code
    //write something else to start the method using the last matrix you generated
    string flag = "first_run";

    const string Seed_File_path=SeedFold+"/muLN_"+as_string(muLN)+"_sigmaLN_"+as_string(sigmaLN)+"/diag_0.csv";

    MatrixXd initialGuessMat;
    int lastIter=0;

    if(flag == "first_run"){

        const string Seed_File_path=SeedFold+"/muLN_"+as_string(muLN)+"_sigmaLN_"+as_string(sigmaLN)+"/diag_0.csv";

        initialGuessMat=readCSV(Seed_File_path);

        lastIter=0;

    } else{

        //RESUME THE CHAIN FROM THE LAST ITERATION SAVED

        Eigen::Matrix<std::string, Eigen::Dynamic, 1> previousIterFileNames = listFiles(expOutFold);  

        vector<string> previousIterMat;
        for (int i = 0; i < previousIterFileNames.size(); ++i) {
            const std::string& fileName = previousIterFileNames(i);
            if (fileName.find("mat_") == 0) {
                previousIterMat.push_back(fileName);
            }
        }
     
        std::vector<int> previousIter;
        for (const auto& matFileName : previousIterMat) {
            std::string iterationStr = matFileName.substr(4);
            int iteration = std::stoi(iterationStr);
            previousIter.push_back(iteration);
        }
   
        std::sort(previousIter.begin(), previousIter.end(), std::greater<int>());
        lastIter = previousIter[0];

        cout << "Last Iteration: " << lastIter << endl;
        cout<<expOutFold+"/"+"mat_"+as_string(lastIter)+".csv"<<endl;
        initialGuessMat=readCSV(expOutFold+"/"+"mat_"+as_string(lastIter)+".csv");

    }

    tuple<double, double, double, double> result=mat_statistics(initialGuessMat,upTrList,lowTrList);

    /************************************************************************/

    /******                   METROPOLIS-HASTINGHS                     ******/

    /************************************************************************/

    // NOISE MATRIX

    //MatrixXd D=EnvironmentalFilterMatrix(S);
    //MatrixXd b=EquivalentEnvFiltMatrix(D);

    //FOR THE SAKE OF SIMPLICITY KEEP THE ENVIRONMENTAL NOISE MATRIX DIAGONAL
    VectorXd sigmaNoiseVector(S);
    sigmaNoiseVector.fill(0.1);
    MatrixXd b = sigmaNoiseVector.asDiagonal();

    // RUN FIRST DYNAMICS

    int anomalies = 1; //number of abundance NaN 
    int infinity = 1; //number of abundance divergences

    MatrixXd abundance_table;
    VectorXd abundancesFlat;

    do{

        // Call the StochasticLotkaVolterraEnviromentalFilter function
        abundance_table = SLVM(initialGuessMat, b, tau, nSim, timeList, abundanceStart, S);

        abundancesFlat = Eigen::Map<VectorXd>(abundance_table.data(), abundance_table.size());

        // Anomalies
        anomalies = (abundancesFlat.array().isNaN()).count();

        // Infinity
        infinity = (abundancesFlat.array().isInf()).count();

        cout<<anomalies<<"   "<<infinity<<endl;

    }while(anomalies>1 && infinity>1);
   
    writeCSV(abundance_table,expOutFold+"/abundance_mat.csv");

    //***********************//

    //     ITERATION 0      //

    //**********************//
    
        // CORRELATION

        VectorXd freqCorrList =  SyntheticCorrelationDistribution(abundance_table, S, upTrList, binCorrList);

        double correlationDistance = EuclideanDistance(freqCorrList,biomeFreqList);

        cout<<"First correlation distance:  "<<correlationDistance<<endl;

        MatrixXd matrix = joinVectors(binCorrList.head(binCorrList.size() - 1), biomeFreqList);

        writeCSV(matrix,expOutFold+"/Empirical.csv"); 

        matrix = joinVectors(binCorrList.head(binCorrList.size() - 1), freqCorrList);

        writeCSV(matrix,expOutFold+"/Iter0.csv"); 

        // MEAN ABUNDANCE DISTRIBUTION (MAD)

        VectorXd Mean_Abundance_Log = abundance_table.rowwise().mean().array().log();    

        VectorXd centered =  Mean_Abundance_Log.array() - Mean_Abundance_Log.mean();
        VectorXd zScoreLogMean = centered.normalized()*sqrt(abundance_table.cols());
        VectorXd initialFreqMadList=histogram(zScoreLogMean,binMadList);

        VectorXd FreqMadMin = initialFreqMadList;

        double madDistance = EuclideanDistance(initialFreqMadList,freqMadIdeal);

        // TAYLOR LAW
        double Slope = Taylor_Slope(abundance_table);
        double SlopesDiff = abs(2-Slope);

        cout<<madDistance<<endl;
        cout<<madDistance/(2*pow(sigmaMad,2))<<"   "<<correlationDistance/(2*pow(sigmaCorr,2))<<"   "<<SlopesDiff/(2*pow(sigmaTaylor,2))<<endl;


        // LIKELIHOOD

        double Q0 = correlationDistance/(2*pow(sigmaCorr,2))+ madDistance/(2*pow(sigmaMad,2))+SlopesDiff/(2*pow(sigmaTaylor,2));

        double Q1;  

                        /************************************************************/

                        /*************             MAIN LOOP            *************/

                        /************************************************************/

        // declare some variables 

        int iNorm,acceptedChanges=0;

        double acceptedRatio;

        double leadingEigenvalue,H,q,epsilon;

        MatrixXd abundance_tableMin,covMat;

        std::string matDir,CovMatDir;

        ofstream parametersDirFile(expOutFold+"/parameters_"+as_string(lastIter)+".csv");
        
        // initialize

        bool save;

        VectorXd freqCorrMin,FreqMadList;

        double Qmin=Q0; 
        MatrixXd M0 = initialGuessMat; // Initial guess mat
        MatrixXd Mmin = initialGuessMat; // Minimum distance matrix inizialization
        MatrixXd M1; // Perturbed matrix


        // Start the timer
        double startTime = omp_get_wtime();

    cout<<"Chain starts!"<<endl;

    for(int i=lastIter;i<numIterations+lastIter+1;i++){
         
        //save=false;

        // ACCEPTANCE 

            iNorm = i-lastIter;
            acceptedRatio = acceptedChanges/((iNorm+1)*1.0);

        // MODIFY NETWORK
        epsilon=epsilonScale*Q0;
        M1 = MatrixElementPerturbation( M0, S, epsilon );

        // DYNAMICS
        
        abundance_table = SLVM(M1, b, tau, nSim, timeList, abundanceStart, S);

        
        // STABILITY AND FEASIBILITY

            // STABILITY: LEADING EIGEN-VALUE      
            leadingEigenvalue=LeadingEigenvalue(M1);
                    
            if(leadingEigenvalue>=0){

                cout<<i<<" unsatable "<<endl;

                continue; //skip current iteration

            }

            //FEASIBILITY:  
            abundancesFlat = Eigen::Map<VectorXd>(abundance_table.data(), abundance_table.size());
            
                int extinctionsNumber = (abundancesFlat.array() <= 1e-6).count();

                if (extinctionsNumber > 0) {
                    cout << "Iteration " << i << " produced unfeasible network" << endl;
                    continue;
                }

            //ANOMALIES 

                // Anomalies
                anomalies = (abundancesFlat.array().isNaN()).count();

                // Infinity
                infinity = (abundancesFlat.array().isInf()).count();

                if(anomalies>0 || infinity>0){

                    cout<<i<<" anomalies or unfeasible "<<endl;

                    continue; //skip current iteration

                }

        // CORRELATION DISTRIBUTION
        freqCorrList = SyntheticCorrelationDistribution(abundance_table, S, upTrList, binCorrList);
        
            // DISCART IF THE CORRELATION MATRIX PRESENTS ANOMALIES

                // Anomalies
                anomalies = (freqCorrList.array().isNaN()).count();

                // Infinity
                infinity = (freqCorrList.array().isInf()).count();

                if(anomalies>0 || infinity>0){

                    cout<<i<<" anomalies in correlations "<<endl;

                     continue; //skip current iteration

                }

            // DISTANCE FROM THE EMPIRICAL CORRELATION DISTRIBUTION
            correlationDistance = EuclideanDistance(freqCorrList,biomeFreqList);
            
            // MAD
            
            Mean_Abundance_Log = abundance_table.colwise().mean().array();   
            centered =  Mean_Abundance_Log.array() - Mean_Abundance_Log.mean();
            zScoreLogMean = centered.normalized()*sqrt(abundance_table.cols());
            FreqMadList=histogram(zScoreLogMean,binMadList);

            madDistance = EuclideanDistance(FreqMadList,freqMadIdeal);
              
            // TAYLOR LAW 
            Slope = Taylor_Slope(abundance_table);
            SlopesDiff = abs(2-Slope);  

        // LIKELIHOOD
        
        Q1 = correlationDistance/(2*pow(sigmaCorr,2))+ madDistance/(2*pow(sigmaMad,2))+SlopesDiff/(2*pow(sigmaTaylor,2));

        // CHECK ACCEPTANCE HASTING FACTOR

        H = exp(-Q1+Q0);
        q = uniform01_R(gen);
        
        //cout<<"Hasting factor of iteration "<<i<<" is "<<H<<endl;
        //cout<<Q0<<"  "<<Q1<<endl;
        //cout<<correlationDistance/(2*pow(sigmaCorr,2))<<"  "<<madDistance/(2*pow(sigmaMad,2))<<"  "<<SlopesDiff/(2*pow(sigmaTaylor,2))<<endl;
        
        if (q < std::min(1.0, H)) {
            
            save=true;
            //cout<<"Change at iteration "<<i<<" has been ACCEPTED "<<acceptedRatio<<endl;

            M0 = M1;
            Q0 = Q1;
            acceptedChanges++;

            cout<<"Iteration: "<<i<<" time "<<1.0*(omp_get_wtime() - startTime)<<" A.R.: "<<acceptedRatio<<" corDist  "<<correlationDistance<<" slope "<<Slope<<"  "<<correlationDistance/(2*pow(sigmaCorr,2))<<"  "<<madDistance/(2*pow(sigmaMad,2))<<"  "<<SlopesDiff/(2*pow(sigmaTaylor,2))<<endl;
            
        } 

        else{
            continue;
        }

        
        /*
        if(Q1<Qmin){

            Qmin = Q1;
            Mmin = M0;
            freqCorrMin = freqCorrList;  
            FreqMadMin = FreqMadList;
            abundance_tableMin = abundance_table;

        }*/

        /************************/
        /*      SAVE DATA       */
        /************************/
        
        if( (i % numItSave) == 0 && save==true){
        
            // INTERACTION MATRIX
            
            matDir = expOutFold+"/mat_"+as_string(i)+".csv";
            writeCSV(M1,matDir);

            // COVARIANCE MATRIX

            covMat=computeCovarianceMatrix(abundance_table);

            CovMatDir = expOutFold+"/"+"covMat_"+as_string(i)+".csv";

            writeCSV(covMat,CovMatDir);
            

            // NETWORK PROPERTIES

            result=mat_statistics(M0,upTrList,lowTrList);

            C=get<0>(result);
            C_eff=get<1>(result);
            muInt=get<2>(result);
            sigmaInt=get<3>(result);

            parametersDirFile << i <<",";
            parametersDirFile << 1.0*(omp_get_wtime() - startTime) <<",";
            parametersDirFile << acceptedRatio <<",";
            parametersDirFile << Q0 << ",";
            parametersDirFile << correlationDistance <<",";
            parametersDirFile << correlationDistance/(2*sigmaCorr*sigmaCorr) <<",";
            parametersDirFile << madDistance/(2*sigmaMad*sigmaMad)<<",";
            parametersDirFile << SlopesDiff/(2*sigmaTaylor*sigmaTaylor)<<",";
            parametersDirFile << Slope << ",";
            parametersDirFile << C <<",";
            parametersDirFile << C_eff <<",";
            parametersDirFile << muInt <<",";
            parametersDirFile << sigmaInt;

            parametersDirFile << "\n";

            parametersDirFile.flush();

        } // Close save data 
        
        

    } //END ITERATION (MAIN) LOOP 
    

/*
    MatrixXd MadDist (FreqMadList.size(),2);

    MadDist << binMadList,FreqMadList;

    writeCSV(MadDist,expOutFold+"/MADDist.csv");

    MatrixXd CorrDist (freqCorrList.size(),2);

    cout<< binCorrList.head(binCorrList.size() - 1).size()<<"   "<<freqCorrList.size()<<endl;

    CorrDist << binCorrList.head(binCorrList.size() - 1),freqCorrList;

    writeCSV(CorrDist,expOutFold+"/CorrDist.csv");

    writeCSV(zScoreLogMean,expOutFold+"/Mean_Log.csv");

    */

    //writeCSV(abundance_table,expOutFold+"/abundance_mat.csv");

    //writeCSV(abundance_table,"/home/jose/MEGA/CPP_MCMC/abundance-table.csv"); 
    


}
