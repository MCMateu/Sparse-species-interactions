#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <tuple>
#include <filesystem>
#include <random>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include<iomanip>
#include <regex>

using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
namespace fs = std::filesystem;

//Function to multuply two matrix
MatrixXd matrixMultiply(const MatrixXd& mat1, const MatrixXd& mat2) {
    // Check if the matrices' dimensions are compatible for multiplication
    if (mat1.cols() != mat2.rows()) {
        cerr << "Error: Incompatible matrix dimensions for multiplication!" << endl;
        exit(1);
    }

    MatrixXd result = mat1 * mat2;
    return result;
}




// Function to split a string into a vector of doubles
VectorXd split(const string& line, char delimiter) {

    /*
    INPUT: 

    -string: each line of the file

    -delimiter: delimiter separating the values in "string", 
                for a .csv it is a ","

    OUTPUT: 

    -tokens: vector whose elements are the numbers in "string" seperated by "delimiter"

    */


    stringstream ss(line);
    string cell;
    VectorXd row;
    while (getline(ss, cell, delimiter)) {
        row.conservativeResize(row.size() + 1);
        try {
            row(row.size() - 1) = stod(cell);
        } catch (const std::invalid_argument& e) {
            // Handle non-numeric values (e.g., replace with NaN)
            row(row.size() - 1) = std::numeric_limits<double>::quiet_NaN();
        }
    }
    return row;
}

// Function to read the CSV file and populate the matrix
MatrixXd readCSV(const string& filename) {
    
    ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    string line;
    vector<VectorXd> rows;

    while (getline(file, line)) {
        VectorXd row = split(line, ',');
        rows.push_back(row);
    }

    file.close();

    // Check if the rows have consistent sizes
    int numColumns = rows[0].size();
    for (int i = 1; i < rows.size(); ++i) {
        if (rows[i].size() != numColumns) {
            throw std::runtime_error("Inconsistent row sizes in the CSV file: " + filename);
        }
    }

    MatrixXd matrix(rows.size(), numColumns);
    for (int i = 0; i < rows.size(); i++) {
        matrix.row(i) = rows[i];
    }

    return matrix;
}

// Function to write a CSV file
void writeCSV(const MatrixXd& data, const string& filename) {

    ofstream file(filename);
    if (!file) {
        cout << "Error opening file: " << filename << endl;
        return;
    }

    for (int i = 0; i < data.rows(); ++i) {
        for (int j = 0; j < data.cols(); ++j) {
            file << data(i, j);
            if (j != data.cols() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
    cout << "Data written to " << filename << " successfully." << endl;
}



// FUNCTION TO RETURN EVENLY SPACED NUMBERS OVER A SPECIFIED INTERVAL
VectorXd linspace(double start, double end, int num) {
    /*
    INPUT: 
    -start: The starting value of the sequence.
    -end: The end value of the sequence
    -num: Number of samples to generate
    OUTPUT: 
    -binCorrList: num evenly spaced samples, calculated over the interval [start, stop].
    */

    
    
    VectorXd binCorrList(num);
    double stepSize = (end - start) / (num - 1);
    
    for (int i = 0; i < num; i++) {
        double binValue = start + i * stepSize;
        binCorrList(i) = binValue;
    }
    
    return binCorrList;
}

// FUNCTION TO RETURN EVENLY SPACED NUMBERS OVER A SPECIFIED INTERVAL WITH A GIVEN SPACING
VectorXd arrange(double start, double end, double step) {
    /*
    INPUT: 
    -start: The starting value of the sequence.
    -end: The end value of the sequence
    -step: Spacing between values.
    OUTPUT: 
    -binCorrList: num evenly spaced samples, calculated over the interval [start, stop].
    */
    
    int num = static_cast<int>((end - start) / step);
    VectorXd binCorrList(num);
    
    for (int i = 0; i < num; i++) {
        double binValue = start + i * step;
        binCorrList(i) = binValue;
    }
    
    return binCorrList;
}

// Function to remove the first row and first column from a MatrixXd object
MatrixXd removeFirstRowAndColumn(MatrixXd& matrix) {
    // Extract the submatrix excluding the first row and column
    return matrix.block(1, 1, matrix.rows() - 1, matrix.cols() - 1);
}


//FUNCTION TO COMPUTE THE ABUNDANCE CORRELATION DISTRIBUTION FOR A GIVEN BIOME
tuple<VectorXd, double, MatrixXd> EmpiricalCorrelationDistribution(const string& biomeDataDir, const VectorXd& binsList, bool header) {
    
    MatrixXd correlationMatrix;

    if(header == true){

        // Read correlation matrix from CSV
        MatrixXd original_correlationMatrix = readCSV(biomeDataDir);

        correlationMatrix = removeFirstRowAndColumn(original_correlationMatrix );

    }else{

         // Read correlation matrix from CSV
        correlationMatrix = readCSV(biomeDataDir);

    }

    // Species number
    double speciesNum = correlationMatrix.rows();

    // Correlation list
    vector<double> correlationList;
    for (int i = 0; i < speciesNum; i++) {
        for (int j = 0; j < i; j++) {
            correlationList.push_back(correlationMatrix(i, j));
        }
    }

    // Histogram
    vector<double> freqList(binsList.size() - 1, 0.0);
    for (double corr : correlationList) {
        size_t binIndex = 0;
        while (binIndex < binsList.size() - 1 && corr > binsList[binIndex + 1]) {
            binIndex++;
        }
        freqList[binIndex] += 1.0;
    }

    // Normalize frequencies
    double h=binsList(1)-binsList(0);
    double total = correlationList.size();
    for (double& freq : freqList) {
        freq /= total*h;
    }

    // Convert freqList to Eigen::VectorXd
    VectorXd eigenFreqList = Eigen::Map<VectorXd>(freqList.data(), freqList.size());

    // Convert speciesNum to Eigen::Index
    Eigen::Index eigenSpeciesNum = static_cast<Eigen::Index>(speciesNum);

    return make_tuple(eigenFreqList, eigenSpeciesNum, correlationMatrix);
}


// Function to calculate the Gaussian distribution
VectorXd gaussianDistribution(const VectorXd& binsList, double mean, double sigma) {
    
    VectorXd distribution(binsList.size());

    for (size_t i = 0; i < binsList.size(); i++) {
        double bin = binsList[i];
        double exponent = -pow(bin - mean, 2) / (2 * pow(sigma, 2));
        double value = exp(exponent) / (sigma * sqrt(2 * M_PI));
        distribution[i] = value;
    }

    return distribution;
}


//FUNCTION TO ROUND THE DECIMALS OF A REAL NUMBER

double roundDecimals(double number, int decimalPlaces) {

    double multiplier = std::pow(10, decimalPlaces);
    double rounded = std::round(number * multiplier) / multiplier;
    return rounded;
}

//FUNCTION TO CONVERT A NUMBER INTO A STRING VARIABLE

string as_string(double number) {
    ostringstream oss;
    oss << fixed << setprecision(8) << number;
    string str = oss.str();

    // Remove trailing zeros
    size_t lastNonZero = str.find_last_not_of('0');
    if (lastNonZero != string::npos) {
        if (str[lastNonZero] == '.') {
            str.resize(lastNonZero);
        } else {
            str.resize(lastNonZero + 1);
        }
    }

    return str;
}


//FUNCTION TO CREATE A RANDOM GAUSSIAN VECTOR
VectorXd generateGaussianNoise(double mu, double sigma, int S) {

    thread_local random_device rd;
    thread_local mt19937 gen(rd() + omp_get_thread_num());
    normal_distribution<double> dist(mu, sigma);

    VectorXd noise(S);
    
    for (int i = 0; i < S; i++) {
        noise(i) = dist(gen);
    }

    return noise;
}

//Function to get the leading eigenvalue
double LeadingEigenvalue(const MatrixXd& M) {

    /*
    Leading eigenvalue of a given matrix;
    
    INPUT
    - "M" (2d array): matrix of real numbers ;
    
    OUTPUT
    - "leadingEigenvalue" (float): maximum real part of the spectrum of M; 
    */

    Eigen::EigenSolver<Eigen::MatrixXd> solver(M);
    const auto& eigenvalues = solver.eigenvalues();
    double maxEigenvalue = eigenvalues.real().maxCoeff();
    return maxEigenvalue;

}



int ExtinctionsNumber(const MatrixXd& M) {

    /*
    
    Number of extinctions in the Lotka-Volterra dynamics associated to a given Network;
    
    INPUT
    - "M" (2d array): adjacency matrix of interactions;
    
    OUTPUT
    - "extinctionsNum" (int): maximum real part of the spectrum of M; 
    
    */
 
    // Species number
    int S = M.rows();
    
    // Stationary state
    VectorXd steadyState = (-M.inverse()) * VectorXd::Ones(S);
    
    // Extinctions number
    int extinctionsNum = (steadyState.array() <= 0).count();

    //cout<<steadyState.array().transpose()<<endl;

    return extinctionsNum;
}

//FUNCTION TO EXTRACT THE NUMBERS IN A STRING
VectorXd get_nums(const string& input) {
    VectorXd numbers;
    regex num_regex(R"([-+]?\d*\.\d+|\d+)");

    auto words_begin = sregex_iterator(input.begin(), input.end(), num_regex);
    auto words_end = sregex_iterator();

    for (sregex_iterator it = words_begin; it != words_end; ++it) {
        double num = stod((*it).str());
        numbers.conservativeResize(numbers.size() + 1);
        numbers(numbers.size() - 1) = num;
    }

    return numbers;
}


Eigen::Matrix<std::string, Eigen::Dynamic, 1> listDir(const string& directory_path) {

    // Get the number of directories first to set the initial size of folder_names
    int num_dirs = 0;
    for (const auto& entry : fs::directory_iterator(directory_path)) {
        if (fs::is_directory(entry.path())) {
            num_dirs++;
        }
    }
   
    // Initialize folder_names with the correct size
    Eigen::Matrix<std::string, Eigen::Dynamic, 1> folder_names(num_dirs, 1);

    int index = 0;
    for (const auto& entry : fs::directory_iterator(directory_path)) {
        if (fs::is_directory(entry.path())) {
            folder_names(index++) = entry.path().filename().string();
        }
    }

    return folder_names;
}


Eigen::Matrix<std::string, Eigen::Dynamic, 1> listFiles(const std::string& directory_path) {

    // Get the number of files first to set the initial size of file_names
    int num_files = 0;
    for (const auto& entry : fs::directory_iterator(directory_path)) {
        if (fs::is_regular_file(entry.path())) {
            num_files++;
        }
    }

    // Initialize file_names with the correct size
    Eigen::Matrix<std::string, Eigen::Dynamic, 1> file_names(num_files, 1);

    int index = 0;
    for (const auto& entry : fs::directory_iterator(directory_path)) {
        if (fs::is_regular_file(entry.path())) {
            file_names(index++) = entry.path().filename().string();
        }
    }

    return file_names;
}

vector<pair<int, int>> triu_indices(int rows, int cols, int k = 0) {
    vector<pair<int, int>> indices;
    for (int j = k; j < cols; ++j) {
        for (int i = j + 1; i < rows; ++i) {
            indices.emplace_back(i, j);
        }
    }
    return indices;
}

vector<pair<int, int>> tril_indices(int rows, int cols, int k = 0) {
    vector<pair<int, int>> indices;
    for (int j = 0; j < cols; ++j) {
        for (int i = k; i < rows; ++i) {
            indices.emplace_back(i, j);
        }
        ++k;
    }
    return indices;
}

MatrixXd PearsCorrMat(const MatrixXd& crossSecAbundances){

    VectorXd Mean = crossSecAbundances.colwise().mean();
    MatrixXd centered = crossSecAbundances.transpose().colwise()-Mean; 
    MatrixXd normalized = centered.rowwise().normalized();

    // compute covariance matrix
    MatrixXd pearsCoeff = normalized * normalized.transpose();
    
    return pearsCoeff;

}

VectorXd SyntheticCorrelationDistribution(const MatrixXd& crossSecAbundances,
                                          int speciesNum,
                                          const vector<pair<int, int>>& upTrList,
                                          const VectorXd& binsList) {

    /*
    Correlation distribution of a synthetic community;

    #INPUT:
    - crossSecAbundances (SxS matrix): matrix with abundances of species i (rows) in the community j (columns);
    - speciesNum (int): number of species S;
    - binsList (VectorXd): list of the distribution bins;
    - upTrList (vector<pair<int, int>>): upper triangle indexes.

    #OUTPUT
    - freqList (VectorXd): list of the distribution frequencies.
    */

    VectorXd Mean = crossSecAbundances.colwise().mean();
    MatrixXd centered = crossSecAbundances.transpose().colwise()-Mean; 
    MatrixXd normalized = centered.rowwise().normalized();

    // compute covariance matrix
    MatrixXd pearsCoeff = normalized * normalized.transpose();

    // extract upper triangle elements (excluding the diagonal)
    VectorXd correlationList(upTrList.size());
    
    //#pragma omp parallel for
    for (int i = 0; i < upTrList.size(); ++i) {
        int row = upTrList[i].first;
        int col = upTrList[i].second;
        correlationList(i) = pearsCoeff(row, col);
    }

    VectorXd freqList = VectorXd::Zero(binsList.size() - 1);

    
    //#pragma omp parallel for
    for (int i = 0; i < correlationList.size(); ++i) {
        double value = correlationList(i);
        for (int j = 0; j < binsList.size()-1; ++j) {
            if (value >= binsList(j) && value < binsList(j + 1)) {
                freqList(j) += 1;
                break;
            }
        }
    }


    // histogram
    /*VectorXd freqList(binsList.size() - 1);
    for (int i = 0; i < correlationList.size(); ++i) {
        double value = correlationList(i);
        auto it = binsList.segment(1, binsList.size() - 1).array() > value;
        int index = it.cast<int>().sum() - 1;
        if (index >= 0 && index < freqList.size()) {
            freqList(index) += 1.0;
        }
    }*/

    // Normalize the histogram by area
    double binWidth = binsList(1) - binsList(0);
    freqList /= roundDecimals(correlationList.size()*binWidth,4);

    return freqList;
}

// FUNCTION TO PERTURBE A RANDOM ELEMENT 

MatrixXd MatrixElementPerturbation(const MatrixXd& matIn, int S, double epsilon) {

    /*
    Random perturbation on a random matrix element;

    INPUT
    - matIn (SxS matrix): input matrix;
    - S (int): number of species;
    - epsilon (double): range of random perturbation distribution;

    OUTPUT
    - matOut (SxS matrix): perturbated matrix;
    */

    // Output matrix

    MatrixXd matOut = matIn;

    // Random element
    random_device rd;
    mt19937 gen(rd()+omp_get_thread_num());
    uniform_int_distribution<> dis(0, S - 1);
    int row1, col1;
    
    // Ensure row1 and col1 are different
    do {
        row1 = dis(gen);
        col1 = dis(gen);
    } while (row1 == col1);


    // Random perturbation
    uniform_real_distribution<> perturb_dist(-epsilon, epsilon);
    double c = perturb_dist(gen);
    matOut(row1, col1) += c;

    return matOut;
}


double EuclideanDistance(const VectorXd& List1, const VectorXd& List2) {

    /*
    Euclidean distance (square of) between two vectors;

    INPUT
    - List1/List2 (1d array): vectors among which measure the distance;

    OUTPUT
    - distance (float): square of the euclidean distance; 
    */

    // Calculate the squared Euclidean distance
    VectorXd diff = List1 - List2;
    double distance = diff.squaredNorm();

    return distance;
}

double EuclideanDistanceLog(const VectorXd& List1, const VectorXd& List2) {

    double epsilon=0.1;

    VectorXd List1_corrected = (List1.array() + epsilon).log();
    VectorXd List2_corrected = (List2.array() + epsilon).log();
    
    VectorXd squaredDifferences = (List1_corrected - List2_corrected).array().square();
    double distance = squaredDifferences.sum();
    
    return distance;
}

// FUNCTION TO COMPUTE THE COVARIANCE MATRIX

MatrixXd computeCovarianceMatrix(const Eigen::MatrixXd& abundance_table) {


    int numSamples = abundance_table.cols();

    VectorXd Mean = abundance_table.transpose().rowwise().mean();
    MatrixXd centered = abundance_table.transpose().colwise()-Mean; 

    // compute covariance matrix
    MatrixXd covMat = centered * centered.transpose() / ( numSamples - 1 );

    return covMat;

}

//FUNCTION TO RETURN LOG SPACED NUMBERS OVER SPECIED INTERVAL
VectorXd logspace(double start, double end, int num) {

    ArrayXd logspaceArray = ArrayXd::LinSpaced(num, log10(start), log10(end));
    return pow(10,logspaceArray);

}

double Taylor_Slope(const MatrixXd& abundance_table){

    
    // MEAN

    VectorXd mean = abundance_table.colwise().mean();

    // VARIANCE

    int num_samples = abundance_table.rows();
    VectorXd var = (abundance_table.rowwise() - mean.transpose()).array().square().colwise().sum() / num_samples;

    // FILTER VALUES 
    double filter=0.1;
    
    if( (mean.array() > filter).count() == 0){
    
    	return 0;
    
    }

    ArrayXd mean_f ((mean.array() > filter).count());
    ArrayXd var_f ((mean.array() > filter).count());  

    int k=0;
    for(int i = 0; i < mean.size() ; ++i){

        if (mean(i)>filter){
            
             mean_f(k)=mean(i);
             var_f(k)=var(i);
             k++;
        }
        
    }
    
    
    
    //ArrayXd mean_f = mean;
    //ArrayXd var_f = var;
    

    // BINNING
    int binsNum=15;
    VectorXd bins = logspace(mean_f.minCoeff(),mean_f.maxCoeff(), binsNum+1);

    VectorXd x_avg = VectorXd::Zero(binsNum);
    VectorXd y_avg = VectorXd::Zero(binsNum);

    for (int i = 0; i < binsNum; ++i) {

        ArrayXd mask = ((mean_f.array() >= bins(i)) && (mean_f.array() < bins(i + 1))).cast<double>(); //boolean mask

       // Compute the sum of elements selected by the mask
       double sum_x = (mean_f.array() * mask).sum();
       double sum_y = (var_f.array() * mask).sum();
        
       // Compute the number of elements selected by the mask (non-zero elements)
       int count = mask.count();
        
       x_avg(i) = sum_x / count;
       y_avg(i) = sum_y / count;

    }


    // PERFORM LINEAR REGRESSION USING EIGEN

    VectorXd lx_avg = (x_avg.array() > 0.0).select(x_avg.array().log10(), -std::numeric_limits<double>::infinity());
    VectorXd ly_avg = (y_avg.array() > 0.0).select(y_avg.array().log10(), -std::numeric_limits<double>::infinity());

    lx_avg = lx_avg.array().isFinite().select(lx_avg, 0);
    ly_avg = ly_avg.array().isFinite().select(ly_avg, 0);

    // Create boolean masks to identify zero values
    VectorXd lx_mask = (lx_avg.array() != 0).cast<double>();
   // Count the non-zero elements
    int numNonZero = lx_mask.count();

    // Filter and reshape the vectors using the masks
    VectorXd lx_filtered(numNonZero);
    VectorXd ly_filtered(numNonZero);
    int count = 0;
    for (int i = 0; i < lx_avg.size(); ++i) {
        if (lx_mask[i]!=0) {
            lx_filtered(count) = lx_avg(i);
            ly_filtered(count) = ly_avg(i);
            count++;
        }
    }

    //cout<<lx_avg.transpose()<<endl<<lx_mask.transpose()<<endl<<lx_filtered.transpose()<<endl;

    VectorXd ones = VectorXd::Ones(lx_filtered.size());
    MatrixXd A(lx_filtered.size(), 2);
    A << ones, lx_filtered;

    VectorXd params = A.colPivHouseholderQr().solve(ly_filtered);

    double slope=params(1);

    /*int numRows = lx_filtered.size();
    int numCols = 2;

    MatrixXd matrix(numRows, numCols);
    for (int i = 0; i < numRows; i++) {
        matrix(i, 0) = lx_filtered(i);
        matrix(i, 1) = ly_filtered(i);
    }

    writeCSV(matrix,"/home/jose/MEGA/CPP_MCMC/Taylor.csv");*/

    return slope;


}

VectorXd histogram(const VectorXd& data, const VectorXd& bins) {

    int binsNum = bins.size();
    VectorXd hist = VectorXd::Zero(binsNum);

    //#pragma omp parallel for
    for (int i = 0; i < data.size(); ++i) {
        double value = data(i);
        for (int j = 0; j < binsNum-1; ++j) {
            if (value >= bins(j) && value < bins(j + 1)) {
                hist(j) += 1;
                break;
            }
        }
    }

    hist /= (data.size() * (bins(1) - bins(0)));

    return hist;
}


MatrixXd joinVectors(const VectorXd& vector1, const VectorXd& vector2) {
    if (vector1.size() != vector2.size()) {
        cout << "Error: Vectors must have the same size." << endl;
        return MatrixXd();
    }

    int numRows = vector1.size();
    int numCols = 2;

    MatrixXd matrix(numRows, numCols);
    for (int i = 0; i < numRows; i++) {
        matrix(i, 0) = vector1(i);
        matrix(i, 1) = vector2(i);
    }

    return matrix;
}

VectorXd flattenMatrix(const MatrixXd& matrix) {
    int rows = matrix.rows();
    int cols = matrix.cols();
    VectorXd flattened(rows * cols);

    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < rows; ++i) {
            flattened(i * cols + j) = matrix(i, j);
        }
    }

    return flattened;
}
