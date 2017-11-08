#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <functional>
#include <vector>
#include <cmath>
#include <thread>
#include <mutex>
#include <condition_variable>

using namespace std;
using namespace std::chrono;

//#define debug
class LU {
public:
    // It reads a matrix from the binary file.
    void readMatrixFromInputFile(const string &inputFile) {
        ifstream bin(inputFile.c_str(), ifstream::in | ifstream::binary);
        if (bin) {
            uint64_t n = 0;
            bin.read((char *) &n, sizeof(uint64_t));
            A.resize(n, vector<double>(n, 0.0));
            L = U = A;
            for (uint64_t r = 0; r < n; ++r)
                bin.read((char *) A[r].data(), n * sizeof(double));
        } else {
            throw invalid_argument("Cannot open the input file!");
        }

        bin.close();
    }

    int numOfWorkers = thread::hardware_concurrency();

    double decompose() {
        vector<thread> workers;
        high_resolution_clock::time_point start = high_resolution_clock::now();
        // TODO: Implement a parallel LU decomposition...\

        for (int jobId = 0; jobId < numOfWorkers; ++jobId)
            workers.push_back(thread(&LU::LU_parallel, this, jobId));
        for (thread &worker : workers)
            worker.join();
        double runtime = duration_cast<duration<double>>(high_resolution_clock::now() - start).count();
#ifdef debug
        checkMatrix();
#endif
        return runtime;
    }

    std::mutex mutex;
    std::condition_variable cv;
    size_t worker = 0;

    void barrier_wait() {
        std::unique_lock<std::mutex> lock{mutex};
        worker++;
        if (worker == numOfWorkers) {
            //cout << "waking by no" << worker << endl;
            cv.notify_all();
            worker = 0;

        } else {
            //cout << "locked " << worker << endl;
            cv.wait(lock);
        }
    }

    double LU_parallel(int id) {
        //parallel//////
        int size = A.size();
        int start, end, pivot;
        start = (id * size) / numOfWorkers;
        end = start + (size / numOfWorkers);//-1;
        //cout << "Thread start id " << id << endl;

        for (int k = 0; k < size; ++k) {
            pivot = max(k, start);//+1
            for (int j = pivot; j <= end; ++j) {
                double sum = 0;
                for (int r = 0; r < k; ++r) {
                    sum += L[k][r] * U[r][j];
                }
                U[k][j] = A[k][j] - sum;
            }
            barrier_wait();

            //cout << "1st barrier passed no " << id << endl;
            L[k][k] = 1;

            for (int i = pivot; i <= end; ++i) {
                if (i >= size) continue;
                double sum = 0;
                for (int r = 0; r < k; ++r) {
                    sum += L[i][r] * U[r][k];
                }
                L[i][k] = (A[i][k] - sum) / U[k][k];
            }
            barrier_wait();

        }
        return 0;
    }

    double LU_sequential() {
        high_resolution_clock::time_point start = high_resolution_clock::now();
        // TODO: Implement a parallel LU decomposition...
        //sequential
        for (int k = 0; k < A.size(); ++k) {
            for (int j = k; j < A.size(); ++j) {
                double sum = 0;
                for (int r = 0; r < k; ++r) {
                    sum += L[k][r] * U[r][j];
                }
                U[k][j] = A[k][j] - sum;
            }
            L[k][k] = 1;
            for (int i = k + 1; i < A.size(); ++i) {
                double sum = 0;
                for (int r = 0; r < k; ++r) {
                    sum += L[i][r] * U[r][k];
                }
                L[i][k] = (A[i][k] - sum) / U[k][k];
            }
            /* for (int i = k+1; i < A.size(); ++i) {
                 for (int j = k+1; j < A.size(); ++j) {
                     double sum=0;
                     for (int r = 0; r < k; ++r) {
                         sum+=L[i][r]*U[r][j];
                     }
                     A[i][j]=(A[i][j]-sum)-L[i][k]*U[k][j];
                 }
             }*/
        }
        double runtime = duration_cast<duration<double>>(high_resolution_clock::now() - start).count();
#ifdef debug
        checkMatrix();
#endif
        return runtime;
    }

#ifdef debug
    void checkMatrix() {
        int n = A.size();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double sum = 0;
                for (int k = 0; k < n; ++k) {
                    sum += L[i][k] * U[k][j];
                }
                double epsilon = 0.001;
                double a = std::fabs(sum - A[i][j]);
                if (a >= epsilon) {
                    cout << "Error at " << i << ", " << j << " is " << sum << " should be " << A[i][j] << endl;
                    cout << "L " << i << ", " << j << " is " << L[i][j] << endl;
                    cout << "U " << i << ", " << j << " is " << U[i][j] << endl;
                }
            }
        }
    }
#endif

    void writeResults(const string &outputFile) {
        ofstream bout(outputFile.c_str(), ofstream::out | ofstream::binary | ofstream::trunc);
        if (bout) {
            uint64_t n = A.size();
            for (uint64_t r = 0; r < n; ++r)
                bout.write((char *) L[r].data(), n * sizeof(double));
            for (uint64_t r = 0; r < n; ++r)
                bout.write((char *) U[r].data(), n * sizeof(double));
        } else {
            throw invalid_argument("Cannot open the input file!");
        }

        bout.close();
    }

private:

    vector<vector<double>> A, L, U, B;

    friend ostream &operator<<(ostream &, const LU &);
};

// Print the matrices A, L, and U in an instance of LU class.
ostream &operator<<(ostream &out, const LU &lu) {
    function<void(const vector<vector<double>> &)> printMatrix = [&](const vector<vector<double>> &M) {
        int n = M.size();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                out << " " << setw(10) << M[i][j];
            out << endl;
        }
    };

    out << "Matrix A:" << endl;
    printMatrix(lu.A);
    out << endl << "Lower matrix:" << endl;
    printMatrix(lu.L);
    out << endl << "Upper matrix:" << endl;
    printMatrix(lu.U);
    return out;
}

int main(int argc, char *argv[]) {
    if (argc <= 1 || argc > 3) {
        cout << "LU decomposition of a square matrix." << endl;
        cout << endl << "Usage:" << endl;
        cout << "\t" << argv[0] << " inputMatrix.bin [output.bin]" << endl;
        return 1;
    }

    string inputFile = argv[1], outputFile;
    if (argc == 3)
        outputFile = argv[2];

    LU lu;
    lu.readMatrixFromInputFile(inputFile);
    double totalDuration = lu.decompose();
    // Decomposition is printed only if the output file is not written.
    if (outputFile.empty())
        //cout<<lu<<endl;

        cout << "Done" << endl;
    else
        lu.writeResults(outputFile);

    cout << "computational time: " << totalDuration << " s" << endl;


    return 0;
}
