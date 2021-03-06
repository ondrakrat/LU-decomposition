#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <atomic>

using namespace std;
using namespace std::chrono;

class LU {
	public:
		std::mutex mutex;
		int workerCount = 1;
		std::condition_variable cv;
		size_t worker = 0;
        std::atomic<int> atomic;
        std::atomic<int> atomic2;
		// It reads a matrix from the binary file.
		void readMatrixFromInputFile(const string& inputFile)	{
			ifstream bin(inputFile.c_str(), ifstream::in | ifstream::binary);
			if (bin)	{
				uint64_t n = 0;
				bin.read((char*) &n, sizeof(uint64_t));
				A.resize(n, vector<double>(n, 0.0));
				L = U = A;
				for (uint64_t r = 0; r < n; ++r)
					bin.read((char*) A[r].data(), n*sizeof(double));
			} else {
				throw invalid_argument("Cannot open the input file!");
			}

			bin.close();
		}

		double decompose()	{
			high_resolution_clock::time_point start = high_resolution_clock::now();

            // LU decomposition
			decompose_parallel();
//            decompose_sequential();

			double runtime = duration_cast<duration<double>>(high_resolution_clock::now()-start).count();

			return runtime;
		}

		void writeResults(const string& outputFile)	{
			ofstream bout(outputFile.c_str(), ofstream::out | ofstream::binary | ofstream::trunc);
			if (bout)	{
				uint64_t n = A.size();
				for (uint64_t r = 0; r < n; ++r)
					bout.write((char*) L[r].data(), n*sizeof(double));
				for (uint64_t r = 0; r < n; ++r) {
//                    double column[n];
                    std::vector<double> column(n);
                    for (uint64_t q = 0; q < n; ++q) {
                        column[q] = U[q][r];
                    }
//                    cout << (char *) column.data() << endl;
                    bout.write((char *) column.data(), n * sizeof(double));
//                    bout.write((char *) U[r].data(), n * sizeof(double));
                }
			} else {
				throw invalid_argument("Cannot open the input file!");
			}

			bout.close();
		}
		void barrier() {
//			std::unique_lock<std::mutex> lock{mutex};
//			worker++;
//			if (worker == workerCount) {
//				cv.notify_all();
//				worker = 0;
//			} else {
//				cv.wait(lock);
//			}

            atomic++;
            while (atomic < workerCount) {
            }
            atomic++;
            if (atomic == (2 * workerCount)) {
                atomic = 0;
            }
		}
    void barrier2() {
        atomic2++;
        while (atomic2 < workerCount) {
        }
        atomic2++;
        if (atomic2 == (2 * workerCount)) {
            atomic2 = 0;
        }
    }

	private:

		vector<vector<double>> A, L, U;
		friend ostream& operator<<(ostream&, const LU&);
        void decompose_sequential() {
            // TODO: handle division by zero by pivoting the A matrix
            for (int k = 0; k < A.size(); ++k) {
                for (int j = k; j < A.size(); ++j) {
                    U[k][j] = A[k][j];
                }
                L[k][k] = 1;
                for (int i = k + 1; i < A.size(); ++i) {
                    L[i][k] = A[i][k] / U[k][k];
                }
                for (int i = k + 1; i < A.size(); ++i) {
                    for (int j = k + 1; j < A.size(); ++j) {
                        A[i][j] = A[i][j] - (L[i][k] * U[k][j]);
                    }
                }
            }
        }
		void decompose_parallel() {
            atomic = 0;
            atomic2 = 0;
            workerCount = min(thread::hardware_concurrency(), (unsigned int) A.size());
			vector<thread> workers;
			for (int i = 0; i < workerCount; ++i) {
				workers.emplace_back(&LU::decompose_parallel_job, this, i);
			}
			for (thread &worker : workers) {
				worker.join();
			}
		}
        double decompose_parallel_job(int threadId) {
			int row_start = (threadId * A.size()) / workerCount;
			int row_end = row_start + (A.size() / workerCount);

            for (int k = 0; k < A.size(); ++k) {
				for (int j = max(k, row_start); j <= row_end; ++j) {
                    if (j < A.size()) {
                        computeU(k, j);
                    }
				}
				barrier();
				L[k][k] = 1;
				for (int i = max(k, row_start); i <= row_end; ++i) {
                    if (i < A.size()) {
                        computeL(k, i);
                    }
				}
				barrier2();
            }
			return 0;
        }
		void computeU(int k, int j) {
            double sum = 0;
            for (int r = 0; r < k; ++r) {
                sum += (L[k][r] * U[j][r]);
            }
            U[j][k] = A[k][j] - sum;
		}
		void computeL(int j, int k) {
			double sum = 0;
			// TODO: j - 1?
			for (int r = 0; r < j; ++r) {
				sum += L[k][r] * U[j][r];
			}
			// TODO: handle division by zero?
			L[k][j] = (1 / U[j][j]) * (A[k][j] - sum);
		}
};

// Print the matrices A, L, and U in an instance of LU class.
ostream& operator<<(ostream& out, const LU& lu)	{
	function<void(const vector<vector<double>>&)> printMatrix = [&](const vector<vector<double>>& M)	{
		int n = M.size();
		for (int i = 0; i < n; ++i)	{
			for (int j = 0; j < n; ++j)
				out<<" "<<setw(10)<<M[i][j];
			out<<endl;
		}
	};

	out<<"Matrix A:"<<endl;
	printMatrix(lu.A);
	out<<endl<<"Lower matrix:"<<endl;
	printMatrix(lu.L);
	out<<endl<<"Upper matrix:"<<endl;
	printMatrix(lu.U);

	return out;
}

int main(int argc, char* argv[])	{
	if (argc <= 1 || argc > 3)	{
		cout<<"LU decomposition of a square matrix."<<endl;
		cout<<endl<<"Usage:"<<endl;
		cout<<"\t"<<argv[0]<<" inputMatrix.bin [output.bin]"<<endl;
		return 1;
	}

	string inputFile = argv[1], outputFile;
	if (argc == 3)
		outputFile = argv[2];

	LU lu;
	lu.readMatrixFromInputFile(inputFile);
	double totalDuration = lu.decompose();
	// Decomposition is printed only if the output file is not written.
	if (outputFile.empty()) {
        cout << lu << endl;
    } else {
//        cout << lu << endl;
        lu.writeResults(outputFile);
    }

	cout<<"computational time: "<<totalDuration<<" s"<<endl;

	return 0;
}

