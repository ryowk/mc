#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include "mc_ising.hpp"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cout << "Input a file name." << std::endl;
        return 0;
    }
    std::string fname = argv[1];
    std::string fin = fname + ".in";
    std::string fout = fname + ".out";
    std::ifstream Fin(fin);
    std::ofstream Fout(fout);

    int N, Nwarm0, Nloop0, NkBT, method;
    double K;
    Fin >> N >> Nwarm0 >> Nloop0 >> NkBT >> K >> method;

    void (*update)(int, std::vector<int>&, double, double, double&);
    switch(method){
        case(0):
            update = updateLocal;
            break;
        case(1):
            update = updateWolff;
            break;
        case(2):
            update = updateSLMC;
            break;
    }

    int NN = N * N;
    int Nwarm = Nwarm0 * NN;
    int Nloop = Nloop0 * NN;
    double dkBT = 3.0 / NkBT;

    std::vector<int> spins(NN);
    initSpins(N, spins);

    for (int ikBT = NkBT; ikBT > 0; ikBT--) {
        double kBT = dkBT * ikBT;
        double beta = 1.0 / kBT;
        double mag;
        double magtot = 0.0;
        double mag2tot = 0.0;
        getQuantities(N, spins, mag);
        // ウォーミングアップ
        for (int iwarm = 0; iwarm < Nwarm; iwarm++) {
            double dmag;
            update(N, spins, beta, K, dmag);
            mag += dmag;
        }
        // サンプリング
        for (int iloop = 0; iloop < Nloop; iloop++) {
            double dmag;
            update(N, spins, beta, K, dmag);
            mag += dmag;
            magtot += fabs(mag);
            mag2tot += mag * mag;
        }
        magtot /= static_cast<double>(Nloop);
        mag2tot /= static_cast<double>(Nloop);
        double var = mag2tot - magtot * magtot;
        double std = sqrt(var);
        Fout << kBT << " " << magtot << " " << std << std::endl;
    }
}
