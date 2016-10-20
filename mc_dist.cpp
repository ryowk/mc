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

    int N, Nwarm0, Nloop0, NkBT;
    double K;
    Fin >> N >> Nwarm0 >> Nloop0 >> NkBT >> K;

    int NN = N * N;
    int Nwarm = Nwarm0 * NN;
    double ukBT = 3.5, bkBT = 1.5;
    double dkBT = (ukBT-bkBT) / NkBT;

    std::vector<int> spins(NN);
    initSpins(N, spins);

    for (int ikBT = NkBT; ikBT >= 0; ikBT--) {
        double kBT = dkBT * ikBT + bkBT;
        double beta = 1.0 / kBT;

        // ウォーミングアップ
        for (int iwarm = 0; iwarm < Nwarm; iwarm++) {
            double dmag;
            updateLocal(N, spins, beta, K, dmag);
        }
        // サンプリング
        for (int iloop = 0; iloop < Nloop0; iloop++) {
            for (int i = 0; i < NN; i++) {
                double dmag;
                updateLocal(N, spins, beta, K, dmag);
            }
            double corr, ene;
            getInOut(N, spins, K, corr, ene);
            Fout << corr << " " << ene << std::endl;
        }
    }
}
