#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include "mc_ising.hpp"

double getCorr(int NN, const std::vector<int> &spins0, const std::vector<int> &spins){
    double res = 0.0;
    for(int i=0; i<NN; i++){
        res += spins0[i] * spins[i];
    }
    return res / NN;
}

int main(int argc, char *argv[]) {
    const int N = 100;
    const int NN = N*N;
    const int Nwarm0 = 50000;
    const int Nwarm = Nwarm0 * NN;
    const int Nloop0 = 1000;
    const double K = 0.2;
    const double Jeff = 1.1064;
    const double kBTc = 2.0 / log(1.0 + sqrt(2.0)) * Jeff;
    const double beta = 1.0 / kBTc;
    double dmag;
    std::vector<int> spins(NN);
    std::vector<int> spins0(NN);
    std::string fout;
    std::ofstream Fout;

    initSpins(N, spins0);
    for(int iwarm = 0; iwarm < Nwarm; iwarm++){
        updateLocal(N, spins0, beta, K, dmag);
    }

    spins = spins0;
    fout = "ac/local_ac";
    Fout.open(fout);
    for(int iloop = 0; iloop < Nloop0; iloop++){
        for(int i=0; i<NN; i++) updateLocal(N, spins, beta, K, dmag);
        Fout << iloop << " " << getCorr(NN, spins0, spins) << std::endl;
    }
    Fout.close();

    spins = spins0;
    fout = "ac/wolff_ac";
    Fout.open(fout);
    for(int iloop = 0; iloop < Nloop0; iloop++){
        for(int i=0; i<NN; i++) updateWolff(N, spins, beta, K, dmag);
        Fout << iloop << " " << getCorr(NN, spins0, spins) << std::endl;
    }
    Fout.close();

    spins = spins0;
    fout = "ac/slmc_ac";
    Fout.open(fout);
    initSpins(N, spins);
    for(int iloop = 0; iloop < Nloop0; iloop++){
        for(int i=0; i<NN; i++) updateSLMC(N, spins, beta, K, dmag);
        Fout << iloop << " " << getCorr(NN, spins0, spins) << std::endl;
    }
    Fout.close();
}
