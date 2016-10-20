#include "mc_ising.hpp"
#include <iostream>

void initSpins(int N, std::vector<int> &spins) {
    static std::random_device seed_gen;
    static std::mt19937 engine(seed_gen());
    static std::uniform_real_distribution<> drand(0.0, 1.0);
    for (int i = 0; i < N * N; i++) {
        if (drand(engine) < 0.5)
            spins[i] = 1;
        else
            spins[i] = -1;
    }
}

void printSpins(int N, const std::vector<int> &spins) {
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            char ch;
            if (spins[N * y + x] == 1)
                ch = 'O';
            else
                ch = 'X';
            std::cout << ch << " ";
        }
        std::cout << std::endl;
    }
}

// ローカルアップデート
void updateLocal(int N, std::vector<int> &spins, double beta, double K,
                 double &dmag) {
    static std::random_device seed_gen;
    static std::mt19937 engine(seed_gen());
    static std::uniform_real_distribution<> drand(0.0, 1.0);

    dmag = 0.0;
    int x = static_cast<int>(drand(engine) * N);
    int y = static_cast<int>(drand(engine) * N);

    const int dx[8] = {1, 0, -1, 0, 1, -1, -1, 1};
    const int dy[8] = {0, 1, 0, -1, 1, 1, -1, -1};

    int xy[8];
    for (int i = 0; i < 8; i++) {
        int x2 = (x + dx[i] + N) % N;
        int y2 = (y + dy[i] + N) % N;
        xy[i] = N * y2 + x2;
    }

    int s = spins[N * y + x];
    int itempJ = 0, itempK = 0;
    for (int i = 0; i < 4; i++) itempJ += spins[xy[i]];
    for (int i = 0; i < 4; i++)
        itempK += spins[xy[i]] * spins[xy[(i + 1) % 4]] * spins[xy[i + 4]];
    double dE = 2.0 * s * (itempJ + K * itempK);

    if (dE <= 0.0 || drand(engine) < exp(-beta * dE)) {
        spins[N * y + x] = -spins[N * y + x];
        dmag = 2.0 * spins[N * y + x] / static_cast<double>(N * N);
    }
}

void getCluster(int N, std::vector<int> &spins, double beta, double K, int x,
                int y, std::vector<bool> &visitedJ, std::vector<bool> &visitedK,
                std::vector<bool> &cluster) {
    static std::random_device seed_gen;
    static std::mt19937 engine(seed_gen());
    static std::uniform_real_distribution<> drand(0.0, 1.0);

    int pos = N * y + x;
    int s = spins[pos];
    const int dx[8] = {1, 0, -1, 0, 1, -1, -1, 1};
    const int dy[8] = {0, 1, 0, -1, 1, 1, -1, -1};
    // 周りのスピンの座標
    int x2[8], y2[8], xy[8];
    for (int i = 0; i < 8; i++) {
        x2[i] = (x + dx[i] + N) % N;
        y2[i] = (y + dy[i] + N) % N;
        xy[i] = N * y2[i] + x2[i];
    }

    // 原点から生えているボンドの座標
    // bondJ: 2*pos(+x方向), 2*pos+1(+y方向)
    // bondK: pos(+xy方向)
    int bondJ[4] = {2 * pos, 2 * pos + 1, 2 * xy[2], 2 * xy[3] + 1};
    int bondK[4] = {pos, xy[2], xy[6], xy[3]};

    for (int i = 0; i < 4; i++) {
        if (visitedJ[bondJ[i]]) continue;
        visitedJ[bondJ[i]] = true;
        if (spins[xy[i]] == s && drand(engine) < 1.0 - exp(-2.0 * beta)) {
            cluster[xy[i]] = true;
            getCluster(N, spins, beta, K, x2[i], y2[i], visitedJ, visitedK,
                       cluster);
        }
    }

    for (int i = 0; i < 4; i++) {
        if (visitedK[bondK[i]]) continue;
        visitedK[bondK[i]] = true;
        // bondK に含まれる原点以外のスピンの座標
        int j[3] = {i, (i + 1) % 4, i + 4};
        int itemp = spins[xy[j[0]]] * spins[xy[j[1]]] *
                    spins[xy[j[2]]];
        if (itemp == s && drand(engine) < 1.0 - exp(-2.0 * K * beta)) {
            for (int k = 0; k < 3; k++) {
                cluster[xy[j[k]]] = true;
                getCluster(N, spins, beta, K, x2[j[k]], y2[j[k]], visitedJ,
                           visitedK, cluster);
            }
        }
    }
}

// Wolffアルゴリズム
void updateWolff(int N, std::vector<int> &spins, double beta, double K,
                 double &dmag) {
    static std::random_device seed_gen;
    static std::mt19937 engine(seed_gen());
    static std::uniform_real_distribution<> drand(0.0, 1.0);

    int x = static_cast<int>(drand(engine) * N);
    int y = static_cast<int>(drand(engine) * N);
    // (x, y) から出発してクラスターを作る
    static std::vector<bool> visitedJ(N * N * 2);
    static std::vector<bool> visitedK(N * N);
    static std::vector<bool> cluster(N * N);
    fill(visitedJ.begin(), visitedJ.end(), false);
    fill(visitedK.begin(), visitedK.end(), false);
    fill(cluster.begin(), cluster.end(), false);
    cluster[N * y + x] = true;

    getCluster(N, spins, beta, K, x, y, visitedJ, visitedK, cluster);

    dmag = 0.0;
    for (int i = 0; i < N * N; i++) {
        if (cluster[i]) {
            spins[i] = -spins[i];
            dmag += 2.0 * spins[i];
        }
    }
    dmag /= static_cast<double>(N * N);
}

// SLMC
void updateSLMC(int N, std::vector<int> &spins, double beta, double K, double &dmag){
    const double Jeff = 1.1064;

    static std::random_device seed_gen;
    static std::mt19937 engine(seed_gen());
    static std::uniform_real_distribution<> drand(0.0, 1.0);

    int x = static_cast<int>(drand(engine) * N);
    int y = static_cast<int>(drand(engine) * N);
    // (x, y) から出発してクラスターを作る
    static std::vector<bool> visitedJ(N * N * 2);
    static std::vector<bool> visitedK(N * N);
    static std::vector<bool> cluster(N * N);
    fill(visitedJ.begin(), visitedJ.end(), false);
    fill(visitedK.begin(), visitedK.end(), false);
    fill(cluster.begin(), cluster.end(), false);
    cluster[N * y + x] = true;

    getCluster(N, spins, Jeff * beta, 0.0, x, y, visitedJ, visitedK, cluster);

    static std::vector<int> spins2(N*N);
    spins2 = spins;
    dmag = 0.0;
    for(int i=0; i<N*N; i++){
        if(cluster[i]) {
            spins2[i] = -spins[i];
            dmag += 2.0 * spins2[i];
        }
    }
    double EA = getEnergy(N, spins, 1.0, K);
    double EB = getEnergy(N, spins2, 1.0, K);
    double EAeff = getEnergy(N, spins, Jeff, 0.0);
    double EBeff = getEnergy(N, spins2, Jeff, 0.0);
    double dE = (EB - EBeff) - (EA - EAeff);

    if(dE <= 0 || drand(engine) < exp(-beta * dE)){
        spins = spins2;
        dmag /= static_cast<double>(N * N);
    }else{
        dmag = 0.0;
    }
}

void getQuantities(int N, const std::vector<int> &spins, double &mag) {
    mag = 0.0;
    for (int i = 0; i < N * N; i++) mag += spins[i];
    mag /= static_cast<double>(N * N);
}

void getInOut(int N, const std::vector<int> &spins, double K, double &corr, double &ene){
    corr = 0.0;
    ene = 0.0;
    for(int y=0; y<N; y++){
        for(int x=0; x<N; x++){
            int s = spins[N*y+x];
            int x2 = (x+1)%N;
            int y2 = (y+1)%N;
            ene -= s * (spins[N*y+x2] + spins[N*y2+x]);
            ene -= K * s * spins[N*y+x2] * spins[N*y2+x] * spins[N*y2+x2];
            corr += s * (spins[N*y+x2] + spins[N*y2+x]);
        }
    }
    corr /= N * N;
    ene /= N * N;
}

double getEnergy(int N, const std::vector<int> &spins, double J, double K){
    double ene = 0.0;
    for(int y=0; y<N; y++){
        for(int x=0; x<N; x++){
            int s = spins[N*y+x];
            int x2 = (x+1)%N;
            int y2 = (y+1)%N;
            ene -= s * (spins[N*y+x2] + spins[N*y2+x]);
            ene -= K * s * spins[N*y+x2] * spins[N*y2+x] * spins[N*y2+x2];
        }
    }
    return ene;
}
