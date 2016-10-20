#pragma once
#include <vector>
#include <random>


void initSpins(int N, std::vector<int> &spins);
void printSpins(int N, const std::vector<int> &spins);
void updateLocal(int N, std::vector<int> &spins, double beta, double K, double &dmag);
void getCluster(int N, std::vector<int> &spins, double beta, double K, int x, int y, std::vector<bool> &visited, std::vector<bool> &cluster); 
void updateWolff(int N, std::vector<int> &spins, double beta, double K, double &dmag);
void updateSLMC(int N, std::vector<int> &spins, double beta, double Jeff, double K, double &dmag);
void getQuantities(int N, const std::vector<int> &spins, double &mag);
void getInOut(int N, const std::vector<int> &spins, double K, double &corr, double &ene);
double getEnergy(int N, const std::vector<int> &spins, double J, double K);
