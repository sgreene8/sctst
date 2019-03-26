// Implementation of direct count algorithm described in Nguyen & Barker (2010), dx.doi.org/10.1021/jp100132s
// Usage: ./doloops.x [file containing frequencies and anharmonicities]
// Output format: (state energy, state coupling to reaction coordinate)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Emax (30 * 0.02)

int nvmax(int nv[], int k, double Eu, double anh[], double freq[], int n);
void findConfigs(double freq[], double anh[], double coupling[], int n, int *NN, int idx, FILE *outFile);

int main(int argc, const char * argv[]) {
    unsigned int i, j;
    
    FILE *inFile = fopen(argv[1], "r");
    
    char *buf = NULL;
    size_t bufLen;
    getline(&buf, &bufLen, inFile);
    
    int nModes;
    fscanf(inFile, "%d", &nModes);
    
    fscanf(inFile, "%s", buf);
    getline(&buf, &bufLen, inFile);
    
    double freq[nModes];
    for (i = 0; i < nModes; i++) {
        fscanf(inFile, "%lf,", &freq[i]);
    }
    
    fscanf(inFile, "%s", buf);
    getline(&buf, &bufLen, inFile);
    
    double coupling[nModes];
    for (i=0; i < nModes; i++) {
        fscanf(inFile, "%lf,", &coupling[i]);
    }
    
    fscanf(inFile, "%s", buf);
    getline(&buf, &bufLen, inFile);
    
    double anh[nModes][nModes];
    for (i=0; i < nModes; i++) {
        for (j=0; j < nModes; j++) {
            fscanf(inFile, "%lf,", &anh[i][j]);
        }
    }
    
    int nv[nModes];
    for (i = 0; i < nModes; i++) {
        nv[i] = 0;
    }
    
    sprintf(buf, "States_Output_%dD.txt", nModes + 1);
    FILE *dataFile = fopen(buf, "w");
    
    findConfigs(freq, (double *)anh, coupling, nModes, nv, 0, dataFile);
}

void findConfigs(double freq[], double anh[], double coupling[], int n, int *NN, int idx, FILE *outFile) {
    double energy = 0;
    unsigned int i, j, k;
    for (i=0; i<n; i++) {
        energy += freq[i] * (NN[i] + 0.5);
        for (j=i; j<n; j++) {
            energy += anh[i * n + j] * (NN[i] + 0.5) * (NN[j] + 0.5);
        }
    }
    double Eu = Emax - energy;
    int maxi = nvmax(NN, idx, Eu, anh, freq, n);
    for (i=0; i <= maxi; i++) {
        NN[idx] = i;
        if (idx == n - 1) {
            energy = 0;
            double nCoupling = 0;
            for (j=0; j<n; j++) {
                energy += freq[j] * (NN[j] + 0.5);
                nCoupling += coupling[j] * (NN[j] + 0.5);
                for (k=j; k<n; k++) {
                    energy += anh[j * n + k] * (NN[j] + 0.5) * (NN[k] + 0.5);
                }
            }
            fprintf(outFile, "%lf, %lf\n", energy, nCoupling);
        }
        else {
            for (j = idx + 1; j < n; j++) {
                NN[j] = 0;
            }
            findConfigs(freq, anh, coupling, n, NN, idx + 1, outFile);
        }
    }
}

int nvmax(int nv[], int k, double Eu, double anh[], double freq[], int n) {
    if (Eu < 0) {
        return -1;
    }
    unsigned int i;
    if (anh[k * n + k] != 0) {
        double sum = 0;
        for (i = 0; i < n; i++) {
            if (i != k) {
                sum += (nv[i] + 0.5) * anh[i * n + k];
            }
        }
        double anhkk = anh[k * n + k];
        double vd = (-freq[k] - sum) / (2 * anhkk) - 0.5;
        double anhksum = 0;
        for (i=0; i < n; i++) {
            if (i != k) {
                anhksum = anh[k * n + i];
            }
        }
        double Dk = -(freq[k] + sum) * (freq[k] + sum) / 4 / anhkk - freq[k] / 2 - anhksum;
        int vmax = vd * (1 - sqrt(1 - Eu / Dk));
        if ((freq[k] + sum) > 0 && anhkk < 0) {
            if (Eu > Dk) {
                return vd;
            }
            else {
                return vmax;
            }
        }
        else {
            if ((freq[k] + sum) > 0 && anhkk > 0) {
                return vmax;
            }
            else {
                return -1;
            }
        }
    }
    
    return 0;
}
