// Implementation of Wang-Landau algorithm described in Nguyen & Barker (2010), dx.doi.org/10.1021/jp100132s
// Usage: ./adensum.x [file containing frequencies and anharmonicities]
// Output format: (interval energy, average energy of states within interval, average coupling for states within interval)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Emax 0.4
#define Estep 0.0001
#define Niter 21
#define Ntrial 10000000

int nvmax(unsigned int nv[], unsigned int k, double Eu, double anh[], double freq[], unsigned int n);
double nNRG(double *freq, double *anh, unsigned int *nv, unsigned int n);
int ckderiv(double *freq, double *anh, unsigned int *nv, unsigned int n); // returns 0 if ok, 1 if fail
double nCoup(double *coup, unsigned int *nv, unsigned int n);

int main(int argc, const char * argv[]) {
    FILE *inFile = fopen(argv[1], "r");
    
    char *buf;
    size_t bufLen;
    getline(&buf, &bufLen, inFile);
    
    unsigned int nModes;
    fscanf(inFile, "%d", &nModes);
    
    fscanf(inFile, "%s", buf);
    getline(&buf, &bufLen, inFile);
    
    double freq[nModes];
    for (unsigned int i = 0; i < nModes; i++) {
        fscanf(inFile, "%lf,", &freq[i]);
    }
    
    double coupling[nModes];
    fscanf(inFile, "%s", buf);
    getline(&buf, &bufLen, inFile);
    for (unsigned int i=0; i < nModes; i++) {
        fscanf(inFile, "%lf,", &coupling[i]);
    }
    
    double anh[nModes][nModes];
    fscanf(inFile, "%s", buf);
    getline(&buf, &bufLen, inFile);
    for (unsigned int i=0; i < nModes; i++) {
        for (unsigned int j=0; j < nModes; j++) {
            fscanf(inFile, "%lf,", &anh[i][j]);
        }
    }
    
    unsigned int nSteps = Emax / Estep;
    double g[nSteps];
    double TTF[nSteps];
    double Ev[nSteps];
    for (unsigned int i = 0; i < nSteps; i++) {
        g[i] = 0;
        TTF[i] = 0;
        Ev[i] = 0;
    }
    
    double lnf = 1;
    
    double ZPE = 0;
    for (unsigned int i = 0; i < nModes; i++) {
        ZPE += freq[i] / 2;
        for (unsigned int j = i; j < nModes; j++) {
            ZPE += anh[i][j] / 4;
        }
    }
    
    time_t t;
    srand((unsigned) time(&t));
    rand(); // This randomly generated number will always be the same, but subsequent calls to rand() will generate random numbers
    
    // Initialize quantum numbers
    unsigned int nvold[nModes];
    int ntest = 1;
    while (ntest == 1) {
        double Estart = (double)rand() / (double)RAND_MAX * Emax;
        
        for (unsigned int i = 0; i < nModes; i++) {
            nvold[i] = 0;
        }
        
        for (unsigned int l = 0; l < nModes; l++) {
            double Eu = Estart - nNRG(freq, (double *)anh, nvold, nModes);
            if (Eu != 0) {
                int nstart = nvmax(nvold, l, Eu, (double *)anh, freq, nModes);
                if (nstart == -1) {
                    break;
                }
                nvold[l] = nstart * (double)rand() / (double)RAND_MAX;
            }
        }
        ntest = ckderiv(freq, (double *)anh, nvold, nModes);
        if (ntest == 1) {
            printf("Failed derivative test\n");
        }
    }
    double Eold = nNRG(freq, (double *)anh, nvold, nModes) - ZPE;
    unsigned int nold = Eold / Estep;
    
    double p = 1.0 / nModes;
    // Begin algorithm
    for (unsigned int i = 0; i < Niter; i++) {
        printf("Iteration %u of %u\n", i + 1, Niter + 1);
        
        for (unsigned int j = 0; j < Ntrial; j++) {
            unsigned int nvnew[nModes];
            int allZero = 1;
            for (unsigned int k = 0; k < nModes; k++) {
                double r1 = (double)rand() / (double)RAND_MAX;
                if (r1 <= p) {
                    int new = nvold[k] - 1;
                    if (new < 0) {
                        new = 0;
                    }
                    nvnew[k] = new;
                }
                else if (r1 <= 2 * p) {
                    nvnew[k] = nvold[k] + 1;
                }
                else {
                    nvnew[k] = nvold[k];
                }
                if (allZero && nvnew[k] != 0) {
                    allZero = 0;
                }
            }
            double newE = nNRG(freq, (double *)anh, nvnew, nModes) - ZPE;
            if (newE < Emax && ckderiv(freq, (double *)anh, nvnew, nModes) == 0) {
                unsigned int nnew = newE / Estep;
                double r2 = (double)rand() / (double)RAND_MAX;
                if (exp(g[nold] - g[nnew]) > r2) {
                    nold = nnew;
                    for (unsigned int k = 0; k < nModes; k++) {
                        nvold[k] = nvnew[k];
                    }
                    Eold = newE;
                }
            }
            g[nold] += lnf;
        }
        
        lnf *= 0.5;
    }
    
    // Last iteration to compute energies (Ev) & couplings (TTF)
    printf("Starting final iteration\n");
    unsigned int H[nSteps];
    for (unsigned int l = 0; l < nSteps; l++) {
        H[l] = 0;
    }
    
    for (unsigned int j = 0; j < Ntrial; j++) {
        unsigned int nvnew[nModes];
        for (unsigned int k = 0; k < nModes; k++) {
            double r1 = (double)rand() / (double)RAND_MAX;
            if (r1 <= p) {
                int new = nvold[k] - 1;
                if (new < 0) {
                    new = 0;
                }
                nvnew[k] = new;
            }
            else if (r1 <= 2 * p) {
                nvnew[k] = nvold[k] + 1;
            }
            else {
                nvnew[k] = nvold[k];
            }
        }
        double newE = nNRG(freq, (double *)anh, nvnew, nModes) - ZPE;
        if (newE < Emax && ckderiv(freq, (double *)anh, nvnew, nModes) == 0) {
            unsigned int nnew = newE / Estep;
            double r2 = (double)rand() / (double)RAND_MAX;
            if (exp(g[nold] - g[nnew]) > r2) {
                nold = nnew;
                for (unsigned int k = 0; k < nModes; k++) {
                    nvold[k] = nvnew[k];
                }
                Eold = newE;
            }
        }
        H[nold]++;
        g[nold] += lnf;
        Ev[nold] += Eold;
        TTF[nold] += nCoup(coupling, nvold, nModes);
    }
    
    FILE *outFile = fopen("States_Output.txt", "w");
    fprintf(outFile, "Estep = %lf; Niter = %u; Ntrial = %u\n", Estep, Niter, Ntrial);
    fprintf(outFile, "Energy, # of states, average Ev, average coupling\n");
    for (unsigned int i = 0; i < nSteps; i++) {
        if (H[i] > 0) {
            Ev[i] = Ev[i] / H[i];
            TTF[i] = TTF[i] / H[i];
        }
        else {
            Ev[i] = -ZPE;
            TTF[i] = -ZPE;
        }
        fprintf(outFile, "%.7lf,%lf,%.7lf,%.8lf\n", ZPE + i * Estep, exp(g[i] - g[0]), Ev[i] + ZPE, TTF[i]);
    }
    
}

// Returns 1 if derivative of energy with respect to any quantum number is <0
int ckderiv(double *freq, double *anh, unsigned int *nv, unsigned int n) {
    int ntest = 0;
    for (unsigned int k = 0; k < n && ntest == 0; k++) {
        double sum = 0;
        for (unsigned int i = 0; i < n; i++) {
            double inTerm = nv[i] + 0.5;
            sum += anh[k * n + i] * inTerm;
        }
        double knTerm = nv[k] + 0.5;
        double deriv = freq[k] + knTerm * anh[k * n + k] + sum;
        if (deriv < 0) {
            ntest = 1;
        }
    }
    return ntest;
}

// Calculates energy from frequencies, anharmonic constants, and quantum numbers
double nNRG(double *freq, double *anh, unsigned int *nv, unsigned int n) {
    double result = 0;
    for (unsigned int i = 0; i < n; i++) {
        result += freq[i] * (nv[i] + 0.5);
        for (unsigned int j = 0; j <= i; j++) {
            result += anh[j * n + i] * (nv[i] + 0.5) * (nv[j] + 0.5);
        }
    }
    return result;
}

// Calculates coupling terms from x_{nF} and quantum numbers
double nCoup(double *coup, unsigned int *nv, unsigned int n) {
    double result = 0;
    for (unsigned int i = 0; i < n; i++) {
        result += coup[i] * (nv[i] + 0.5);
    }
    return result;
}

// Calculates maximum allowed quantum number for a mode k given an energy Eu and a set of quantum numbers nv
int nvmax(unsigned int nv[], unsigned int k, double Eu, double anh[], double freq[], unsigned int n) {
    if (Eu < 0) {
        return -1;
    }
    if (anh[k * n + k] != 0) {
        double sum = 0;
        for (unsigned int i = 0; i < n; i++) {
            if (i != k) {
                sum += (nv[i] + 0.5) * anh[i * n + k];
            }
        }
        double anhkk = anh[k * n + k];
        double zpe = anhkk * 0.25 + freq[k] * 0.5 + sum * 0.5;
        double vd = (-freq[k] - sum) / (2 * anhkk) - 0.5;
        
        double anhksum = 0;
        for (unsigned int i=0; i < n; i++) {
            if (i != k) {
                anhksum = anh[k * n + i];
            }
        }
        
        double Dk = -(freq[k] + sum) * (freq[k] + sum) / 4 / anhkk - zpe;
        
        unsigned int vmax = vd * (1 - sqrt(1 - Eu / Dk));
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
    
    return -1;
}
