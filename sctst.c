// Calculate cumulative reaction probability from density of vibrational states and their coupling to the reaction coordinate, using either the Miller (M) potential barrier form or Wagner (W) form
// Usage: ./sctst.x [SCTST input file] [States file] [M or W]

#include <stdio.h>
#include <math.h>
#define MaxBins 5000
#define CRPLen 5000

double wagProb(double barHeight, double revBarHeight, double xFF, double omegaF, double energy);
double thetaIrreg(double zlo, double zhi, double eta, double rho, double aa, double bb, double cc);
double thetaReg(int iseg, double zlo, double zhi, double eta, double rho);

int main(int argc, const char * argv[]) {
    FILE *inputF = fopen(argv[1], "r");
    
    char *buf;
    size_t bufLen;
    getline(&buf, &bufLen, inputF);
    
    // Read in imaginary frequency (a.u.)
    double rxnFreq;
    fscanf(inputF, "%lf", &rxnFreq);
    
    fscanf(inputF, "%s", buf);
    getline(&buf, &bufLen, inputF);
    
    // Read in rxn barrier height w/ reac ZPE, reacG0, TS G0 (a.u.)
    double barHeight;
    fscanf(inputF, "%lf", &barHeight);
    
    fscanf(inputF, "%s", buf);
    getline(&buf, &bufLen, inputF);
    
    // Read in reverse barrier height w/ prod ZPE, TS G0 (a.u.)
    double revBarHeight;
    fscanf(inputF, "%lf", &revBarHeight);
    
    fscanf(inputF, "%s", buf);
    getline(&buf, &bufLen, inputF);
    
    // Read in xFF
    double xFF;
    fscanf(inputF, "%lf", &xFF);
    
    fscanf(inputF, "%s", buf);
    getline(&buf, &bufLen, inputF);
    
    // Read in ZPE of reactants (a.u.)
    double reacZPE;
    fscanf(inputF, "%lf", &reacZPE);
    
    fscanf(inputF, "%s", buf);
    getline(&buf, &bufLen, inputF);
    
    // Read in ZPE of products (a.u.)
    double prodZPE;
    fscanf(inputF, "%lf", &prodZPE);
    
    FILE *statesF = fopen(argv[2], "r");
    
    size_t len;
    getline(&buf, &len, statesF);
    getline(&buf, &len, statesF);
    
    double energies[MaxBins];
    double nStates[MaxBins];
    double Ev[MaxBins];
    double coupling[MaxBins];
    unsigned int nBins = 0;
    
    int retVal = 4;
    while (retVal == 4) {
        retVal = fscanf(statesF, "%lf,%lf,%lf,%lf", &energies[nBins], &nStates[nBins], &Ev[nBins], &coupling[nBins]);
        nBins++;
    }
    nBins--;
    
    
    
    double Estep = (energies[nBins - 1] + (energies[1] - energies[0])) / CRPLen;
    double CRP[CRPLen];
    for (unsigned int i = 0; i < CRPLen; i++) {
        CRP[0] = 0;
    }
    
    double deltaE;
    double tmp;
    double omegaF;
    double prob;
    
    double ergicity = barHeight - revBarHeight;
    int ergjMin = ergicity / Estep;
    double wagForward; // Wagner requires rxn to be exoergic
    double wagReverse;
    if (ergjMin < 0) {
        ergjMin = 0;
        wagForward = barHeight;
        wagReverse = revBarHeight;
    }
    else {
        wagForward = revBarHeight;
        wagReverse = barHeight;
    }
    
    
    for (unsigned int i = 0; i < nBins; i++) {
        omegaF = rxnFreq + coupling[i];
        int jMin = (energies[i] - reacZPE) / Estep;
        if (jMin < ergjMin) {
            jMin = ergjMin;
        }
        
        for (unsigned int j = jMin; j < CRPLen; j++) {
            if (argv[3][0] == 'M') {
                deltaE = barHeight + energies[i] - j * Estep;
                tmp = omegaF * omegaF + 4 * xFF * deltaE;
                if (tmp >= 0) {
                    prob = 1 / (1 + exp(2 * M_PI * (-omegaF + sqrt(tmp)) / (2 * xFF)));
                }
                else {
                    prob = 0;
                }
            }
            else if (argv[3][0] == 'W') {
                prob = wagProb(wagForward + energies[i], wagReverse + energies[i], xFF, omegaF, j * Estep);
            }
            else {
                printf("Please specify W or M for the method.\n");
                return 0;
            }
            
            CRP[j] += nStates[i] * prob;
        }
        
    }
    
    char outName[15];
    sprintf(outName, "CRP_%c.txt", argv[3][0]);
    
    FILE *outFile = fopen(outName, "w");
    for (unsigned int i = 0; i < CRPLen; i++) {
        fprintf(outFile, "%lf,%.15lf\n", i * Estep, CRP[i]);
    }
    
    return 0;
}

double thetaReg(int iseg, double zlo, double zhi, double eta, double rho) {
    double c0 = rho * rho + eta - 1;
    double c3 = 1 + rho;
    double c2 = sqrt(c0);
    double c1 = sqrt(eta);
    
    if (eta > 1) {
        return M_PI * (c3 - c2 - c1);
    }
    else {
        double c4 = sqrt(1 - eta);
        double c5 = c3 * c4;
        double c6 = c3 - eta;
        double thtal, thtah;
        if (iseg >= 2) {
            double arg1l = (c6 * zlo - eta) / zlo / c5;
            arg1l = arg1l < -1 ? -1 : arg1l;
            arg1l = arg1l > 1 ? 1 : arg1l;
            double arg2l = (c6 - c0 * zlo) / c5;
            arg2l = arg2l < -1 ? -1 :arg2l;
            arg2l = arg2l > 1 ? 1 : arg2l;
            double arg3l = (rho * zlo - 1) / (1 + zlo) / c4;
            arg3l = arg3l < -1 ? -1 : arg3l;
            arg3l = arg3l > 1 ? 1 : arg3l;
            thtal = c3 * asin(arg3l) + c2 * asin(arg2l) - c1 * asin(arg1l);
        }
        else {
            thtal = -M_PI / 2 * (c3 - c2 - c1);
        }
        if (iseg <= 2) {
            double arg1h = (c6 * zhi - eta) / zhi / c5;
            arg1h = arg1h < -1 ? -1 : arg1h;
            arg1h = arg1h > 1 ? 1 : arg1h;
            double arg2h = (c6 - c0 * zhi) / c5;
            arg2h = arg2h < -1 ? -1 : arg2h;
            arg2h = arg2h > 1 ? 1 : arg2h;
            double arg3h = (rho * zhi - 1) / (1 + zhi) / c4;
            arg3h = arg3h < -1 ? -1 : arg3h;
            arg3h = arg3h > 1 ? 1 : arg3h;
            thtah = c3 * asin(arg3h) + c2 * asin(arg2h) - c1 * asin(arg1h);
        }
        else {
            thtah = M_PI / 2 * (c3 - c2 - c1);
        }
        return thtah - thtal;
    }
}

double thetaIrreg(double zlo, double zhi, double eta, double rho, double aa, double bb, double cc) {
    double c1 = sqrt(1 - eta);
    double c2 = 1 - rho * rho - eta;
    double c3 = 1 + rho - eta;
    double arg3h = (rho * zhi - 1) / (1 + zhi) / c1;
    arg3h = arg3h < -1 ? -1 : arg3h;
    arg3h = arg3h > 1 ? 1 : arg3h;
    double arg3l = (rho * zlo - 1) / (1 + zlo) / c1;
    arg3l = arg3l < -1 ? -1 : arg3l;
    arg3l = arg3l > 1 ? 1 : arg3l;
    double thta = (1 + rho) * (asin(arg3h) - asin(arg3l));
    
    double rhi = cc + bb * zhi + aa * zhi * zhi;
    rhi = rhi < 0 ? 0 : rhi;
    double rlo  = cc + bb * zlo + aa * zlo * zlo;
    rlo = rlo < 0 ? 0 : rlo;
    double arg1 = (2 * cc + bb * zhi + 2 * sqrt(cc * rhi)) / zhi;
    double arg2 = (2 * cc + bb * zlo + 2 * sqrt(cc * rlo)) / zlo;
    thta -= sqrt(-eta) * log(arg1 / arg2);
    
    if (-c2 >= 0) {
        double arg2h = (c2 * zhi + c3) / (1 + rho) / c1;
        arg2h = arg2h < -1 ? -1 : arg2h;
        arg2h = arg2h > 1 ? 1 : arg2h;
        double arg2l = (c2 * zlo + c3) / (1 + rho) / c1;
        arg2l = arg2l < -1 ? -1 : arg2l;
        arg2l = arg2l > 1 ? 1 : arg2l;
        return thta + sqrt(-c2) * (asin(arg2h) - asin(arg2l));
    }
    else {
        double arg1 = 2 * sqrt(aa * rhi) + 2 * aa * zhi + bb;
        double arg2 = 2 * sqrt(aa * rlo) + 2 * aa * zlo + bb;
        return thta + sqrt(c2) * log(arg1 / arg2);
    }
}

double wagProb(double barHeight, double revBarHeight, double xFF, double omegaF, double energy) {
    double rho = sqrt(revBarHeight / barHeight);
    
    double Db = -omegaF * omegaF / 4.0 / xFF * (1.0 - 1.0 / rho + 1.0 / rho / rho);
    double alphab = fabs(omegaF) * (1.0 + rho) / rho * sqrt(0.5 / Db);
    double f = sqrt(Db / barHeight);
    double ur = (sqrt(rho * rho + rho + 1) - (rho - 1)) / (f > 1 ? 3 * f : 3);
    double Rr = (210 * rho * rho + 168 * (1 + f) * rho * (1 - rho) * ur + 140 * (f - rho * (1 + f) * (1 + f) + f * rho * rho) * ur * ur - 120 * f * (1 + f) * (1 - rho) * ur * ur * ur + 105 * f * f * ur * ur * ur * ur) / (210 * rho * rho + 336 * rho * (1 - rho) * f * ur + 140 * (1 - 4 * rho + rho * rho) * f * f * ur * ur - 240 * (1 - rho) * f * f * f * ur * ur * ur + 105 * pow(f * ur, 4));
    double alphar = alphab * f * Rr;
    double umr1 = (sqrt((1 - rho) * (1 - rho) * pow(1 - Rr * f, 2) - 4 * rho * (1 - Rr * f * f) * (Rr - 1)) + (1 - rho) * (1 - Rr * f)) / (2 * (1 - Rr * f * f));
    double umr2 = -(sqrt(pow(1 - rho, 2) * pow(1 - Rr * f, 2) - 4 * rho * (1 - Rr * f * f) * (Rr - 1)) + (1 - rho) * (1 - Rr * f)) / (2 * (1 - Rr * f * f));
    double umr = (umr1 < 0 || umr1 > ur) ? umr2 : umr1;
    double ybr = log((1 - umr) / (rho + umr)) / alphab;
    double dr = -ybr + log((1 - f * umr) / (rho + f * umr)) / alphar;
    double up = (sqrt(rho * rho + rho + 1) + rho - 1) / (f > 1 ? 3 * f : 3);
    double Rp = (210 * rho * rho + 168 * (1 + f) * rho * (1 - rho) * up + 140 * (f - rho * (1 + f) * (1 + f) + f * rho * rho) * up * up - 120 * f * (1 + f) * (1 - rho) * up * up * up + 105 * f * f * pow(up, 4)) / (210 * rho * rho + 336 * rho * (1 - rho) * f * up + 140 * (1 - 4 * rho + rho * rho) * f * f * up * up - 240 * (1 - rho) * f * f * f * up * up * up + 105 * pow(f * up, 4));
    double alphap = alphab * f * Rp;
    double ump1 = (sqrt((1 - rho) * (1 - rho) * pow(1 - Rp * f, 2) - 4 * rho * (1 - Rp * f * f) * (Rp - 1)) + (1 - rho) * (1 - Rp * f)) / (2 * (1 - Rp * f * f));
    double ump2 = -(sqrt(pow(1 - rho, 2) * pow(1 - Rp * f, 2) - 4 * rho * (1 - Rp * f * f) * (Rp - 1)) + (1 - rho) * (1 - Rp * f)) / (2 * (1 - Rp * f * f));
    double ump = (ump1 < 0 || ump1 > up) ? ump2 : ump1;
    double ybp = log((-1 - ump) / (-rho + ump)) / alphab;
    double dp = log((-1 - f * ump) / (-rho + f * ump)) / alphap - ybp;
    double zrb = exp(alphar * (ybr + dr));
    double zbr = exp(alphab * ybr);
    double zpb = exp(alphap * (ybp + dp));
    double zbp = exp(alphab * ybp);
    double headb = 2 * Db / omegaF * rho / (1 + rho);
    double headr = headb * barHeight / Db / Rr;
    double headp = headr * Rr / Rp;
    
    double thtaEn = energy < barHeight ? energy : 2 * barHeight - energy;
    double ep = thtaEn - barHeight + Db;
    double eta = ep / Db;
    double thta;
    
    double zlo;
    double zhi;
    
    if (thtaEn > (barHeight - Db)) {
        double sqrtdelld2a = (1 + rho) * sqrt(1 - eta) / (1 - rho * rho - eta);
        double zzlo = -(1 + rho - eta) / (1 - rho * rho - eta) + sqrtdelld2a;
        double zzhi = zzlo - 2 * sqrtdelld2a;
        zlo = zbr > zzlo ? zbr : zzlo;
        zhi = zbp < zzhi ? zbp : zzhi;
        thta = headb * thetaReg(2, zlo, zhi, eta, rho);
    }
    else {
        zlo = zbr;
        zhi = zbp;
        double bb = 2 * Db * (1 + rho - eta);
        double aa = Db * (1 - rho * rho - eta);
        double cc = -ep;
        thta = headb * thetaIrreg(zlo, zhi, eta, rho, aa, bb, cc);
    }
    if (zlo == zbr) {
        eta = thtaEn / barHeight;
        double sqrtdelld2a = (1 + rho) * sqrt(1 - eta) / (1 - rho * rho - eta);
        zlo = -(1 + rho - eta) / (1 - rho * rho - eta) + sqrtdelld2a;
        double zzhi = zlo - 2 * sqrtdelld2a;
        thta += headr * thetaReg(1, zlo, zrb, eta, rho);
        if (zhi == zbp) {
            thta += headp * thetaReg(3, zpb, zzhi, eta, rho);
        }
    }
    
    if (energy < barHeight) {
        return 1 / (1 + exp(2 * thta));
    }
    else if ((energy - barHeight) < barHeight) {
        return 1 - 1 / (1 + exp(2 * thta));
    }
    else {
        return 1;
    }
}
