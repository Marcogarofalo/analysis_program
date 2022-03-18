
double ZAl(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double mu = file_head.musea;

    double num = in[j][1][t][0];
    double den = in[j][0][t + 1][0] - in[j][0][t][0];
    return (2 * mu * num / (den));

}


double ZVl(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double mu = file_head.musea;

    double num = in[j][4][t][0];
    double den = in[j][3][t + 1][0] - in[j][3][t][0];
    return (2 * mu * num / (den));

}