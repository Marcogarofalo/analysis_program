#define read_nissa_C
#include "read_nissa.hpp"

void read_single_nissa(double** out, std::string contraction, std::string gamma, std::string namefile) {
    std::fstream newfile;
    bool found_contraction = false;
    bool found_gamma = false;

    newfile.open(namefile, std::ios::in); // open a file to perform read operation using file object
    if (newfile.is_open()) { // checking whether the file is open
        std::string tp;
        while (getline(newfile, tp)) { // read data from file object and put it into string.
            // std::cout << "lin: " << tp << ":  len:" << tp.length() << "\n";
            // printf("is lenght enough %d \n", tp.length() >= 16);
            if (tp.length() >= 18) {
                if (tp.compare(1/*pos*/, 1/*leng*/, "#") == 0) {
                    // printf("is there # %d \n", tp.substr(1, 2).compare("#"));
                    // printf("is it a contraction %d \n", tp.compare(3,11,"Contraction")==0);
                    // printf("%s \n", tp.substr(3, 11).c_str());
                    if (tp.compare(3, 11, "Contraction") == 0) {

                        if (tp.compare(18, tp.size() - 18, contraction) == 0) {
                            // std::cout << "contraction found: " << tp << " \n";
                            // read inside contraction 
                            while (getline(newfile, tp)) {
                                // std::cout << "in contraction  " << tp << "\n";
                                if (tp.length() >= 3)  if (tp.compare(1/*pos*/, 1/*leng*/, "#") == 0) {
                                    if (tp.compare(3, tp.length() - 3, gamma) == 0) {
                                        // std::cout << "gamma found  " << tp << " "<< gamma<< "\n";
                                        getline(newfile, tp);
                                        int t = 0;
                                        while (!tp.empty()) {
                                            // std::cout << "in gamma " << gamma << ": t=" << t << "  " << tp << "\n";
                                            out[t][0] = std::stod(tp.substr(0, 22));
                                            out[t][1] = std::stod(tp.substr(24));
                                            t++;
                                            getline(newfile, tp);
                                        }
                                        printf("tmax=%d\n", t);
                                        return;
                                    }
                                }
                                if (tp.length() >= 18) if (tp.compare(1/*pos*/, 1/*leng*/, "#") == 0) if (tp.compare(18, tp.size() - 18, contraction)) {
                                    printf("error in %s:\n gamma structure not find in contraction\n", __func__);
                                    exit(1);
                                }

                            }
                        }
                    }
                }
            }
        }
        printf("error in %s:\n Contraction not find in file %s\n", __func__, namefile.c_str());
        exit(1);
    }
    else {
        printf("error in %s:\n unable to open %s\n", __func__, namefile.c_str());
        exit(1);
    }
    exit(1);
}


double*** return_all_nissa(std::string namefile, int Ncorr, int T, bool check ) {
    std::fstream newfile;
    double*** out;

    newfile.open(namefile, std::ios::in); // open a file to perform read operation using file object
    int Tc = 0;
    int Ncorrc = 0;
    if (newfile.is_open()) { // checking whether the file is open
        std::string tp;
        if (check) {
            while (getline(newfile, tp)) { // read data from file object and put it into string.  
                if (!tp.empty()) {
                    if (tp.compare(0, 1, "+") == 0 || tp.compare(0, 1, "-") == 0) {
                        Tc = 0;
                        while (!tp.empty()) {
                            Tc++;
                            getline(newfile, tp);
                        }
                        Ncorrc++;
                    }
                }
            }
            if (Ncorr != Ncorrc) printf("error in %s: given Ncorr=%d does nto match the read one=%d ", __func__, Ncorr, Ncorrc);
            if (T != Tc) printf("error in %s: given T=%d does nto match the read one=%d ", __func__, T, Tc);
        }

        printf("file %s: Ncorr=%d  T=%d\n", namefile.c_str(), Ncorr, T);
        out = malloc_3<double>(Ncorr, T, 2);
        newfile.clear();
        newfile.seekg(0);
        int corr = 0;
        while (getline(newfile, tp)) { // read data from file object and put it into string.
            if (!tp.empty()) {
                if (tp.compare(0, 1, "+") == 0 || tp.compare(0, 1, "-") == 0) {
                    int t = 0;
                    out[corr][t][0] = std::stod(tp.substr(0, 22));
                    out[corr][t][1] = std::stod(tp.substr(24));
                    while (!tp.empty()) {
                        // std::cout << "in gamma " << gamma << ": t=" << t << "  " << tp << "\n";
                        out[corr][t][0] = std::stod(tp.substr(0, 22));
                        out[corr][t][1] = std::stod(tp.substr(24));
                        t++;
                        getline(newfile, tp);
                    }
                    corr++;
                }
            }
        }
    }
    else {
        printf("error in %s:\n unable to open %s\n", __func__, namefile.c_str());
        exit(1);
    }
    return out;
}
