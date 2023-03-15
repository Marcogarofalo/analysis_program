#define read_nissa_C
#include "read_nissa.hpp"
#include "mutils.hpp"

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


double*** return_all_nissa(std::string namefile, int Ncorr, int T, bool check) {
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
            error(Ncorr != Ncorrc, 1, "return_all_nissa", " given Ncorr=%d does not match the read one=%d\n file:%s", Ncorr, Ncorrc, namefile.c_str());
            error(T != Tc, 1, "return_all_nissa", "error in %s: given T=%d does not match the read one=%d\n file:%s", T, Tc, namefile.c_str());
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

void read_all_nissa(double*** out, std::string namefile, int Ncorr, int T, bool check) {
    std::fstream newfile;
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
            error(Ncorr != Ncorrc, 1, "read_all_nissa", " given Ncorr=%d does not match the read one=%d\n file:%s", Ncorr, Ncorrc, namefile.c_str());
            error(T != Tc, 1, "read_all_nissa", "error in %s: given T=%d does not match the read one=%d\n file:%s", T, Tc, namefile.c_str());
        }

        printf("file %s: Ncorr=%d  T=%d\n", namefile.c_str(), Ncorr, T);
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
    
    return;
}

void read_all_nissa_gamma(double*** out, std::string namefile, int Ncorr, int T, std::vector<int> id_gamma, bool check) {
    std::fstream newfile;
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
            error(Ncorr != Ncorrc, 1, "read_all_nissa", " given Ncorr=%d does not match the read one=%d\n file:%s", Ncorr, Ncorrc, namefile.c_str());
            error(T != Tc, 1, "read_all_nissa", "error in %s: given T=%d does not match the read one=%d\n file:%s", T, Tc, namefile.c_str());
        }

        printf("file %s: Ncorr=%d  T=%d\n", namefile.c_str(), Ncorr, T);
        newfile.clear();
        newfile.seekg(0);
        int corr = 0;
        int id=0;
        while (getline(newfile, tp)) { // read data from file object and put it into string.
            if (!tp.empty()) {
                if (tp.compare(0, 1, "+") == 0 || tp.compare(0, 1, "-") == 0) {
                    if (corr==id_gamma[id]){
                        id++;
                        int t = 0;
                        out[id][t][0] = std::stod(tp.substr(0, 22));
                        out[id][t][1] = std::stod(tp.substr(24));
                        while (!tp.empty()) {
                            // std::cout << "in gamma " << gamma << ": t=" << t << "  " << tp << "\n";
                            out[id][t][0] = std::stod(tp.substr(0, 22));
                            out[id][t][1] = std::stod(tp.substr(24));
                            t++;
                            getline(newfile, tp);
                        }
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
    
    return;
}

struct_nissa_out_info::struct_nissa_out_info(std::string namefile) {
    std::fstream newfile;

    newfile.open(namefile, std::ios::in); // open a file to perform read operation using file object
    T = 0;
    Ncorr = 0;
    Ncontr = 0;
    Ngamma = 0;
    int tmp_gamma;
    std::string tp;
    if (newfile.is_open()) { // checking whether the file is open
        while (getline(newfile, tp)) {
            if (tp.length() >= 18) {
                if (tp.compare(1/*pos*/, 1/*leng*/, "#") == 0) {
                    if (tp.compare(3, 11, "Contraction") == 0) {
                        Ncontr++;
                        tmp_gamma = 0;
                        while (getline(newfile, tp)) {
                            if (tp.length() >= 3) {
                                if (tp.compare(3, 11, "Contraction") == 0) {
                                    if (Ncontr == 1) Ngamma = tmp_gamma;
                                    else if (Ngamma != tmp_gamma) {
                                        printf("error in %s:\n", __func__);
                                        printf("number of gamma in first contraction %d\n", Ngamma);
                                        printf("number of gamma in %d th contraction %d\n", tmp_gamma, Ncontr);
                                        exit(1);
                                    }
                                    Ncontr++;
                                    tmp_gamma = 0;
                                }
                                else if (tp.compare(1/*pos*/, 1/*leng*/, "#") == 0) {
                                    if (Ncontr == 1) {
                                        std::vector<std::string> x = split(tp, ' ');
                                        gamma.emplace_back(x[1]);
                                    }
                                    getline(newfile, tp);
                                    int t = 0;
                                    while (!tp.empty()) {
                                        t++;
                                        getline(newfile, tp);
                                    }

                                    tmp_gamma++;
                                    Ncorr++;
                                    if (Ncontr == 1 && tmp_gamma == 1) {
                                        T = t;
                                    }
                                    else if (T != t) {
                                        printf("error in %s:\n", __func__);
                                        printf("T in first contraction %d\n", T);
                                        printf("T in %d th correlator %d\n", Ncorr, t);
                                        exit(1);
                                    }
                                }
                            }
                        }
                        if (Ncontr == 1) Ngamma = tmp_gamma;
                        else if (Ngamma != tmp_gamma) {
                            printf("error in %s:\n", __func__);
                            printf("number of gamma in first contraction %d\n", Ngamma);
                            printf("number of gamma in %d th contraction %d\n", tmp_gamma, Ncontr);
                            exit(1);
                        }

                    }
                }
            }
        }

    }
    else {
        printf("error in %s:\n unable to open %s\n", __func__, namefile.c_str());
        exit(1);
    }
}

void struct_nissa_out_info::print() {
    printf("nissa mes file:\n");
    printf("T %i\n", T);
    printf("Ncorr %i\n", Ncorr);
    printf("Ncontr %i\n", Ncontr);
    printf("Ngamma %i\n", Ngamma);
    for (std::string s : gamma)
        printf("%s\t", s.c_str());
    printf("\n");
}

std::vector<int> struct_nissa_out_info::inidices_of_gamma(std::vector<std::string> mygammas) {
    std::vector<int> out(Ncontr * mygammas.size());
    std::vector<int> ii(mygammas.size());
    for (int mg = 0; mg < mygammas.size();mg++) {
        int i = -1;
        for (i = 0; i < Ngamma; i++) {
            if (gamma[i].compare(mygammas[mg]) == 0)
                break;
        }
        error(i == Ngamma, 1, "inidices_of_gamma:", "%s gamma not found", mygammas[mg].c_str());
        error(i == -1, 1, "inidices_of_gamma:", "%s gamma not found", mygammas[mg].c_str());
        ii[mg] = i;
    }
    for (int i = 0; i < Ncontr; i++) {
        for (int mg = 0; mg < mygammas.size();mg++) {
            out[mg + i * mygammas.size()] = ii[mg] + i * Ngamma;
        }
    }
    return out;
}