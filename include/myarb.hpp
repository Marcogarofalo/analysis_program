#ifndef MYARB
#define MYARB

#include "arb.h"
#include <stdexcept>

class myarb {
public:
    arb_t a;
    slong prec;
    myarb(slong prec_) {
        arb_init(a);
        prec = prec_;

    };
    myarb(arb_t b, slong prec_) {
        arb_init(a);
        arb_set(a, b);
        prec = prec_;
    };
    myarb(double d, slong prec_) {
        arb_init(a);
        arb_set_d(a, d);
        prec = prec_;
    };
    myarb(int i, slong prec_) {
        arb_init(a);
        arb_set_ui(a, i);
        prec = prec_;
    };

    ~myarb() {
        arb_clear(a);
    };


    // copy constructor
    myarb(const myarb& obj) {
        arb_init(a);
        arb_set(a, obj.a);
        prec = obj.prec;
    }
    // move constructor
    myarb(const myarb&& obj) {
        arb_init(a);
        arb_set(a, obj.a);
        prec = obj.prec;

    }
    myarb operator=(myarb const& obj) {
        arb_set(a, obj.a);
        prec = obj.prec;
        return *this;
    }
    // add
    myarb operator+(myarb const& obj) {
        myarb res(obj);
        arb_add(res.a, a, obj.a, prec);
        return res;
    }
    myarb operator+(int const& i) {
        myarb res(prec);
        arb_add_ui(res.a, a, i, prec);
        return res;
    }
    myarb operator+(double const& d) {
        myarb res(prec);
        arb_t ad;
        arb_init(ad);
        arb_set_d(ad, d);
        arb_add(res.a, a, ad, prec);
        arb_clear(ad);
        return res;
    }
    myarb operator+=(myarb const& obj) {
        arb_add(a, a, obj.a, prec);
        return *this;
    }
    myarb operator+=(int const& i) {
        arb_add_ui(a, a, i, prec);
        return *this;
    }
    myarb operator+=(double const& d) {
        arb_t ad;
        arb_init(ad);
        arb_set_d(ad, d);
        arb_add(a, a, ad, prec);
        arb_clear(ad);
        return *this;
    }
    // minus
    myarb operator-() const {
        myarb res(*this);
        arb_neg(res.a, this->a);
        return res;
    }
    myarb operator-(myarb const& obj) {
        myarb res(obj);
        arb_sub(res.a, a, obj.a, prec);
        return res;
    }
    myarb operator-(int const& i) {
        myarb res(prec);
        arb_sub_ui(res.a, a, i, prec);
        return res;
    }
    myarb operator-(double const& d) {
        myarb res(prec);
        arb_t ad;
        arb_init(ad);
        arb_set_d(ad, d);
        arb_sub(res.a, a, ad, prec);
        arb_clear(ad);
        return res;
    }
    myarb operator-=(myarb const& obj) {
        arb_sub(a, a, obj.a, prec);
        return *this;
    }
    myarb operator-=(int const& i) {
        arb_sub_ui(a, a, i, prec);
        return *this;
    }
    myarb operator-=(double const& d) {
        arb_t ad;
        arb_init(ad);
        arb_set_d(ad, d);
        arb_sub(a, a, ad, prec);
        arb_clear(ad);
        return *this;
    }
    // mult 
    myarb operator*(myarb const& obj) {
        myarb res(obj);
        arb_mul(res.a, a, obj.a, prec);
        return res;
    }
    myarb operator*(arb_t b) {
        myarb res(prec);
        arb_mul(res.a, a, b, prec);
        return res;
    }
    myarb operator*(int const& i) {
        myarb res(prec);
        arb_mul_ui(res.a, a, i, prec);
        return res;
    }
    myarb operator*(double const& d) {
        myarb res(prec);
        arb_t ad;
        arb_init(ad);
        arb_set_d(ad, d);
        arb_mul(res.a, a, ad, prec);
        arb_clear(ad);
        return res;
    }
    myarb operator*=(myarb const& obj) {
        arb_mul(a, a, obj.a, prec);
        return *this;
    }
    myarb operator*=(int const& i) {
        arb_mul_ui(a, a, i, prec);
        return *this;
    }
    myarb operator*=(double const& d) {
        arb_t ad;
        arb_init(ad);
        arb_set_d(ad, d);
        arb_mul(a, a, ad, prec);
        arb_clear(ad);
        return *this;
    }
    // div
    myarb operator/(myarb const& obj) {
        myarb res(obj);
        arb_div(res.a, a, obj.a, prec);
        return res;
    }
    myarb operator/(int const& i) {
        myarb res(prec);
        arb_div_ui(res.a, a, i, prec);
        return res;
    }
    myarb operator/(double const& d) {
        myarb res(prec);
        arb_t ad;
        arb_init(ad);
        arb_set_d(ad, d);
        arb_div(res.a, a, ad, prec);
        arb_clear(ad);
        return res;
    }
    myarb operator/=(myarb const& obj) {
        arb_div(a, a, obj.a, prec);
        return *this;
    }
    myarb operator/=(int const& i) {
        arb_div_ui(a, a, i, prec);
        return *this;
    }
    myarb operator/=(double const& d) {
        arb_t ad;
        arb_init(ad);
        arb_set_d(ad, d);
        arb_div(a, a, ad, prec);
        arb_clear(ad);
        return *this;
    }
};


class myacb {
public:
    slong prec;
    acb_t a;
    myacb(slong prec_) : prec{ prec_ } {
        acb_init(a);
    };
    myacb(const acb_t b, slong prec_) : prec{ prec_ } {
        acb_init(a);
        acb_set(a, b);
    };
    myacb(double d, slong prec_) : prec{ prec_ } {
        acb_init(a);
        acb_set_d(a, d);
    };
    myacb(int i, slong prec_) : prec{ prec_ } {
        acb_init(a);
        acb_set_ui(a, i);
    };

    ~myacb() {
        acb_clear(a);
    };


    // copy constructor
    myacb(const myacb& obj) {
        acb_init(a);
        acb_set(a, obj.a);
        prec = obj.prec;
    }
    // move constructor
    myacb(myacb&& obj) {
        // acb_init(a);
        // acb_set(a, obj.a);
        acb_swap(a, obj.a);
        prec = obj.prec;

    }
    // move assignment
    myacb& operator=(myacb&& obj) {
        // acb_set(a, obj.a);
        acb_swap(a, obj.a);
        prec = obj.prec;
        return *this;
    }
    // copy assignment
    myacb operator=(const myacb& obj) {
        acb_set(a, obj.a);
        prec = obj.prec;
        return *this;
    }
    // add
    myacb operator+(myacb const& obj) {
        myacb res(obj);
        acb_add(res.a, a, obj.a, prec);
        return res;
    }
    myacb operator+(int const& i) {
        myacb res(prec);
        acb_add_ui(res.a, a, i, prec);
        return res;
    }
    myacb operator+(double const& d) {
        myacb res(prec);
        acb_t ad;
        acb_init(ad);
        acb_set_d(ad, d);
        acb_add(res.a, a, ad, prec);
        acb_clear(ad);
        return res;
    }
    myacb operator+=(myacb const& obj) {
        acb_add(a, a, obj.a, prec);
        return *this;
    }
    myacb operator+=(int const& i) {
        acb_add_ui(a, a, i, prec);
        return *this;
    }
    myacb operator+=(double const& d) {
        acb_t ad;
        acb_init(ad);
        acb_set_d(ad, d);
        acb_add(a, a, ad, prec);
        acb_clear(ad);
        return *this;
    }
    // minus
    myacb operator-() const {
        myacb res(*this);
        acb_neg(res.a, this->a);
        return res;
    }
    myacb operator-(myacb const& obj) {
        myacb res(obj);
        acb_sub(res.a, a, obj.a, prec);
        return res;
    }
    myacb operator-(int const& i) {
        myacb res(prec);
        acb_sub_ui(res.a, a, i, prec);
        return res;
    }
    myacb operator-(double const& d) {
        myacb res(prec);
        acb_t ad;
        acb_init(ad);
        acb_set_d(ad, d);
        acb_sub(res.a, a, ad, prec);
        acb_clear(ad);
        return res;
    }
    myacb operator-=(myacb const& obj) {
        acb_sub(a, a, obj.a, prec);
        return *this;
    }
    myacb operator-=(int const& i) {
        acb_sub_ui(a, a, i, prec);
        return *this;
    }
    myacb operator-=(double const& d) {
        acb_t ad;
        acb_init(ad);
        acb_set_d(ad, d);
        acb_sub(a, a, ad, prec);
        acb_clear(ad);
        return *this;
    }
    // mult 
    myacb operator*(myacb const& obj) {
        myacb res(obj);
        acb_mul(res.a, a, obj.a, prec);
        return res;
    }
    myacb operator*(int const& i) {
        myacb res(prec);
        acb_mul_ui(res.a, a, i, prec);
        return res;
    }
    myacb operator*(double const& d) {
        myacb res(prec);
        acb_t ad;
        acb_init(ad);
        acb_set_d(ad, d);
        acb_mul(res.a, a, ad, prec);
        acb_clear(ad);
        return res;
    }
    myacb operator*=(myacb const& obj) {
        acb_mul(a, a, obj.a, prec);
        return *this;
    }
    myacb operator*=(int const& i) {
        acb_mul_ui(a, a, i, prec);
        return *this;
    }
    myacb operator*=(double const& d) {
        acb_t ad;
        acb_init(ad);
        acb_set_d(ad, d);
        acb_mul(a, a, ad, prec);
        acb_clear(ad);
        return *this;
    }
    // div
    myacb operator/(myacb const& obj) {
        myacb res(obj);
        acb_div(res.a, a, obj.a, prec);
        return res;
    }
    myacb operator/(int const& i) {
        myacb res(prec);
        acb_div_ui(res.a, a, i, prec);
        return res;
    }
    myacb operator/(double const& d) {
        myacb res(prec);
        acb_t ad;
        acb_init(ad);
        acb_set_d(ad, d);
        acb_div(res.a, a, ad, prec);
        acb_clear(ad);
        return res;
    }
    myacb operator/=(myacb const& obj) {
        acb_div(a, a, obj.a, prec);
        return *this;
    }
    myacb operator/=(int const& i) {
        acb_div_ui(a, a, i, prec);
        return *this;
    }
    myacb operator/=(double const& d) {
        acb_t ad;
        acb_init(ad);
        acb_set_d(ad, d);
        acb_div(a, a, ad, prec);
        acb_clear(ad);
        return *this;
    }
    // 
};

myacb exp(myacb b) {
    myacb res(b.prec);
    acb_exp(res.a, b.a, b.prec);
    return res;
}
myarb exp(myarb b) {
    myarb res(b.prec);
    arb_exp(res.a, b.a, b.prec);
    return res;
}

double arbtod(arb_t a) {
    double d, r;
    // arb_t ra;
    // arb_init(ra);
    try {
        d = std::stod(arb_get_str(a, 20, ARB_STR_NO_RADIUS));
        // arb_get_rad_arb(ra, a);
        // r = std::stod(arb_get_str(ra, 20, ARB_STR_NO_RADIUS));
        // if (r > fabs(d)) {
        //     printf("arbtod:\n radius is bigger than the error \n");
        //     arb_printn(a, 50/ 3.33, 0); flint_printf("\n");
        // };
    }
    catch (std::out_of_range) {
        char* c = arb_get_str(a, 20, ARB_STR_NO_RADIUS);
        // printf("arbtod out of range  %s\n",c);
        while (strcmp(c, "") != 0) {
            // printf("%s\n",c);
            if (*c == 'e') { c++; break; }
            else c++;
        }
        int i = std::stoi(c);
        if (i < 0)
            d = 0;
        else
            d = 1.0f / 0.0f;

    }
    return d;
}

#endif // !MYARB
