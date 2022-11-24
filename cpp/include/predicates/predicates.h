#pragma once

extern "C" {
void exactinit();
void initstaticfilter(double, double, double);

double estimate(const int elen, const double* e);
void invert_expansion(int elen, double* e);
void exact_scale(const int elen, double* e, const double b);

int fast_expansion_sum_zeroelim(const int elen, const double* e, const int flen, const double* f, double* h);
int fast_expansion_diff_zeroelim(const int elen, const double* e, const int flen, const double* f, double* h);
int scale_expansion_zeroelim(const int elen, const double* e, const double b, double* h);
int product_expansion_zeroelim(const int elen, const double* e, const int flen, const double* f, double* h);

double orient2d(double*, double*, double*);
double incircle(double* p1, double* p2, double* p3, double* p4);

double orient3d(const double* pa, const double* pb, const double* pc, const double* pd);
}

#define Quick_Two_Sum(a, b, x, y)                                                                                      \
    x = a + b;                                                                                                         \
    y = b - (x - a)
#define Two_Sum(a, b, x, y)                                                                                            \
    x = a + b;                                                                                                         \
    _bv = x - a;                                                                                                       \
    y = (a - (x - _bv)) + (b - _bv)
#define Two_One_Sum(a1, a0, b, x2, x1, x0)                                                                             \
    Two_Sum(a0, b, _i, x0);                                                                                            \
    Two_Sum(a1, _i, x2, x1)

#define Two_Diff(a, b, x, y)                                                                                           \
    x = a - b;                                                                                                         \
    _bv = a - x;                                                                                                       \
    y = (a - (x + _bv)) + (_bv - b)
#define Two_One_Diff(a1, a0, b, x2, x1, x0)                                                                            \
    Two_Diff(a0, b, _i, x0);                                                                                           \
    Two_Sum(a1, _i, x2, x1)

#define Split(a, _ah, _al)                                                                                             \
    _c = 1.3421772800000003e+008 * a;                                                                                  \
    _ah = _c - (_c - a);                                                                                               \
    _al = a - _ah
#define Two_Prod_PreSplit(a, b, _bh, _bl, x, y)                                                                        \
    x = a * b;                                                                                                         \
    Split(a, _ah, _al);                                                                                                \
    y = (_al * _bl) - (((x - (_ah * _bh)) - (_al * _bh)) - (_ah * _bl))
#define Two_Product_2Presplit(a, _ah, _al, b, _bh, _bl, x, y)                                                          \
    x = a * b;                                                                                                         \
    y = (_al * _bl) - (((x - _ah * _bh) - (_al * _bh)) - (_ah * _bl))

inline void two_sum(const double a, const double b, double* xy) {
    double _bv;
    Two_Sum(a, b, xy[1], xy[0]);
}

inline void two_diff(const double a, const double b, double* xy) {
    double _bv;
    Two_Diff(a, b, xy[1], xy[0]);
}
inline void two_prod(const double a, const double b, double* xy) {
    xy[1] = a * b;
    double _c, _ah, _al, _bh, _bl;
    Split(a, _ah, _al);
    Split(b, _bh, _bl);
    xy[0] = ((_ah * _bh - xy[1]) + _ah * _bl + _al * _bh) + _al * _bl;
}
inline void two_two_sum(const double* a, const double* b, double* x) {
    double _i, _j, _0, _bv;
    Two_One_Sum(a[1], a[0], b[0], _j, _0, x[0]);
    Two_One_Sum(_j, _0, b[1], x[3], x[2], x[1]);
}

inline void two_two_diff(const double* a, const double* b, double* x) {
    double _i, _j, _0, _bv, _u3;
    Two_One_Diff(a[1], a[0], b[0], _j, _0, x[0]);
    Two_One_Diff(_j, _0, b[1], _u3, x[2], x[1]);
    x[3] = _u3;
}

 