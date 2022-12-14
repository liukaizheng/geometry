#include <predicates/generic_point.h>
#include <predicates/predicates.h>

#include <cmath>
#include <vector>

inline int orient2d_EEE(const explicitPoint2D& p1, const explicitPoint2D& p2, const explicitPoint2D& p3) {
    const double ret = orient2d(p1.ptr(), p2.ptr(), p3.ptr());
    return (ret > 0) - (ret < 0);
}

inline bool same_point(const double* p, const double* q) { return p[0] == q[0] && p[1] == q[1] && p[2] == q[2]; }

inline int orient2d_SEE_filtered(const implicitPoint2D_SSI& p1, const double* p2, const double* p3) {
    double max_var = 0;
    if (!p1.getFilteredLambda(max_var)) {
        return IPSign::ZERO;
    }

    double t1x = p2[1] - p3[1];
    double t1y = p3[0] - p2[0];
    double e2 = p1.ssfilter[0] * t1x;
    double e3 = p1.ssfilter[1] * t1y;
    double e = e2 + e3;
    double pr1 = p2[0] * p3[1];
    double pr2 = p2[1] * p3[0];
    double pr = pr1 - pr2;
    double dpr = p1.ssfilter[2] * pr;
    double det = dpr + e;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(p2[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p2[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(t1x)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(t1y)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= 1.110439865059655e-14;
    if (det > epsilon) return IPSign::POSITIVE;
    if (-det > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient2d_SEE_interval(const implicitPoint2D_SSI& p1, const double* p2, const double* p3) {

    if (!p1.getIntervalLambda()) return IPSign::ZERO;

    IntervalNumber t1x(p2[1] - p3[1]);
    IntervalNumber t1y(p3[0] - p2[0]);
    IntervalNumber e2(p1.dfilter[0] * t1x);
    IntervalNumber e3(p1.dfilter[1] * t1y);
    IntervalNumber e(e2 + e3);
    IntervalNumber pr1(p2[0] * p3[1]);
    IntervalNumber pr2(p2[1] * p3[0]);
    IntervalNumber pr(pr1 - pr2);
    IntervalNumber dpr(p1.dfilter[2] * pr);
    IntervalNumber det(dpr + e);

    if (!det.signIsReliable()) return IPSign::ZERO;
    return det.sign();
}

inline int orient2d_SEE_exact(const implicitPoint2D_SSI& p1, const double* p2, const double* p3) {

    double l1x[32], l1y[32], d1[16];
    int l1x_len, l1y_len, d1_len;
    p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, d1, d1_len);
    if (d1[d1_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    double t1x[2];
    two_diff(p2[1], p3[1], t1x);
    double t1y[2];
    two_diff(p3[0], p2[0], t1y);
    double e2[128];
    int e2_len = product_expansion_zeroelim(l1x_len, l1x, 2, t1x, e2);
    double e3[128];
    int e3_len = product_expansion_zeroelim(l1y_len, l1y, 2, t1y, e3);
    std::vector<double> e(static_cast<uint32_t>((e2_len * e3_len) << 1));
    int e_len = product_expansion_zeroelim(e2_len, e2, e3_len, e3, e.data());
    double pr1[2];
    two_prod(p2[0], p3[1], pr1);
    double pr2[2];
    two_prod(p2[1], p3[0], pr2);
    double pr[4];
    two_two_diff(pr1, pr2, pr);
    double dpr[128];
    int dpr_len = product_expansion_zeroelim(d1_len, d1, 4, pr, dpr);
    std::vector<double> det(static_cast<uint32_t>((dpr_len * e_len) << 1));
    uint32_t det_len = static_cast<uint32_t>(product_expansion_zeroelim(dpr_len, dpr, e_len, e.data(), det.data()));

    if (det[det_len - 1] > 0.0) {
        return IPSign::POSITIVE;
    } else if (det[det_len - 1] < 0.0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient2d_SEE(const implicitPoint2D_SSI& p1, const double* p2, const double* p3) {
    int ret;
    if ((ret = orient2d_SEE_filtered(p1, p2, p3)) != 0) {
        return ret;
    }
    if ((ret = orient2d_SEE_interval(p1, p2, p3)) != 0) {
        return ret;
    }
    return orient2d_SEE_exact(p1, p2, p3);
}

inline int orient2d_SSE_filtered(const implicitPoint2D_SSI& p1, const implicitPoint2D_SSI& p2, const double* p3) {
    double max_var = 0;
    if (!p1.getFilteredLambda(max_var) || !p2.getFilteredLambda(max_var)) return IPSign::ZERO;

    double a = p1.ssfilter[2] * p2.ssfilter[0];
    double b = p2.ssfilter[2] * p1.ssfilter[0];
    double c = p1.ssfilter[2] * p3[1];
    double e = p1.ssfilter[2] * p2.ssfilter[1];
    double f = p2.ssfilter[2] * p1.ssfilter[1];
    double g = p1.ssfilter[2] * p3[0];
    double ab = a - b;
    double cd = c - p1.ssfilter[1];
    double ef = e - f;
    double gh = g - p1.ssfilter[0];
    double abcd = ab * cd;
    double efgh = ef * gh;
    double L = abcd - efgh;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(p3[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3[1])) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= 3.837902218251096e-13;
    if (L > epsilon) return IPSign::POSITIVE;
    if (-L > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient2d_SSE_interval(const implicitPoint2D_SSI& p1, const implicitPoint2D_SSI& p2, const double* p3) {
    if (!p1.getIntervalLambda() || !p2.getIntervalLambda()) return IPSign::ZERO;

    IntervalNumber a(p1.dfilter[2] * p2.dfilter[0]);
    IntervalNumber b(p2.dfilter[2] * p1.dfilter[0]);
    IntervalNumber c(p1.dfilter[2] * p3[1]);
    IntervalNumber e(p1.dfilter[2] * p2.dfilter[1]);
    IntervalNumber f(p2.dfilter[2] * p1.dfilter[1]);
    IntervalNumber g(p1.dfilter[2] * p3[0]);
    IntervalNumber ab(a - b);
    IntervalNumber cd(c - p1.dfilter[1]);
    IntervalNumber ef(e - f);
    IntervalNumber gh(g - p1.dfilter[0]);
    IntervalNumber abcd(ab * cd);
    IntervalNumber efgh(ef * gh);
    IntervalNumber L(abcd - efgh);

    if (!L.signIsReliable()) return IPSign::ZERO;
    return L.sign();
}

inline int orient2d_SSE_exact(const implicitPoint2D_SSI& p1, const implicitPoint2D_SSI& p2, const double* p3) {

    double l1x[32], l1y[32], d1[16], l2x[32], l2y[32], d2[16];
    int l1x_len, l1y_len, d1_len, l2x_len, l2y_len, d2_len;
    p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, d1, d1_len);
    p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, d2, d2_len);
    if (d1[d1_len - 1] == 0.0 || d2[d2_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> a(static_cast<uint32_t>((d1_len * l2x_len) << 1));
    int a_len = product_expansion_zeroelim(d1_len, d1, l2x_len, l2x, a.data());
    std::vector<double> b(static_cast<uint32_t>((d2_len * l1x_len) << 1));
    int b_len = product_expansion_zeroelim(d2_len, d2, l1x_len, l1x, b.data());
    double c[32];
    int c_len = scale_expansion_zeroelim(d1_len, d1, p3[1], c);
    std::vector<double> e(static_cast<uint32_t>((d1_len * l2y_len) << 1));
    int e_len = product_expansion_zeroelim(d1_len, d1, l2y_len, l2y, e.data());
    std::vector<double> f(static_cast<uint32_t>((d2_len * l1y_len) << 1));
    int f_len = product_expansion_zeroelim(d2_len, d2, l1y_len, l1y, f.data());
    double g[32];
    int g_len = scale_expansion_zeroelim(d1_len, d1, p3[0], g);
    std::vector<double> ab(static_cast<uint32_t>((a_len * b_len) << 1));
    int ab_len = fast_expansion_diff_zeroelim(a_len, a.data(), b_len, b.data(), ab.data());
    double cd[64];
    int cd_len = fast_expansion_diff_zeroelim(c_len, c, l1y_len, l1y, cd);
    std::vector<double> ef(static_cast<uint32_t>((e_len * f_len) << 1));
    int ef_len = fast_expansion_diff_zeroelim(e_len, e.data(), f_len, f.data(), ef.data());
    double gh[64];
    int gh_len = fast_expansion_diff_zeroelim(g_len, g, l1x_len, l1x, gh);
    std::vector<double> abcd(static_cast<uint32_t>((ab_len * cd_len) << 1));
    int abcd_len = product_expansion_zeroelim(ab_len, ab.data(), cd_len, cd, abcd.data());
    std::vector<double> efgh(static_cast<uint32_t>((ef_len * gh_len) << 1));
    int efgh_len = product_expansion_zeroelim(ef_len, ef.data(), gh_len, gh, efgh.data());
    std::vector<double> L(static_cast<uint32_t>((abcd_len * efgh_len) << 1));
    const uint32_t L_len =
        static_cast<uint32_t>(fast_expansion_diff_zeroelim(abcd_len, abcd.data(), efgh_len, efgh.data(), L.data()));
    if (L[L_len - 1] > 0.0) {
        return IPSign::POSITIVE;
    } else if (L[L_len - 1] < 0.0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient2d_SSE(const implicitPoint2D_SSI& p1, const implicitPoint2D_SSI& p2, const double* p3) {
    int ret = 0;
    if ((ret = orient2d_SSE_filtered(p1, p2, p3)) != 0) {
        return ret;
    }
    if ((ret = orient2d_SSE_interval(p1, p2, p3)) != 0) {
        return ret;
    }
    return orient2d_SSE_exact(p1, p2, p3);
}

inline int
orient2d_SSS_filtered(const implicitPoint2D_SSI& p1, const implicitPoint2D_SSI& p2, const implicitPoint2D_SSI& p3) {
    double max_var = 0.0;
    if (!p1.getFilteredLambda(max_var) || !p2.getFilteredLambda(max_var) || !p3.getFilteredLambda(max_var)) {
        return IPSign::ZERO;
    }

    const double* sf1 = p1.ssfilter.data();
    const double* sf2 = p2.ssfilter.data();
    const double* sf3 = p3.ssfilter.data();

    double a = sf1[2] * sf2[0];
    double b = sf2[3] * sf1[0];
    double c = sf1[2] * sf3[1];
    double d = sf3[3] * sf1[1];
    double e = sf1[2] * sf2[1];
    double f = sf2[3] * sf1[1];
    double g = sf1[2] * sf3[0];
    double h = sf3[3] * sf1[0];
    double ab = a - b;
    double cd = c - d;
    double ef = e - f;
    double gh = g - h;
    double abcd = ab * cd;
    double efgh = ef * gh;
    double L = abcd - efgh;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 1.364575119566784e-12;
    if (L > epsilon) return IPSign::POSITIVE;
    if (-L > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int
orient2d_SSS_interval(const implicitPoint2D_SSI& p1, const implicitPoint2D_SSI& p2, const implicitPoint2D_SSI& p3) {
    if (!p1.getIntervalLambda() || !p2.getIntervalLambda() || !p3.getIntervalLambda()) return IPSign::ZERO;

    IntervalNumber a(p1.dfilter[2] * p2.dfilter[0]);
    IntervalNumber b(p2.dfilter[2] * p1.dfilter[0]);
    IntervalNumber c(p1.dfilter[2] * p3.dfilter[1]);
    IntervalNumber d(p3.dfilter[2] * p1.dfilter[1]);
    IntervalNumber e(p1.dfilter[2] * p2.dfilter[1]);
    IntervalNumber f(p2.dfilter[2] * p1.dfilter[1]);
    IntervalNumber g(p1.dfilter[2] * p3.dfilter[0]);
    IntervalNumber h(p3.dfilter[2] * p1.dfilter[0]);
    IntervalNumber ab(a - b);
    IntervalNumber cd(c - d);
    IntervalNumber ef(e - f);
    IntervalNumber gh(g - h);
    IntervalNumber abcd(ab * cd);
    IntervalNumber efgh(ef * gh);
    IntervalNumber L(abcd - efgh);

    if (!L.signIsReliable()) return IPSign::ZERO;
    return L.sign();
}

inline int
orien2d_SSS_exact(const implicitPoint2D_SSI& p1, const implicitPoint2D_SSI& p2, const implicitPoint2D_SSI& p3) {
    double l1x[32], l1y[32], d1[16], l2x[32], l2y[32], d2[16], l3x[32], l3y[32], d3[16];
    int l1x_len, l1y_len, d1_len, l2x_len, l2y_len, d2_len, l3x_len, l3y_len, d3_len;
    p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, d1, d1_len);
    p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, d2, d2_len);
    p3.getExactLambda(l3x, l3x_len, l3y, l3y_len, d3, d3_len);

    if (d1[d1_len - 1] == 0.0 || d2[d2_len - 1] == 0.0 || d3[d3_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> a(static_cast<uint32_t>((d1_len * l2x_len) << 1));
    int a_len = product_expansion_zeroelim(d1_len, d1, l2x_len, l2x, a.data());
    std::vector<double> b(static_cast<uint32_t>((d2_len * l1x_len) << 1));
    int b_len = product_expansion_zeroelim(d2_len, d2, l1x_len, l1x, b.data());
    std::vector<double> c(static_cast<uint32_t>((d1_len * l3y_len) << 1));
    int c_len = product_expansion_zeroelim(d1_len, d1, l3y_len, l3y, c.data());
    std::vector<double> d(static_cast<uint32_t>((d3_len * l1y_len) << 1));
    int d_len = product_expansion_zeroelim(d3_len, d3, l1y_len, l1y, d.data());
    std::vector<double> e(static_cast<uint32_t>((d1_len * l2y_len) << 1));
    int e_len = product_expansion_zeroelim(d1_len, d1, l2y_len, l2y, e.data());
    std::vector<double> f(static_cast<uint32_t>((d2_len * l1y_len) << 1));
    int f_len = product_expansion_zeroelim(d2_len, d2, l1y_len, l1y, f.data());
    std::vector<double> g(static_cast<uint32_t>((d1_len * l3x_len) << 1));
    int g_len = product_expansion_zeroelim(d1_len, d1, l3x_len, l3x, g.data());
    std::vector<double> h(static_cast<uint32_t>((d3_len * l1x_len) << 1));
    int h_len = product_expansion_zeroelim(d3_len, d3, l1x_len, l1x, h.data());
    std::vector<double> ab(static_cast<uint32_t>(a_len + b_len));
    int ab_len = fast_expansion_diff_zeroelim(a_len, a.data(), b_len, b.data(), ab.data());
    std::vector<double> cd(static_cast<uint32_t>(c_len + d_len));
    int cd_len = fast_expansion_diff_zeroelim(c_len, c.data(), d_len, d.data(), cd.data());
    std::vector<double> ef(static_cast<uint32_t>(e_len + f_len));
    int ef_len = fast_expansion_diff_zeroelim(e_len, e.data(), f_len, f.data(), ef.data());
    std::vector<double> gh(static_cast<uint32_t>(g_len + h_len));
    int gh_len = fast_expansion_diff_zeroelim(g_len, g.data(), h_len, h.data(), gh.data());
    std::vector<double> abcd(static_cast<uint32_t>((ab_len * cd_len) << 1));
    int abcd_len = product_expansion_zeroelim(ab_len, ab.data(), cd_len, cd.data(), abcd.data());
    std::vector<double> efgh(static_cast<uint32_t>((ef_len * gh_len) << 1));
    int efgh_len = product_expansion_zeroelim(ef_len, ef.data(), gh_len, gh.data(), efgh.data());
    std::vector<double> L(static_cast<uint32_t>(abcd_len + efgh_len));
    const uint32_t L_len =
        static_cast<uint32_t>(fast_expansion_diff_zeroelim(abcd_len, abcd.data(), efgh_len, efgh.data(), L.data()));

    if (L[L_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (L[L_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient2d_SSS(const implicitPoint2D_SSI& a, const implicitPoint2D_SSI& b, const implicitPoint2D_SSI& c) {
    int ret;
    if ((ret = orient2d_SSS_filtered(a, b, c)) != IPSign::ZERO) {
        return ret;
    }
    if ((ret = orient2d_SSS_interval(a, b, c)) != IPSign::ZERO) {
        return ret;
    }
    return orien2d_SSS_exact(a, b, c);
}

int GenericPoint2D::orient2D(const GenericPoint2D& a, const GenericPoint2D& b, const GenericPoint2D& c) {
    switch (a.isExplicit() << 2 | b.isExplicit() << 1 | c.isExplicit()) {
    case 0:
        return orient2d_SSS(a.toSSI(), b.toSSI(), c.toSSI());
    case 1:
        return orient2d_SSE(a.toSSI(), b.toSSI(), c.toExplicit2D().ptr());
    case 2:
        return orient2d_SSE(c.toSSI(), a.toSSI(), b.toExplicit2D().ptr());
    case 3:
        return orient2d_SEE(a.toSSI(), b.toExplicit2D().ptr(), c.toExplicit2D().ptr());
    case 4:
        return orient2d_SSE(b.toSSI(), c.toSSI(), a.toExplicit2D().ptr());
    case 5:
        return orient2d_SEE(b.toSSI(), c.toExplicit2D().ptr(), a.toExplicit2D().ptr());
    case 6:
        return orient2d_SEE(c.toSSI(), a.toExplicit2D().ptr(), b.toExplicit2D().ptr());
    case 7:
        return orient2d_EEE(a.toExplicit2D(), b.toExplicit2D(), c.toExplicit2D());
    }
    return IPSign::UNDEFINED;
}
int GenericPoint2D::sign_orient2d(const double* p, const double* q, const double* r) {
    const double ret = orient2d(p, q, r);
    return (ret > .0) - (ret < .0);
}


inline bool
lambda2d_SSI_filtered(const double* ea1, const double* ea2, const double* eb1, const double* eb2, double* filter) {
    double t1a = ea1[0] * ea2[1];
    double t1b = ea2[0] * ea1[1];
    double t1 = t1a - t1b;
    double tx2 = eb1[0] - eb2[0];
    double t3a = eb1[0] * eb2[1];
    double t3b = eb2[0] * eb1[1];
    double t3 = t3a - t3b;
    double tx4 = ea1[0] - ea2[0];
    double ty2 = eb1[1] - eb2[1];
    double ty4 = ea1[1] - ea2[1];
    double lxa = t1 * tx2;
    double lxb = t3 * tx4;
    filter[0] = lxa - lxb;
    double lya = t1 * ty2;
    double lyb = t3 * ty4;
    filter[1] = lya - lyb;
    double deta = tx4 * ty2;
    double detb = tx2 * ty4;
    filter[2] = deta - detb;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(ea1[0])) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ea1[1])) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ea2[0])) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ea2[1])) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(eb1[0])) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(eb1[1])) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(eb2[0])) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(eb2[1])) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(tx2)) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(tx4)) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ty2)) > filter[3]) filter[3] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ty4)) > filter[3]) filter[3] = _tmp_fabs;
    double lambda_det_eps = filter[3];
    lambda_det_eps *= lambda_det_eps;
    lambda_det_eps *= 8.881784197001253e-16;

    return ((filter[2] > lambda_det_eps || filter[2] < -lambda_det_eps));
}

bool implicitPoint2D_SSI::getFilteredLambda(double& mv) const {
    if (needsFilteredLambda()) {
        if (!lambda2d_SSI_filtered(l1_1.ptr(), l1_2.ptr(), l2_1.ptr(), l2_2.ptr(), ssfilter.data())) {
            ssfilter[2] = 0.0;
        }
        if (ssfilter[2] < 0.0) {
            ssfilter[0] = -ssfilter[0];
            ssfilter[1] = -ssfilter[1];
            ssfilter[2] = -ssfilter[2];
        }
    }
    if (mv < ssfilter[3]) {
        mv = ssfilter[3];
    }
    return ssfilter[2] != 0.0;
}

inline bool lambda2d_SSI_interval(
    const IntervalNumber& ea1x, const IntervalNumber& ea1y, const IntervalNumber& ea2x, const IntervalNumber& ea2y,
    const IntervalNumber& eb1x, const IntervalNumber& eb1y, const IntervalNumber& eb2x, const IntervalNumber& eb2y,
    IntervalNumber* filter
) {
    IntervalNumber t1a(ea1x * ea2y);
    IntervalNumber t1b(ea2x * ea1y);
    IntervalNumber t1(t1a - t1b);
    IntervalNumber tx2(eb1x - eb2x);
    IntervalNumber t3a(eb1x * eb2y);
    IntervalNumber t3b(eb2x * eb1y);
    IntervalNumber t3(t3a - t3b);
    IntervalNumber tx4(ea1x - ea2x);
    IntervalNumber ty2(eb1y - eb2y);
    IntervalNumber ty4(ea1y - ea2y);
    IntervalNumber lxa(t1 * tx2);
    IntervalNumber lxb(t3 * tx4);
    filter[0] = lxa - lxb;
    IntervalNumber lya(t1 * ty2);
    IntervalNumber lyb(t3 * ty4);
    filter[1] = lya - lyb;
    IntervalNumber deta(tx4 * ty2);
    IntervalNumber detb(tx2 * ty4);
    filter[2] = deta - detb;

    return filter[2].signIsReliable();
}

bool implicitPoint2D_SSI::getIntervalLambda() const {
    if (needsIntervalLambda()) {

        lambda2d_SSI_interval(
            l1_1.X(), l1_1.Y(), l1_2.X(), l1_2.Y(), l2_1.X(), l2_1.Y(), l2_2.X(), l2_2.Y(), dfilter.data()
        );
        if (dfilter[2].isNegative()) {
            dfilter[0].invert();
            dfilter[1].invert();
            dfilter[2].invert();
        }
    }
    return dfilter[2].signIsReliable();
}

inline void lambda2d_SSI_exact(
    const double* ea, const double* eb, const double* ec, const double* ed, double* lambda_x, int& lambda_x_len,
    double* lambda_y, int& lambda_y_len, double* lambda_det, int& lambda_det_len
) {
    double t1a[2];
    two_prod(ea[0], eb[1], t1a);
    double t1b[2];
    two_prod(eb[0], ea[1], t1b);
    double t1[4];
    two_two_diff(t1a, t1b, t1);
    double tx2[2];
    two_diff(ec[0], ed[0], tx2);
    double t3a[2];
    two_prod(ec[0], ed[1], t3a);
    double t3b[2];
    two_prod(ed[0], ec[1], t3b);
    double t3[4];
    two_two_diff(t3a, t3b, t3);
    double tx4[2];
    two_diff(ea[0], eb[0], tx4);
    double ty2[2];
    two_diff(ec[1], ed[1], ty2);
    double ty4[2];
    two_diff(ea[1], eb[1], ty4);
    double lxa[16];
    int lxa_len = product_expansion_zeroelim(4, t1, 2, tx2, lxa);
    double lxb[16];
    int lxb_len = product_expansion_zeroelim(4, t3, 2, tx4, lxb);
    lambda_x_len = fast_expansion_diff_zeroelim(lxa_len, lxa, lxb_len, lxb, lambda_x);
    double lya[16];
    int lya_len = product_expansion_zeroelim(4, t1, 2, ty2, lya);
    double lyb[16];
    int lyb_len = product_expansion_zeroelim(4, t3, 2, ty4, lyb);
    lambda_y_len = fast_expansion_diff_zeroelim(lya_len, lya, lyb_len, lyb, lambda_y);
    double deta[8];
    int deta_len = product_expansion_zeroelim(2, tx4, 2, ty2, deta);
    double detb[8];
    int detb_len = product_expansion_zeroelim(2, tx2, 2, ty4, detb);
    lambda_det_len = fast_expansion_diff_zeroelim(deta_len, deta, detb_len, detb, lambda_det);
}

void implicitPoint2D_SSI::getExactLambda(double* lx, int& lxl, double* ly, int& lyl, double* d, int& dl) const {
    lambda2d_SSI_exact(l1_1.ptr(), l1_2.ptr(), l2_1.ptr(), l2_2.ptr(), lx, lxl, ly, lyl, d, dl);
    if (d[dl - 1] < 0) {
        invert_expansion(lxl, lx);
        invert_expansion(lyl, ly);
        invert_expansion(dl, d);
    }
}

inline int orient3d_LEEE_filtered(const ImplicitPointLPI& p1, const double* a, const double* b, const double* c) {
    double max_var;
    if (!p1.get_filtered_lambda(max_var)) return 0;

    double dcx = p1.ssfilter[3] * c[0];
    double dcy = p1.ssfilter[3] * c[1];
    double dcz = p1.ssfilter[3] * c[2];
    double ix_cx = p1.ssfilter[0] - dcx;
    double iy_cy = p1.ssfilter[1] - dcy;
    double ax_cx = a[0] - c[0];
    double ay_cy = a[1] - c[1];
    double az_cz = a[2] - c[2];
    double iz_cz = p1.ssfilter[2] - dcz;
    double bx_cx = b[0] - c[0];
    double by_cy = b[1] - c[1];
    double bz_cz = b[2] - c[2];
    double tmc_a = ix_cx * ay_cy;
    double tmc_b = iy_cy * ax_cx;
    double m01 = tmc_a - tmc_b;
    double tmi_a = ix_cx * az_cz;
    double tmi_b = iz_cz * ax_cx;
    double m02 = tmi_a - tmi_b;
    double tma_a = iy_cy * az_cz;
    double tma_b = iz_cz * ay_cy;
    double m12 = tma_a - tma_b;
    double mt1 = m01 * bz_cz;
    double mt2 = m02 * by_cy;
    double mt3 = m12 * bx_cx;
    double mtt = mt2 - mt1;
    double m012 = mtt - mt3;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(c[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(c[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(c[2])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(ax_cx)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(ay_cy)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(az_cz)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(bx_cx)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(by_cy)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(bz_cz)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 1.861039534284405e-13;
    if (m012 > epsilon) return IPSign::POSITIVE;
    if (-m012 > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient3d_LEEE_interval(const ImplicitPointLPI& p1, const double* a, const double* b, const double* c) {
    if (!p1.get_interval_lambda()) return 0;

    IntervalNumber dcx(p1.dfilter[3] * c[0]);
    IntervalNumber dcy(p1.dfilter[3] * c[1]);
    IntervalNumber dcz(p1.dfilter[3] * c[2]);
    IntervalNumber ix_cx(p1.dfilter[0] - dcx);
    IntervalNumber iy_cy(p1.dfilter[1] - dcy);
    IntervalNumber ax_cx(IntervalNumber(a[0]) - c[0]);
    IntervalNumber ay_cy(IntervalNumber(a[1]) - c[1]);
    IntervalNumber az_cz(IntervalNumber(a[2]) - c[2]);
    IntervalNumber iz_cz(p1.dfilter[2] - dcz);
    IntervalNumber bx_cx(IntervalNumber(b[0]) - c[0]);
    IntervalNumber by_cy(IntervalNumber(b[1]) - c[1]);
    IntervalNumber bz_cz(IntervalNumber(b[2]) - c[2]);
    IntervalNumber tmc_a(ix_cx * ay_cy);
    IntervalNumber tmc_b(iy_cy * ax_cx);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(ix_cx * az_cz);
    IntervalNumber tmi_b(iz_cz * ax_cx);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(iy_cy * az_cz);
    IntervalNumber tma_b(iz_cz * ay_cy);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * bz_cz);
    IntervalNumber mt2(m02 * by_cy);
    IntervalNumber mt3(m12 * bx_cx);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);

    return m012.sign();
}

inline int orient3d_LEEE_exact(const ImplicitPointLPI& p1, const double* a, const double* b, const double* c) {
    std::vector<double> l1x, l1y, l1z, d1;
    int l1x_len, l1y_len, l1z_len, d1_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    if (d1.data()[d1_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    std::vector<double> dcx(static_cast<uint32_t>(d1_len << 1));
    int dcx_len = scale_expansion_zeroelim(d1_len, d1.data(), c[0], dcx.data());
    std::vector<double> dcy(static_cast<uint32_t>(d1_len << 1));
    int dcy_len = scale_expansion_zeroelim(d1_len, d1.data(), c[1], dcy.data());
    std::vector<double> dcz(static_cast<uint32_t>(d1_len << 1));
    int dcz_len = scale_expansion_zeroelim(d1_len, d1.data(), c[2], dcz.data());
    std::vector<double> ix_cx(static_cast<uint32_t>(l1x_len + dcx_len));
    int ix_cx_len = fast_expansion_diff_zeroelim(l1x_len, l1x.data(), dcx_len, dcx.data(), ix_cx.data());
    std::vector<double> iy_cy(static_cast<uint32_t>(l1y_len + dcy_len));
    int iy_cy_len = fast_expansion_diff_zeroelim(l1y_len, l1y.data(), dcy_len, dcy.data(), iy_cy.data());
    double ax_cx[2];
    two_diff(a[0], c[0], ax_cx);
    double ay_cy[2];
    two_diff(a[1], c[1], ay_cy);
    double az_cz[2];
    two_diff(a[2], c[2], az_cz);
    std::vector<double> iz_cz(static_cast<uint32_t>(l1z_len + dcz_len));
    int iz_cz_len = fast_expansion_diff_zeroelim(l1z_len, l1z.data(), dcz_len, dcz.data(), iz_cz.data());
    double bx_cx[2];
    two_diff(b[0], c[0], bx_cx);
    double by_cy[2];
    two_diff(b[1], c[1], by_cy);
    double bz_cz[2];
    two_diff(b[2], c[2], bz_cz);
    std::vector<double> tmc_a(static_cast<uint32_t>(ix_cx_len << 2));
    int tmc_a_len = product_expansion_zeroelim(ix_cx_len, ix_cx.data(), 2, ay_cy, tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>(iy_cy_len << 2));
    int tmc_b_len = product_expansion_zeroelim(iy_cy_len, iy_cy.data(), 2, ax_cx, tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>(ix_cx_len << 2));
    int tmi_a_len = product_expansion_zeroelim(ix_cx_len, ix_cx.data(), 2, az_cz, tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>(iz_cz_len << 2));
    int tmi_b_len = product_expansion_zeroelim(iz_cz_len, iz_cz.data(), 2, ax_cx, tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>(iy_cy_len << 2));
    int tma_a_len = product_expansion_zeroelim(iy_cy_len, iy_cy.data(), 2, az_cz, tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>(iz_cz_len << 2));
    int tma_b_len = product_expansion_zeroelim(iz_cz_len, iz_cz.data(), 2, ay_cy, tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>(m01_len << 2));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), 2, bz_cz, mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>(m02_len << 2));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), 2, by_cy, mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>(m12_len << 2));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), 2, bx_cx, mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());

    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient3d_LEEE(const ImplicitPointLPI& p1, const double* p2, const double* p3, const double* p4) {
    int ret;
    if ((ret = orient3d_LEEE_filtered(p1, p2, p3, p4)) != IPSign::ZERO) {
        return ret;
    }
    if ((ret = orient3d_LEEE_interval(p1, p2, p3, p4)) != IPSign::ZERO) {
        return ret;
    }
    return orient3d_LEEE_exact(p1, p2, p3, p4);
}

inline int orient3d_TEEE_filtered(const ImplicitPointTPI& p1, const double* a, const double* b, const double* c) {
    double max_var;
    if (!p1.get_filtered_lambda(max_var)) return 0;

    double dcx = p1.ssfilter[3] * c[0];
    double dcy = p1.ssfilter[3] * c[1];
    double dcz = p1.ssfilter[3] * c[2];
    double ix_cx = p1.ssfilter[0] - dcx;
    double iy_cy = p1.ssfilter[1] - dcy;
    double ax_cx = a[0] - c[0];
    double ay_cy = a[1] - c[1];
    double az_cz = a[2] - c[2];
    double iz_cz = p1.ssfilter[2] - dcz;
    double bx_cx = b[0] - c[0];
    double by_cy = b[1] - c[1];
    double bz_cz = b[2] - c[2];
    double tmc_a = ix_cx * ay_cy;
    double tmc_b = iy_cy * ax_cx;
    double m01 = tmc_a - tmc_b;
    double tmi_a = ix_cx * az_cz;
    double tmi_b = iz_cz * ax_cx;
    double m02 = tmi_a - tmi_b;
    double tma_a = iy_cy * az_cz;
    double tma_b = iz_cz * ay_cy;
    double m12 = tma_a - tma_b;
    double mt1 = m01 * bz_cz;
    double mt2 = m02 * by_cy;
    double mt3 = m12 * bx_cx;
    double mtt = mt2 - mt1;
    double m012 = mtt - mt3;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(c[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(c[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(c[2])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(ax_cx)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(ay_cy)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(az_cz)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(bx_cx)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(by_cy)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= 3.070283610684406e-12;
    if ((_tmp_fabs = fabs(bz_cz)) > max_var) max_var = _tmp_fabs;
    if (m012 > epsilon) return IPSign::POSITIVE;
    if (-m012 > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient3d_TEEE_interval(const ImplicitPointTPI& p1, const double* a, const double* b, const double* c) {
    if (!p1.get_interval_lambda()) return 0;

    IntervalNumber dcx(p1.dfilter[3] * c[0]);
    IntervalNumber dcy(p1.dfilter[3] * c[1]);
    IntervalNumber dcz(p1.dfilter[3] * c[2]);
    IntervalNumber ix_cx(p1.dfilter[0] - dcx);
    IntervalNumber iy_cy(p1.dfilter[1] - dcy);
    IntervalNumber ax_cx(IntervalNumber(a[0]) - c[0]);
    IntervalNumber ay_cy(IntervalNumber(a[1]) - c[1]);
    IntervalNumber az_cz(IntervalNumber(a[2]) - c[2]);
    IntervalNumber iz_cz(p1.dfilter[2] - dcz);
    IntervalNumber bx_cx(IntervalNumber(b[0]) - c[0]);
    IntervalNumber by_cy(IntervalNumber(b[1]) - c[1]);
    IntervalNumber bz_cz(IntervalNumber(b[2]) - c[2]);
    IntervalNumber tmc_a(ix_cx * ay_cy);
    IntervalNumber tmc_b(iy_cy * ax_cx);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(ix_cx * az_cz);
    IntervalNumber tmi_b(iz_cz * ax_cx);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(iy_cy * az_cz);
    IntervalNumber tma_b(iz_cz * ay_cy);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * bz_cz);
    IntervalNumber mt2(m02 * by_cy);
    IntervalNumber mt3(m12 * bx_cx);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);

    return m012.sign();
}

inline int orient3d_TEEE_exact(const ImplicitPointTPI& p1, const double* a, const double* b, const double* c) {
    std::vector<double> l1x, l1y, l1z, d1;
    int l1x_len, l1y_len, l1z_len, d1_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    if (d1.data()[d1_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    std::vector<double> dcx(static_cast<uint32_t>(d1_len << 1));
    int dcx_len = scale_expansion_zeroelim(d1_len, d1.data(), c[0], dcx.data());
    std::vector<double> dcy(static_cast<uint32_t>(d1_len << 1));
    int dcy_len = scale_expansion_zeroelim(d1_len, d1.data(), c[1], dcy.data());
    std::vector<double> dcz(static_cast<uint32_t>(d1_len << 1));
    int dcz_len = scale_expansion_zeroelim(d1_len, d1.data(), c[2], dcz.data());
    std::vector<double> ix_cx(static_cast<uint32_t>(l1x_len + dcx_len));
    int ix_cx_len = fast_expansion_diff_zeroelim(l1x_len, l1x.data(), dcx_len, dcx.data(), ix_cx.data());
    std::vector<double> iy_cy(static_cast<uint32_t>(l1y_len + dcy_len));
    int iy_cy_len = fast_expansion_diff_zeroelim(l1y_len, l1y.data(), dcy_len, dcy.data(), iy_cy.data());
    double ax_cx[2];
    two_diff(a[0], c[0], ax_cx);
    double ay_cy[2];
    two_diff(a[1], c[1], ay_cy);
    double az_cz[2];
    two_diff(a[2], c[2], az_cz);
    std::vector<double> iz_cz(static_cast<uint32_t>(l1z_len + dcz_len));
    int iz_cz_len = fast_expansion_diff_zeroelim(l1z_len, l1z.data(), dcz_len, dcz.data(), iz_cz.data());
    double bx_cx[2];
    two_diff(b[0], c[0], bx_cx);
    double by_cy[2];
    two_diff(b[1], c[1], by_cy);
    double bz_cz[2];
    two_diff(b[2], c[2], bz_cz);
    std::vector<double> tmc_a(static_cast<uint32_t>(ix_cx_len << 2));
    int tmc_a_len = product_expansion_zeroelim(ix_cx_len, ix_cx.data(), 2, ay_cy, tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>(iy_cy_len << 2));
    int tmc_b_len = product_expansion_zeroelim(iy_cy_len, iy_cy.data(), 2, ax_cx, tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>(ix_cx_len << 2));
    int tmi_a_len = product_expansion_zeroelim(ix_cx_len, ix_cx.data(), 2, az_cz, tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>(iz_cz_len << 2));
    int tmi_b_len = product_expansion_zeroelim(iz_cz_len, iz_cz.data(), 2, ax_cx, tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>(iy_cy_len << 2));
    int tma_a_len = product_expansion_zeroelim(iy_cy_len, iy_cy.data(), 2, az_cz, tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>(iz_cz_len << 2));
    int tma_b_len = product_expansion_zeroelim(iz_cz_len, iz_cz.data(), 2, ay_cy, tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>(m01_len << 2));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), 2, bz_cz, mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>(m02_len << 2));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), 2, by_cy, mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>(m12_len << 2));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), 2, bx_cx, mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());

    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient3d_TEEE(const ImplicitPointTPI& p1, const double* p2, const double* p3, const double* p4) {
    int ret;
    if ((ret = orient3d_TEEE_filtered(p1, p2, p3, p4)) != IPSign::ZERO) {
        return ret;
    }
    if ((ret = orient3d_TEEE_interval(p1, p2, p3, p4)) != IPSign::ZERO) {
        return ret;
    }
    return orient3d_TEEE_exact(p1, p2, p3, p4);
}

inline int
orient3d_LLEE_filtered(const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const double* p3, const double* p4) {
    double max_var = 0;
    if (!p1.get_filtered_lambda(max_var) || !p2.get_filtered_lambda(max_var)) return 0;

    double d1p4x = p1.ssfilter[3] * p4[0];
    double d1p4y = p1.ssfilter[3] * p4[1];
    double d1p4z = p1.ssfilter[3] * p4[2];
    double d2p4x = p2.ssfilter[3] * p4[0];
    double d2p4y = p2.ssfilter[3] * p4[1];
    double d2p4z = p2.ssfilter[3] * p4[2];
    double p1p4x = p1.ssfilter[0] - d1p4x;
    double p1p4y = p1.ssfilter[1] - d1p4y;
    double p1p4z = p1.ssfilter[2] - d1p4z;
    double p2p4x = p2.ssfilter[0] - d2p4x;
    double p2p4y = p2.ssfilter[1] - d2p4y;
    double p2p4z = p2.ssfilter[2] - d2p4z;
    double p3p4x = p3[0] - p4[0];
    double p3p4y = p3[1] - p4[1];
    double p3p4z = p3[2] - p4[2];
    double tmc_a = p1p4x * p2p4y;
    double tmc_b = p1p4y * p2p4x;
    double m01 = tmc_a - tmc_b;
    double tmi_a = p1p4x * p2p4z;
    double tmi_b = p1p4z * p2p4x;
    double m02 = tmi_a - tmi_b;
    double tma_a = p1p4y * p2p4z;
    double tma_b = p1p4z * p2p4y;
    double m12 = tma_a - tma_b;
    double mt1 = m01 * p3p4z;
    double mt2 = m02 * p3p4y;
    double mt3 = m12 * p3p4x;
    double mtt = mt2 - mt1;
    double m012 = mtt - mt3;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(p4[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p4[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p4[2])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3p4x)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3p4y)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3p4z)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= 5.12855469897434e-12;
    if (m012 > epsilon) return IPSign::POSITIVE;
    if (-m012 > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int
orient3d_LLEE_interval(const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const double* p3, const double* p4) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda()) return IPSign::ZERO;

    IntervalNumber d1p4x(p1.dfilter[3] * p4[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4[2]);
    IntervalNumber p1p4x(p1.dfilter[0] - d1p4x);
    IntervalNumber p1p4y(p1.dfilter[1] - d1p4y);
    IntervalNumber p1p4z(p1.dfilter[2] - d1p4z);
    IntervalNumber p2p4x(p2.dfilter[0] - d2p4x);
    IntervalNumber p2p4y(p2.dfilter[1] - d2p4y);
    IntervalNumber p2p4z(p2.dfilter[2] - d2p4z);
    IntervalNumber p3p4x(p3[0] - p4[0]);
    IntervalNumber p3p4y(p3[1] - p4[1]);
    IntervalNumber p3p4z(p3[2] - p4[2]);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);

    return m012.sign();
}

inline int
orient3d_LLEE_exact(const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const double* p3, const double* p4) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, d1, d2;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, d1_len, d2_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> d1p4x(static_cast<uint32_t>(d1_len << 1));
    int d1p4x_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[0], d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>(d1_len << 1));
    int d1p4y_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[1], d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>(d1_len << 1));
    int d1p4z_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[2], d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>(d2_len << 1));
    int d2p4x_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[0], d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>(d2_len << 1));
    int d2p4y_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[1], d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>(d2_len << 1));
    int d2p4z_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[2], d2p4z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(l1x_len, l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(l1y_len, l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(l1z_len, l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(l2x_len, l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(l2y_len, l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(l2y_len + d2p4y_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(l2z_len, l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    double p3p4x[2];
    two_diff(p3[0], p4[0], p3p4x);
    double p3p4y[2];
    two_diff(p3[1], p4[1], p3p4y);
    double p3p4z[2];
    two_diff(p3[2], p4[2], p3p4z);
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>(m01_len << 2));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), 2, p3p4z, mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>(m02_len << 2));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), 2, p3p4y, mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>(m12_len << 2));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), 2, p3p4x, mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient3d_LLEE(const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const double* p3, const double* p4) {
    int ret;
    if ((ret = orient3d_LLEE_filtered(p1, p2, p3, p4)) != IPSign::ZERO) {
        return ret;
    }
    if ((ret = orient3d_LLEE_interval(p1, p2, p3, p4)) != IPSign::ZERO) {
        return ret;
    }
    return orient3d_LLEE_exact(p1, p2, p3, p4);
}

inline int
orient3d_LTEE_filtered(const ImplicitPointLPI& p1, const ImplicitPointTPI& p2, const double* p3, const double* p4) {
    double max_var = 0;
    if (!p1.get_filtered_lambda(max_var) || !p2.get_filtered_lambda(max_var)) return 0;

    double d1p4x = p1.ssfilter[3] * p4[0];
    double d1p4y = p1.ssfilter[3] * p4[1];
    double d1p4z = p1.ssfilter[3] * p4[2];
    double d2p4x = p2.ssfilter[3] * p4[0];
    double d2p4y = p2.ssfilter[3] * p4[1];
    double d2p4z = p2.ssfilter[3] * p4[2];
    double p1p4x = p1.ssfilter[0] - d1p4x;
    double p1p4y = p1.ssfilter[1] - d1p4y;
    double p1p4z = p1.ssfilter[2] - d1p4z;
    double p2p4x = p2.ssfilter[0] - d2p4x;
    double p2p4y = p2.ssfilter[1] - d2p4y;
    double p2p4z = p2.ssfilter[2] - d2p4z;
    double p3p4x = p3[0] - p4[0];
    double p3p4y = p3[1] - p4[1];
    double p3p4z = p3[2] - p4[2];
    double tmc_a = p1p4x * p2p4y;
    double tmc_b = p1p4y * p2p4x;
    double m01 = tmc_a - tmc_b;
    double tmi_a = p1p4x * p2p4z;
    double tmi_b = p1p4z * p2p4x;
    double m02 = tmi_a - tmi_b;
    double tma_a = p1p4y * p2p4z;
    double tma_b = p1p4z * p2p4y;
    double m12 = tma_a - tma_b;
    double mt1 = m01 * p3p4z;
    double mt2 = m02 * p3p4y;
    double mt3 = m12 * p3p4x;
    double mtt = mt2 - mt1;
    double m012 = mtt - mt3;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(p4[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p4[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p4[2])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3p4x)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3p4y)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 7.437036403379365e-11;
    if ((_tmp_fabs = fabs(p3p4z)) > max_var) max_var = _tmp_fabs;
    if (m012 > epsilon) return IPSign::POSITIVE;
    if (-m012 > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int
orient3d_LTEE_interval(const ImplicitPointLPI& p1, const ImplicitPointTPI& p2, const double* p3, const double* p4) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda()) return IPSign::ZERO;

    IntervalNumber d1p4x(p1.dfilter[3] * p4[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4[2]);
    IntervalNumber p1p4x(p1.dfilter[0] - d1p4x);
    IntervalNumber p1p4y(p1.dfilter[1] - d1p4y);
    IntervalNumber p1p4z(p1.dfilter[2] - d1p4z);
    IntervalNumber p2p4x(p2.dfilter[0] - d2p4x);
    IntervalNumber p2p4y(p2.dfilter[1] - d2p4y);
    IntervalNumber p2p4z(p2.dfilter[2] - d2p4z);
    IntervalNumber p3p4x(p3[0] - p4[0]);
    IntervalNumber p3p4y(p3[1] - p4[1]);
    IntervalNumber p3p4z(p3[2] - p4[2]);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);

    return m012.sign();
}

inline int
orient3d_LTEE_exact(const ImplicitPointLPI& p1, const ImplicitPointTPI& p2, const double* p3, const double* p4) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, d1, d2;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, d1_len, d2_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> d1p4x(static_cast<uint32_t>(d1_len << 1));
    int d1p4x_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[0], d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>(d1_len << 1));
    int d1p4y_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[1], d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>(d1_len << 1));
    int d1p4z_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[2], d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>(d2_len << 1));
    int d2p4x_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[0], d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>(d2_len << 1));
    int d2p4y_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[1], d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>(d2_len << 1));
    int d2p4z_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[2], d2p4z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(l1x_len, l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(l1y_len, l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(l1z_len, l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(l2x_len, l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(l2y_len, l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(l2y_len + d2p4y_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(l2z_len, l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    double p3p4x[2];
    two_diff(p3[0], p4[0], p3p4x);
    double p3p4y[2];
    two_diff(p3[1], p4[1], p3p4y);
    double p3p4z[2];
    two_diff(p3[2], p4[2], p3p4z);
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>(m01_len << 2));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), 2, p3p4z, mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>(m02_len << 2));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), 2, p3p4y, mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>(m12_len << 2));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), 2, p3p4x, mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient3d_LTEE(const ImplicitPointLPI& p1, const ImplicitPointTPI& p2, const double* p3, const double* p4) {
    int ret;
    if ((ret = orient3d_LTEE_filtered(p1, p2, p3, p4)) != IPSign::ZERO) {
        return ret;
    }
    if ((ret = orient3d_LTEE_interval(p1, p2, p3, p4)) != IPSign::ZERO) {
        return ret;
    }
    return orient3d_LTEE_exact(p1, p2, p3, p4);
}

inline int
orient3d_TTEE_filtered(const ImplicitPointTPI& p1, const ImplicitPointTPI& p2, const double* p3, const double* p4) {
    double max_var = 0;
    if (!p1.get_filtered_lambda(max_var) || !p2.get_filtered_lambda(max_var)) return 0;

    double d1p4x = p1.ssfilter[3] * p4[0];
    double d1p4y = p1.ssfilter[3] * p4[1];
    double d1p4z = p1.ssfilter[3] * p4[2];
    double d2p4x = p2.ssfilter[3] * p4[0];
    double d2p4y = p2.ssfilter[3] * p4[1];
    double d2p4z = p2.ssfilter[3] * p4[2];
    double p1p4x = p1.ssfilter[0] - d1p4x;
    double p1p4y = p1.ssfilter[1] - d1p4y;
    double p1p4z = p1.ssfilter[2] - d1p4z;
    double p2p4x = p2.ssfilter[0] - d2p4x;
    double p2p4y = p2.ssfilter[1] - d2p4y;
    double p2p4z = p2.ssfilter[2] - d2p4z;
    double p3p4x = p3[0] - p4[0];
    double p3p4y = p3[1] - p4[1];
    double p3p4z = p3[2] - p4[2];
    double tmc_a = p1p4x * p2p4y;
    double tmc_b = p1p4y * p2p4x;
    double m01 = tmc_a - tmc_b;
    double tmi_a = p1p4x * p2p4z;
    double tmi_b = p1p4z * p2p4x;
    double m02 = tmi_a - tmi_b;
    double tma_a = p1p4y * p2p4z;
    double tma_b = p1p4z * p2p4y;
    double m12 = tma_a - tma_b;
    double mt1 = m01 * p3p4z;
    double mt2 = m02 * p3p4y;
    double mt3 = m12 * p3p4x;
    double mtt = mt2 - mt1;
    double m012 = mtt - mt3;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(p4[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p4[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p4[2])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3p4x)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3p4y)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 1.036198238324465e-09;
    if ((_tmp_fabs = fabs(p3p4z)) > max_var) max_var = _tmp_fabs;
    if (m012 > epsilon) return IPSign::POSITIVE;
    if (-m012 > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int
orient3d_TTEE_interval(const ImplicitPointTPI& p1, const ImplicitPointTPI& p2, const double* p3, const double* p4) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda()) return IPSign::ZERO;

    IntervalNumber d1p4x(p1.dfilter[3] * p4[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4[2]);
    IntervalNumber p1p4x(p1.dfilter[0] - d1p4x);
    IntervalNumber p1p4y(p1.dfilter[1] - d1p4y);
    IntervalNumber p1p4z(p1.dfilter[2] - d1p4z);
    IntervalNumber p2p4x(p2.dfilter[0] - d2p4x);
    IntervalNumber p2p4y(p2.dfilter[1] - d2p4y);
    IntervalNumber p2p4z(p2.dfilter[2] - d2p4z);
    IntervalNumber p3p4x(p3[0] - p4[0]);
    IntervalNumber p3p4y(p3[1] - p4[1]);
    IntervalNumber p3p4z(p3[2] - p4[2]);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);

    return m012.sign();
}

inline int
orient3d_TTEE_exact(const ImplicitPointTPI& p1, const ImplicitPointTPI& p2, const double* p3, const double* p4) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, d1, d2;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, d1_len, d2_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> d1p4x(static_cast<uint32_t>(d1_len << 1));
    int d1p4x_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[0], d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>(d1_len << 1));
    int d1p4y_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[1], d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>(d1_len << 1));
    int d1p4z_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[2], d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>(d2_len << 1));
    int d2p4x_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[0], d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>(d2_len << 1));
    int d2p4y_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[1], d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>(d2_len << 1));
    int d2p4z_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[2], d2p4z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(l1x_len, l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(l1y_len, l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(l1z_len, l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(l2x_len, l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(l2y_len, l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(l2y_len + d2p4y_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(l2z_len, l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    double p3p4x[2];
    two_diff(p3[0], p4[0], p3p4x);
    double p3p4y[2];
    two_diff(p3[1], p4[1], p3p4y);
    double p3p4z[2];
    two_diff(p3[2], p4[2], p3p4z);
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>(m01_len << 2));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), 2, p3p4z, mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>(m02_len << 2));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), 2, p3p4y, mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>(m12_len << 2));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), 2, p3p4x, mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient3d_TTEE(const ImplicitPointTPI& p1, const ImplicitPointTPI& p2, const double* p3, const double* p4) {
    int ret;
    if ((ret = orient3d_TTEE_filtered(p1, p2, p3, p4)) != IPSign::ZERO) {
        return ret;
    }
    if ((ret = orient3d_TTEE_interval(p1, p2, p3, p4)) != IPSign::ZERO) {
        return ret;
    }
    return orient3d_TTEE_exact(p1, p2, p3, p4);
}

inline int orient3d_LLLE_interval(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointLPI& p3, const double*& p4
) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda() || !p3.get_interval_lambda()) {
        return IPSign::ZERO;
    }

    IntervalNumber d1p4x(p1.dfilter[3] * p4[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4[2]);
    IntervalNumber d3p4x(p3.dfilter[3] * p4[0]);
    IntervalNumber d3p4y(p3.dfilter[3] * p4[1]);
    IntervalNumber d3p4z(p3.dfilter[3] * p4[2]);
    IntervalNumber p1p4x(p1.dfilter[0] - d1p4x);
    IntervalNumber p1p4y(p1.dfilter[1] - d1p4y);
    IntervalNumber p1p4z(p1.dfilter[2] - d1p4z);
    IntervalNumber p2p4x(p2.dfilter[0] - d2p4x);
    IntervalNumber p2p4y(p2.dfilter[1] - d2p4y);
    IntervalNumber p2p4z(p2.dfilter[2] - d2p4z);
    IntervalNumber p3p4x(p3.dfilter[0] - d3p4x);
    IntervalNumber p3p4y(p3.dfilter[1] - d3p4y);
    IntervalNumber p3p4z(p3.dfilter[2] - d3p4z);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);
    return m012.sign();
}

inline int orient3d_LLLE_exact(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointLPI& p3, const double* p4
) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, l3x, l3y, l3z, d1, d2, d3;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, l3x_len, l3y_len, l3z_len, d1_len, d2_len, d3_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    p3.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0 || d3.data()[d3_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> d1p4x(static_cast<uint32_t>(d1_len << 1));
    int d1p4x_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[0], d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>(d1_len << 1));
    int d1p4y_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[1], d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>(d1_len << 1));
    int d1p4z_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[2], d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>(d2_len << 1));
    int d2p4x_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[0], d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>(d2_len << 1));
    int d2p4y_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[1], d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>(d2_len << 1));
    int d2p4z_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[2], d2p4z.data());
    std::vector<double> d3p4x(static_cast<uint32_t>(d3_len << 1));
    int d3p4x_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[0], d3p4x.data());
    std::vector<double> d3p4y(static_cast<uint32_t>(d3_len << 1));
    int d3p4y_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[1], d3p4y.data());
    std::vector<double> d3p4z(static_cast<uint32_t>(d3_len << 1));
    int d3p4z_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[2], d3p4z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(l1x_len, l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(l1y_len, l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(l1z_len, l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(l2x_len, l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(l2y_len, l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(l2z_len + d2p4z_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(l2z_len, l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    std::vector<double> p3p4x(static_cast<uint32_t>(l3x_len + d3p4x_len));
    int p3p4x_len = fast_expansion_diff_zeroelim(l3x_len, l3x.data(), d3p4x_len, d3p4x.data(), p3p4x.data());
    std::vector<double> p3p4y(static_cast<uint32_t>(l3y_len + d3p4y_len));
    int p3p4y_len = fast_expansion_diff_zeroelim(l3y_len, l3y.data(), d3p4y_len, d3p4y.data(), p3p4y.data());
    std::vector<double> p3p4z(static_cast<uint32_t>(l3z_len + d3p4z_len));
    int p3p4z_len = fast_expansion_diff_zeroelim(l3z_len, l3z.data(), d3p4z_len, d3p4z.data(), p3p4z.data());
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>((m01_len * p3p4z_len) << 1));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), p3p4z_len, p3p4z.data(), mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>((m02_len * p3p4y_len) << 1));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), p3p4y_len, p3p4y.data(), mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>((m12_len * p3p4x_len) << 1));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), p3p4x_len, p3p4x.data(), mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int
orient3d_LLLE(const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointLPI& p3, const double* p4) {
    int ret = orient3d_LLLE_interval(p1, p2, p3, p4);
    if (ret != 0) {
        return ret;
    }
    return orient3d_LLLE_exact(p1, p2, p3, p4);
}

inline int orient3d_LLTE_interval(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointTPI& p3, const double* p4
) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda() || !p3.get_interval_lambda()) {
        return IPSign::ZERO;
    }
    IntervalNumber d1p4x(p1.dfilter[3] * p4[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4[2]);
    IntervalNumber d3p4x(p3.dfilter[3] * p4[0]);
    IntervalNumber d3p4y(p3.dfilter[3] * p4[1]);
    IntervalNumber d3p4z(p3.dfilter[3] * p4[2]);
    IntervalNumber p1p4x(p1.dfilter[0] - d1p4x);
    IntervalNumber p1p4y(p1.dfilter[1] - d1p4y);
    IntervalNumber p1p4z(p1.dfilter[2] - d1p4z);
    IntervalNumber p2p4x(p2.dfilter[0] - d2p4x);
    IntervalNumber p2p4y(p2.dfilter[1] - d2p4y);
    IntervalNumber p2p4z(p2.dfilter[2] - d2p4z);
    IntervalNumber p3p4x(p3.dfilter[0] - d3p4x);
    IntervalNumber p3p4y(p3.dfilter[1] - d3p4y);
    IntervalNumber p3p4z(p3.dfilter[2] - d3p4z);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);
    return m012.sign();
}

inline int orient3d_LLTE_exact(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointTPI& p3, const double* p4
) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, l3x, l3y, l3z, d1, d2, d3;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, l3x_len, l3y_len, l3z_len, d1_len, d2_len, d3_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    p3.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0 || d3.data()[d3_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> d1p4x(static_cast<uint32_t>(d1_len << 1));
    int d1p4x_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[0], d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>(d1_len << 1));
    int d1p4y_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[1], d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>(d1_len << 1));
    int d1p4z_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[2], d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>(d2_len << 1));
    int d2p4x_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[0], d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>(d2_len << 1));
    int d2p4y_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[1], d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>(d2_len << 1));
    int d2p4z_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[2], d2p4z.data());
    std::vector<double> d3p4x(static_cast<uint32_t>(d3_len << 1));
    int d3p4x_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[0], d3p4x.data());
    std::vector<double> d3p4y(static_cast<uint32_t>(d3_len << 1));
    int d3p4y_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[1], d3p4y.data());
    std::vector<double> d3p4z(static_cast<uint32_t>(d3_len << 1));
    int d3p4z_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[2], d3p4z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(l1x_len, l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(l1y_len, l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(l1z_len, l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(l2x_len, l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(l2y_len, l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(l2z_len + d2p4z_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(l2z_len, l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    std::vector<double> p3p4x(static_cast<uint32_t>(l3x_len + d3p4x_len));
    int p3p4x_len = fast_expansion_diff_zeroelim(l3x_len, l3x.data(), d3p4x_len, d3p4x.data(), p3p4x.data());
    std::vector<double> p3p4y(static_cast<uint32_t>(l3y_len + d3p4y_len));
    int p3p4y_len = fast_expansion_diff_zeroelim(l3y_len, l3y.data(), d3p4y_len, d3p4y.data(), p3p4y.data());
    std::vector<double> p3p4z(static_cast<uint32_t>(l3z_len + d3p4z_len));
    int p3p4z_len = fast_expansion_diff_zeroelim(l3z_len, l3z.data(), d3p4z_len, d3p4z.data(), p3p4z.data());
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>((m01_len * p3p4z_len) << 1));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), p3p4z_len, p3p4z.data(), mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>((m02_len * p3p4y_len) << 1));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), p3p4y_len, p3p4y.data(), mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>((m12_len * p3p4x_len) << 1));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), p3p4x_len, p3p4x.data(), mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int
orient3d_LLTE(const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointTPI& p3, const double* p4) {

    int ret = orient3d_LLTE_interval(p1, p2, p3, p4);
    if (ret != 0) {
        return ret;
    }
    return orient3d_LLTE_exact(p1, p2, p3, p4);
}

inline int orient3d_LTTE_interval(
    const ImplicitPointLPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const double* p4
) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda() || !p3.get_interval_lambda()) {
        return IPSign::ZERO;
    }
    IntervalNumber d1p4x(p1.dfilter[3] * p4[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4[2]);
    IntervalNumber d3p4x(p3.dfilter[3] * p4[0]);
    IntervalNumber d3p4y(p3.dfilter[3] * p4[1]);
    IntervalNumber d3p4z(p3.dfilter[3] * p4[2]);
    IntervalNumber p1p4x(p1.dfilter[0] - d1p4x);
    IntervalNumber p1p4y(p1.dfilter[1] - d1p4y);
    IntervalNumber p1p4z(p1.dfilter[2] - d1p4z);
    IntervalNumber p2p4x(p2.dfilter[0] - d2p4x);
    IntervalNumber p2p4y(p2.dfilter[1] - d2p4y);
    IntervalNumber p2p4z(p2.dfilter[2] - d2p4z);
    IntervalNumber p3p4x(p3.dfilter[0] - d3p4x);
    IntervalNumber p3p4y(p3.dfilter[1] - d3p4y);
    IntervalNumber p3p4z(p3.dfilter[2] - d3p4z);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);
    return m012.sign();
}

inline int orient3d_LTTE_exact(
    const ImplicitPointLPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const double* p4
) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, l3x, l3y, l3z, d1, d2, d3;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, l3x_len, l3y_len, l3z_len, d1_len, d2_len, d3_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    p3.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0 || d3.data()[d3_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> d1p4x(static_cast<uint32_t>(d1_len << 1));
    int d1p4x_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[0], d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>(d1_len << 1));
    int d1p4y_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[1], d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>(d1_len << 1));
    int d1p4z_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[2], d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>(d2_len << 1));
    int d2p4x_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[0], d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>(d2_len << 1));
    int d2p4y_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[1], d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>(d2_len << 1));
    int d2p4z_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[2], d2p4z.data());
    std::vector<double> d3p4x(static_cast<uint32_t>(d3_len << 1));
    int d3p4x_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[0], d3p4x.data());
    std::vector<double> d3p4y(static_cast<uint32_t>(d3_len << 1));
    int d3p4y_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[1], d3p4y.data());
    std::vector<double> d3p4z(static_cast<uint32_t>(d3_len << 1));
    int d3p4z_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[2], d3p4z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(l1x_len, l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(l1y_len, l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(l1z_len, l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(l2x_len, l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(l2y_len, l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(l2z_len + d2p4z_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(l2z_len, l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    std::vector<double> p3p4x(static_cast<uint32_t>(l3x_len + d3p4x_len));
    int p3p4x_len = fast_expansion_diff_zeroelim(l3x_len, l3x.data(), d3p4x_len, d3p4x.data(), p3p4x.data());
    std::vector<double> p3p4y(static_cast<uint32_t>(l3y_len + d3p4y_len));
    int p3p4y_len = fast_expansion_diff_zeroelim(l3y_len, l3y.data(), d3p4y_len, d3p4y.data(), p3p4y.data());
    std::vector<double> p3p4z(static_cast<uint32_t>(l3z_len + d3p4z_len));
    int p3p4z_len = fast_expansion_diff_zeroelim(l3z_len, l3z.data(), d3p4z_len, d3p4z.data(), p3p4z.data());
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>((m01_len * p3p4z_len) << 1));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), p3p4z_len, p3p4z.data(), mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>((m02_len * p3p4y_len) << 1));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), p3p4y_len, p3p4y.data(), mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>((m12_len * p3p4x_len) << 1));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), p3p4x_len, p3p4x.data(), mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int
orient3d_LTTE(const ImplicitPointLPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const double* p4) {
    int ret = orient3d_LTTE_interval(p1, p2, p3, p4);
    if (ret != 0) {
        return ret;
    }
    return orient3d_LTTE_exact(p1, p2, p3, p4);
}

inline int orient3d_TTTE_interval(
    const ImplicitPointTPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const double* p4
) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda() || !p3.get_interval_lambda()) {
        return IPSign::ZERO;
    }
    IntervalNumber d1p4x(p1.dfilter[3] * p4[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4[2]);
    IntervalNumber d3p4x(p3.dfilter[3] * p4[0]);
    IntervalNumber d3p4y(p3.dfilter[3] * p4[1]);
    IntervalNumber d3p4z(p3.dfilter[3] * p4[2]);
    IntervalNumber p1p4x(p1.dfilter[0] - d1p4x);
    IntervalNumber p1p4y(p1.dfilter[1] - d1p4y);
    IntervalNumber p1p4z(p1.dfilter[2] - d1p4z);
    IntervalNumber p2p4x(p2.dfilter[0] - d2p4x);
    IntervalNumber p2p4y(p2.dfilter[1] - d2p4y);
    IntervalNumber p2p4z(p2.dfilter[2] - d2p4z);
    IntervalNumber p3p4x(p3.dfilter[0] - d3p4x);
    IntervalNumber p3p4y(p3.dfilter[1] - d3p4y);
    IntervalNumber p3p4z(p3.dfilter[2] - d3p4z);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);
    return m012.sign();
}

inline int orient3d_TTTE_exact(
    const ImplicitPointTPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const double* p4
) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, l3x, l3y, l3z, d1, d2, d3;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, l3x_len, l3y_len, l3z_len, d1_len, d2_len, d3_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    p3.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0 || d3.data()[d3_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> d1p4x(static_cast<uint32_t>(d1_len << 1));
    int d1p4x_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[0], d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>(d1_len << 1));
    int d1p4y_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[1], d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>(d1_len << 1));
    int d1p4z_len = scale_expansion_zeroelim(d1_len, d1.data(), p4[2], d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>(d2_len << 1));
    int d2p4x_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[0], d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>(d2_len << 1));
    int d2p4y_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[1], d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>(d2_len << 1));
    int d2p4z_len = scale_expansion_zeroelim(d2_len, d2.data(), p4[2], d2p4z.data());
    std::vector<double> d3p4x(static_cast<uint32_t>(d3_len << 1));
    int d3p4x_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[0], d3p4x.data());
    std::vector<double> d3p4y(static_cast<uint32_t>(d3_len << 1));
    int d3p4y_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[1], d3p4y.data());
    std::vector<double> d3p4z(static_cast<uint32_t>(d3_len << 1));
    int d3p4z_len = scale_expansion_zeroelim(d3_len, d3.data(), p4[2], d3p4z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(l1x_len, l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(l1y_len, l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(l1z_len, l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(l2x_len, l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(l2y_len, l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(l2z_len + d2p4z_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(l2z_len, l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    std::vector<double> p3p4x(static_cast<uint32_t>(l3x_len + d3p4x_len));
    int p3p4x_len = fast_expansion_diff_zeroelim(l3x_len, l3x.data(), d3p4x_len, d3p4x.data(), p3p4x.data());
    std::vector<double> p3p4y(static_cast<uint32_t>(l3y_len + d3p4y_len));
    int p3p4y_len = fast_expansion_diff_zeroelim(l3y_len, l3y.data(), d3p4y_len, d3p4y.data(), p3p4y.data());
    std::vector<double> p3p4z(static_cast<uint32_t>(l3z_len + d3p4z_len));
    int p3p4z_len = fast_expansion_diff_zeroelim(l3z_len, l3z.data(), d3p4z_len, d3p4z.data(), p3p4z.data());
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>((m01_len * p3p4z_len) << 1));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), p3p4z_len, p3p4z.data(), mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>((m02_len * p3p4y_len) << 1));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), p3p4y_len, p3p4y.data(), mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>((m12_len * p3p4x_len) << 1));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), p3p4x_len, p3p4x.data(), mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int
orient3d_TTTE(const ImplicitPointTPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const double* p4) {
    int ret = orient3d_TTTE_interval(p1, p2, p3, p4);
    if (ret != 0) {
        return ret;
    }
    return orient3d_TTTE_exact(p1, p2, p3, p4);
}

inline int orient3d_LLLL_interval(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointLPI& p3, const ImplicitPointLPI& p4
) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda() || !p3.get_interval_lambda() ||
        !p4.get_interval_lambda()) {
        return IPSign::ZERO;
    }

    IntervalNumber d1p4x(p1.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d3p4x(p3.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d3p4y(p3.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d3p4z(p3.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d4l1x(p4.dfilter[3] * p1.dfilter[0]);
    IntervalNumber d4l1y(p4.dfilter[3] * p1.dfilter[1]);
    IntervalNumber d4l1z(p4.dfilter[3] * p1.dfilter[2]);
    IntervalNumber d4l2x(p4.dfilter[3] * p2.dfilter[0]);
    IntervalNumber d4l2y(p4.dfilter[3] * p2.dfilter[1]);
    IntervalNumber d4l2z(p4.dfilter[3] * p2.dfilter[2]);
    IntervalNumber d4l3x(p4.dfilter[3] * p3.dfilter[0]);
    IntervalNumber d4l3y(p4.dfilter[3] * p3.dfilter[1]);
    IntervalNumber d4l3z(p4.dfilter[3] * p3.dfilter[2]);
    IntervalNumber p1p4x(d4l1x - d1p4x);
    IntervalNumber p1p4y(d4l1y - d1p4y);
    IntervalNumber p1p4z(d4l1z - d1p4z);
    IntervalNumber p2p4x(d4l2x - d2p4x);
    IntervalNumber p2p4y(d4l2y - d2p4y);
    IntervalNumber p2p4z(d4l2z - d2p4z);
    IntervalNumber p3p4x(d4l3x - d3p4x);
    IntervalNumber p3p4y(d4l3y - d3p4y);
    IntervalNumber p3p4z(d4l3z - d3p4z);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);
    if (!m012.signIsReliable()) {
        return IPSign::ZERO;
    }
    return m012.sign();
}

inline int orient3d_LLLL_exact(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointLPI& p3, const ImplicitPointLPI& p4
) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, l3x, l3y, l3z, l4x, l4y, l4z, d1, d2, d3, d4;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, l3x_len, l3y_len, l3z_len, l4x_len, l4y_len, l4z_len,
        d1_len, d2_len, d3_len, d4_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    p3.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    p4.get_exact_lambda(l4x, l4x_len, l4y, l4y_len, l4z, l4z_len, d4, d4_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0 || d3.data()[d3_len - 1] == 0.0 ||
        d4.data()[d4_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> d1p4x(static_cast<uint32_t>((d1_len * l4x_len) << 1));
    int d1p4x_len = product_expansion_zeroelim(d1_len, d1.data(), l4x_len, l4x.data(), d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>((d1_len * l4y_len) << 1));
    int d1p4y_len = product_expansion_zeroelim(d1_len, d1.data(), l4y_len, l4y.data(), d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>((d1_len * l4z_len) << 1));
    int d1p4z_len = product_expansion_zeroelim(d1_len, d1.data(), l4z_len, l4z.data(), d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>((d2_len * l4x_len) << 1));
    int d2p4x_len = product_expansion_zeroelim(d2_len, d2.data(), l4x_len, l4x.data(), d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>((d2_len * l4y_len) << 1));
    int d2p4y_len = product_expansion_zeroelim(d2_len, d2.data(), l4y_len, l4y.data(), d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>((d2_len * l4z_len) << 1));
    int d2p4z_len = product_expansion_zeroelim(d2_len, d2.data(), l4z_len, l4z.data(), d2p4z.data());
    std::vector<double> d3p4x(static_cast<uint32_t>((d3_len * l4x_len) << 1));
    int d3p4x_len = product_expansion_zeroelim(d3_len, d3.data(), l4x_len, l4x.data(), d3p4x.data());
    std::vector<double> d3p4y(static_cast<uint32_t>((d3_len * l4y_len) << 1));
    int d3p4y_len = product_expansion_zeroelim(d3_len, d3.data(), l4y_len, l4y.data(), d3p4y.data());
    std::vector<double> d3p4z(static_cast<uint32_t>((d3_len * l4z_len) << 1));
    int d3p4z_len = product_expansion_zeroelim(d3_len, d3.data(), l4z_len, l4z.data(), d3p4z.data());
    std::vector<double> d4l1x(static_cast<uint32_t>((d4_len * l1x_len) << 1));
    int d4l1x_len = product_expansion_zeroelim(d4_len, d4.data(), l1x_len, l1x.data(), d4l1x.data());
    std::vector<double> d4l1y(static_cast<uint32_t>((d4_len * l1y_len) << 1));
    int d4l1y_len = product_expansion_zeroelim(d4_len, d4.data(), l1y_len, l1y.data(), d4l1y.data());
    std::vector<double> d4l1z(static_cast<uint32_t>((d4_len * l1z_len) << 1));
    int d4l1z_len = product_expansion_zeroelim(d4_len, d4.data(), l1z_len, l1z.data(), d4l1z.data());
    std::vector<double> d4l2x(static_cast<uint32_t>((d4_len * l2x_len) << 1));
    int d4l2x_len = product_expansion_zeroelim(d4_len, d4.data(), l2x_len, l2x.data(), d4l2x.data());
    std::vector<double> d4l2y(static_cast<uint32_t>((d4_len * l2y_len) << 1));
    int d4l2y_len = product_expansion_zeroelim(d4_len, d4.data(), l2y_len, l2y.data(), d4l2y.data());
    std::vector<double> d4l2z(static_cast<uint32_t>((d4_len * l2z_len) << 1));
    int d4l2z_len = product_expansion_zeroelim(d4_len, d4.data(), l2z_len, l2z.data(), d4l2z.data());
    std::vector<double> d4l3x(static_cast<uint32_t>((d4_len * l3x_len) << 1));
    int d4l3x_len = product_expansion_zeroelim(d4_len, d4.data(), l3x_len, l3x.data(), d4l3x.data());
    std::vector<double> d4l3y(static_cast<uint32_t>((d4_len * l3y_len) << 1));
    int d4l3y_len = product_expansion_zeroelim(d4_len, d4.data(), l3y_len, l3y.data(), d4l3y.data());
    std::vector<double> d4l3z(static_cast<uint32_t>((d4_len * l3z_len) << 1));
    int d4l3z_len = product_expansion_zeroelim(d4_len, d4.data(), l3z_len, l3z.data(), d4l3z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(d4l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(d4l1x_len, d4l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(d4l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(d4l1y_len, d4l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(d4l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(d4l1z_len, d4l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(d4l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(d4l2x_len, d4l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(d4l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(d4l2y_len, d4l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(d4l2z_len + d2p4z_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(d4l2z_len, d4l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    std::vector<double> p3p4x(static_cast<uint32_t>(d4l3x_len + d3p4x_len));
    int p3p4x_len = fast_expansion_diff_zeroelim(d4l3x_len, d4l3x.data(), d3p4x_len, d3p4x.data(), p3p4x.data());
    std::vector<double> p3p4y(static_cast<uint32_t>(d4l3y_len + d3p4y_len));
    int p3p4y_len = fast_expansion_diff_zeroelim(d4l3y_len, d4l3y.data(), d3p4y_len, d3p4y.data(), p3p4y.data());
    std::vector<double> p3p4z(static_cast<uint32_t>(d4l3z_len + d3p4z_len));
    int p3p4z_len = fast_expansion_diff_zeroelim(d4l3z_len, d4l3z.data(), d3p4z_len, d3p4z.data(), p3p4z.data());
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>((m01_len * p3p4z_len) << 1));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), p3p4z_len, p3p4z.data(), mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>((m02_len * p3p4y_len) << 1));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), p3p4y_len, p3p4y.data(), mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>((m12_len * p3p4x_len) << 1));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), p3p4x_len, p3p4x.data(), mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient3d_LLLL(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointLPI& p3, const ImplicitPointLPI& p4
) {
    int ret = orient3d_LLLL_interval(p1, p2, p3, p4);
    if (ret != 0) {
        return ret;
    }
    return orient3d_LLLL_exact(p1, p2, p3, p4);
}

inline int orient3d_LLLT_interval(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointLPI& p3, const ImplicitPointTPI& p4
) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda() || !p3.get_interval_lambda() ||
        !p4.get_interval_lambda()) {
        return IPSign::ZERO;
    }

    IntervalNumber d1p4x(p1.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d3p4x(p3.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d3p4y(p3.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d3p4z(p3.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d4l1x(p4.dfilter[3] * p1.dfilter[0]);
    IntervalNumber d4l1y(p4.dfilter[3] * p1.dfilter[1]);
    IntervalNumber d4l1z(p4.dfilter[3] * p1.dfilter[2]);
    IntervalNumber d4l2x(p4.dfilter[3] * p2.dfilter[0]);
    IntervalNumber d4l2y(p4.dfilter[3] * p2.dfilter[1]);
    IntervalNumber d4l2z(p4.dfilter[3] * p2.dfilter[2]);
    IntervalNumber d4l3x(p4.dfilter[3] * p3.dfilter[0]);
    IntervalNumber d4l3y(p4.dfilter[3] * p3.dfilter[1]);
    IntervalNumber d4l3z(p4.dfilter[3] * p3.dfilter[2]);
    IntervalNumber p1p4x(d4l1x - d1p4x);
    IntervalNumber p1p4y(d4l1y - d1p4y);
    IntervalNumber p1p4z(d4l1z - d1p4z);
    IntervalNumber p2p4x(d4l2x - d2p4x);
    IntervalNumber p2p4y(d4l2y - d2p4y);
    IntervalNumber p2p4z(d4l2z - d2p4z);
    IntervalNumber p3p4x(d4l3x - d3p4x);
    IntervalNumber p3p4y(d4l3y - d3p4y);
    IntervalNumber p3p4z(d4l3z - d3p4z);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);
    return m012.sign();
}

inline int orient3d_LLLT_exact(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointLPI& p3, const ImplicitPointTPI& p4
) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, l3x, l3y, l3z, l4x, l4y, l4z, d1, d2, d3, d4;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, l3x_len, l3y_len, l3z_len, l4x_len, l4y_len, l4z_len,
        d1_len, d2_len, d3_len, d4_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    p3.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    p4.get_exact_lambda(l4x, l4x_len, l4y, l4y_len, l4z, l4z_len, d4, d4_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0 || d3.data()[d3_len - 1] == 0.0 ||
        d4.data()[d4_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    std::vector<double> d1p4x(static_cast<uint32_t>((d1_len * l4x_len) << 1));
    int d1p4x_len = product_expansion_zeroelim(d1_len, d1.data(), l4x_len, l4x.data(), d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>((d1_len * l4y_len) << 1));
    int d1p4y_len = product_expansion_zeroelim(d1_len, d1.data(), l4y_len, l4y.data(), d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>((d1_len * l4z_len) << 1));
    int d1p4z_len = product_expansion_zeroelim(d1_len, d1.data(), l4z_len, l4z.data(), d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>((d2_len * l4x_len) << 1));
    int d2p4x_len = product_expansion_zeroelim(d2_len, d2.data(), l4x_len, l4x.data(), d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>((d2_len * l4y_len) << 1));
    int d2p4y_len = product_expansion_zeroelim(d2_len, d2.data(), l4y_len, l4y.data(), d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>((d2_len * l4z_len) << 1));
    int d2p4z_len = product_expansion_zeroelim(d2_len, d2.data(), l4z_len, l4z.data(), d2p4z.data());
    std::vector<double> d3p4x(static_cast<uint32_t>((d3_len * l4x_len) << 1));
    int d3p4x_len = product_expansion_zeroelim(d3_len, d3.data(), l4x_len, l4x.data(), d3p4x.data());
    std::vector<double> d3p4y(static_cast<uint32_t>((d3_len * l4y_len) << 1));
    int d3p4y_len = product_expansion_zeroelim(d3_len, d3.data(), l4y_len, l4y.data(), d3p4y.data());
    std::vector<double> d3p4z(static_cast<uint32_t>((d3_len * l4z_len) << 1));
    int d3p4z_len = product_expansion_zeroelim(d3_len, d3.data(), l4z_len, l4z.data(), d3p4z.data());
    std::vector<double> d4l1x(static_cast<uint32_t>((d4_len * l1x_len) << 1));
    int d4l1x_len = product_expansion_zeroelim(d4_len, d4.data(), l1x_len, l1x.data(), d4l1x.data());
    std::vector<double> d4l1y(static_cast<uint32_t>((d4_len * l1y_len) << 1));
    int d4l1y_len = product_expansion_zeroelim(d4_len, d4.data(), l1y_len, l1y.data(), d4l1y.data());
    std::vector<double> d4l1z(static_cast<uint32_t>((d4_len * l1z_len) << 1));
    int d4l1z_len = product_expansion_zeroelim(d4_len, d4.data(), l1z_len, l1z.data(), d4l1z.data());
    std::vector<double> d4l2x(static_cast<uint32_t>((d4_len * l2x_len) << 1));
    int d4l2x_len = product_expansion_zeroelim(d4_len, d4.data(), l2x_len, l2x.data(), d4l2x.data());
    std::vector<double> d4l2y(static_cast<uint32_t>((d4_len * l2y_len) << 1));
    int d4l2y_len = product_expansion_zeroelim(d4_len, d4.data(), l2y_len, l2y.data(), d4l2y.data());
    std::vector<double> d4l2z(static_cast<uint32_t>((d4_len * l2z_len) << 1));
    int d4l2z_len = product_expansion_zeroelim(d4_len, d4.data(), l2z_len, l2z.data(), d4l2z.data());
    std::vector<double> d4l3x(static_cast<uint32_t>((d4_len * l3x_len) << 1));
    int d4l3x_len = product_expansion_zeroelim(d4_len, d4.data(), l3x_len, l3x.data(), d4l3x.data());
    std::vector<double> d4l3y(static_cast<uint32_t>((d4_len * l3y_len) << 1));
    int d4l3y_len = product_expansion_zeroelim(d4_len, d4.data(), l3y_len, l3y.data(), d4l3y.data());
    std::vector<double> d4l3z(static_cast<uint32_t>((d4_len * l3z_len) << 1));
    int d4l3z_len = product_expansion_zeroelim(d4_len, d4.data(), l3z_len, l3z.data(), d4l3z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(d4l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(d4l1x_len, d4l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(d4l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(d4l1y_len, d4l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(d4l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(d4l1z_len, d4l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(d4l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(d4l2x_len, d4l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(d4l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(d4l2y_len, d4l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(d4l2z_len + d2p4z_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(d4l2z_len, d4l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    std::vector<double> p3p4x(static_cast<uint32_t>(d4l3x_len + d3p4x_len));
    int p3p4x_len = fast_expansion_diff_zeroelim(d4l3x_len, d4l3x.data(), d3p4x_len, d3p4x.data(), p3p4x.data());
    std::vector<double> p3p4y(static_cast<uint32_t>(d4l3y_len + d3p4y_len));
    int p3p4y_len = fast_expansion_diff_zeroelim(d4l3y_len, d4l3y.data(), d3p4y_len, d3p4y.data(), p3p4y.data());
    std::vector<double> p3p4z(static_cast<uint32_t>(d4l3z_len + d3p4z_len));
    int p3p4z_len = fast_expansion_diff_zeroelim(d4l3z_len, d4l3z.data(), d3p4z_len, d3p4z.data(), p3p4z.data());
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>((m01_len * p3p4z_len) << 1));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), p3p4z_len, p3p4z.data(), mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>((m02_len * p3p4y_len) << 1));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), p3p4y_len, p3p4y.data(), mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>((m12_len * p3p4x_len) << 1));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), p3p4x_len, p3p4x.data(), mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient3d_LLLT(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointLPI& p3, const ImplicitPointTPI& p4
) {
    int ret = orient3d_LLLT_interval(p1, p2, p3, p4);
    if (ret != 0) {
        return ret;
    }
    return orient3d_LLLT_exact(p1, p2, p3, p4);
}

inline int orient3d_LLTT_interval(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointTPI& p3, const ImplicitPointTPI& p4
) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda() || !p3.get_interval_lambda() ||
        !p4.get_interval_lambda()) {
        return IPSign::ZERO;
    }

    IntervalNumber d1p4x(p1.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d3p4x(p3.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d3p4y(p3.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d3p4z(p3.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d4l1x(p4.dfilter[3] * p1.dfilter[0]);
    IntervalNumber d4l1y(p4.dfilter[3] * p1.dfilter[1]);
    IntervalNumber d4l1z(p4.dfilter[3] * p1.dfilter[2]);
    IntervalNumber d4l2x(p4.dfilter[3] * p2.dfilter[0]);
    IntervalNumber d4l2y(p4.dfilter[3] * p2.dfilter[1]);
    IntervalNumber d4l2z(p4.dfilter[3] * p2.dfilter[2]);
    IntervalNumber d4l3x(p4.dfilter[3] * p3.dfilter[0]);
    IntervalNumber d4l3y(p4.dfilter[3] * p3.dfilter[1]);
    IntervalNumber d4l3z(p4.dfilter[3] * p3.dfilter[2]);
    IntervalNumber p1p4x(d4l1x - d1p4x);
    IntervalNumber p1p4y(d4l1y - d1p4y);
    IntervalNumber p1p4z(d4l1z - d1p4z);
    IntervalNumber p2p4x(d4l2x - d2p4x);
    IntervalNumber p2p4y(d4l2y - d2p4y);
    IntervalNumber p2p4z(d4l2z - d2p4z);
    IntervalNumber p3p4x(d4l3x - d3p4x);
    IntervalNumber p3p4y(d4l3y - d3p4y);
    IntervalNumber p3p4z(d4l3z - d3p4z);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);
    return m012.sign();
}

inline int orient3d_LLTT_exact(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointTPI& p3, const ImplicitPointTPI& p4
) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, l3x, l3y, l3z, l4x, l4y, l4z, d1, d2, d3, d4;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, l3x_len, l3y_len, l3z_len, l4x_len, l4y_len, l4z_len,
        d1_len, d2_len, d3_len, d4_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    p3.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    p4.get_exact_lambda(l4x, l4x_len, l4y, l4y_len, l4z, l4z_len, d4, d4_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0 || d3.data()[d3_len - 1] == 0.0 ||
        d4.data()[d4_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    std::vector<double> d1p4x(static_cast<uint32_t>((d1_len * l4x_len) << 1));
    int d1p4x_len = product_expansion_zeroelim(d1_len, d1.data(), l4x_len, l4x.data(), d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>((d1_len * l4y_len) << 1));
    int d1p4y_len = product_expansion_zeroelim(d1_len, d1.data(), l4y_len, l4y.data(), d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>((d1_len * l4z_len) << 1));
    int d1p4z_len = product_expansion_zeroelim(d1_len, d1.data(), l4z_len, l4z.data(), d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>((d2_len * l4x_len) << 1));
    int d2p4x_len = product_expansion_zeroelim(d2_len, d2.data(), l4x_len, l4x.data(), d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>((d2_len * l4y_len) << 1));
    int d2p4y_len = product_expansion_zeroelim(d2_len, d2.data(), l4y_len, l4y.data(), d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>((d2_len * l4z_len) << 1));
    int d2p4z_len = product_expansion_zeroelim(d2_len, d2.data(), l4z_len, l4z.data(), d2p4z.data());
    std::vector<double> d3p4x(static_cast<uint32_t>((d3_len * l4x_len) << 1));
    int d3p4x_len = product_expansion_zeroelim(d3_len, d3.data(), l4x_len, l4x.data(), d3p4x.data());
    std::vector<double> d3p4y(static_cast<uint32_t>((d3_len * l4y_len) << 1));
    int d3p4y_len = product_expansion_zeroelim(d3_len, d3.data(), l4y_len, l4y.data(), d3p4y.data());
    std::vector<double> d3p4z(static_cast<uint32_t>((d3_len * l4z_len) << 1));
    int d3p4z_len = product_expansion_zeroelim(d3_len, d3.data(), l4z_len, l4z.data(), d3p4z.data());
    std::vector<double> d4l1x(static_cast<uint32_t>((d4_len * l1x_len) << 1));
    int d4l1x_len = product_expansion_zeroelim(d4_len, d4.data(), l1x_len, l1x.data(), d4l1x.data());
    std::vector<double> d4l1y(static_cast<uint32_t>((d4_len * l1y_len) << 1));
    int d4l1y_len = product_expansion_zeroelim(d4_len, d4.data(), l1y_len, l1y.data(), d4l1y.data());
    std::vector<double> d4l1z(static_cast<uint32_t>((d4_len * l1z_len) << 1));
    int d4l1z_len = product_expansion_zeroelim(d4_len, d4.data(), l1z_len, l1z.data(), d4l1z.data());
    std::vector<double> d4l2x(static_cast<uint32_t>((d4_len * l2x_len) << 1));
    int d4l2x_len = product_expansion_zeroelim(d4_len, d4.data(), l2x_len, l2x.data(), d4l2x.data());
    std::vector<double> d4l2y(static_cast<uint32_t>((d4_len * l2y_len) << 1));
    int d4l2y_len = product_expansion_zeroelim(d4_len, d4.data(), l2y_len, l2y.data(), d4l2y.data());
    std::vector<double> d4l2z(static_cast<uint32_t>((d4_len * l2z_len) << 1));
    int d4l2z_len = product_expansion_zeroelim(d4_len, d4.data(), l2z_len, l2z.data(), d4l2z.data());
    std::vector<double> d4l3x(static_cast<uint32_t>((d4_len * l3x_len) << 1));
    int d4l3x_len = product_expansion_zeroelim(d4_len, d4.data(), l3x_len, l3x.data(), d4l3x.data());
    std::vector<double> d4l3y(static_cast<uint32_t>((d4_len * l3y_len) << 1));
    int d4l3y_len = product_expansion_zeroelim(d4_len, d4.data(), l3y_len, l3y.data(), d4l3y.data());
    std::vector<double> d4l3z(static_cast<uint32_t>((d4_len * l3z_len) << 1));
    int d4l3z_len = product_expansion_zeroelim(d4_len, d4.data(), l3z_len, l3z.data(), d4l3z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(d4l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(d4l1x_len, d4l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(d4l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(d4l1y_len, d4l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(d4l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(d4l1z_len, d4l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(d4l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(d4l2x_len, d4l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(d4l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(d4l2y_len, d4l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(d4l2z_len + d2p4z_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(d4l2z_len, d4l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    std::vector<double> p3p4x(static_cast<uint32_t>(d4l3x_len + d3p4x_len));
    int p3p4x_len = fast_expansion_diff_zeroelim(d4l3x_len, d4l3x.data(), d3p4x_len, d3p4x.data(), p3p4x.data());
    std::vector<double> p3p4y(static_cast<uint32_t>(d4l3y_len + d3p4y_len));
    int p3p4y_len = fast_expansion_diff_zeroelim(d4l3y_len, d4l3y.data(), d3p4y_len, d3p4y.data(), p3p4y.data());
    std::vector<double> p3p4z(static_cast<uint32_t>(d4l3z_len + d3p4z_len));
    int p3p4z_len = fast_expansion_diff_zeroelim(d4l3z_len, d4l3z.data(), d3p4z_len, d3p4z.data(), p3p4z.data());
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>((m01_len * p3p4z_len) << 1));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), p3p4z_len, p3p4z.data(), mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>((m02_len * p3p4y_len) << 1));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), p3p4y_len, p3p4y.data(), mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>((m12_len * p3p4x_len) << 1));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), p3p4x_len, p3p4x.data(), mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient3d_LLTT(
    const ImplicitPointLPI& p1, const ImplicitPointLPI& p2, const ImplicitPointTPI& p3, const ImplicitPointTPI& p4
) {
    int ret = orient3d_LLTT_interval(p1, p2, p3, p4);
    if (ret != 0) {
        return ret;
    }
    return orient3d_LLTT_exact(p1, p2, p3, p4);
}

inline int orient3d_LTTT_interval(
    const ImplicitPointLPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const ImplicitPointTPI& p4
) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda() || !p3.get_interval_lambda() ||
        !p4.get_interval_lambda()) {
        return IPSign::ZERO;
    }

    IntervalNumber d1p4x(p1.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d3p4x(p3.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d3p4y(p3.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d3p4z(p3.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d4l1x(p4.dfilter[3] * p1.dfilter[0]);
    IntervalNumber d4l1y(p4.dfilter[3] * p1.dfilter[1]);
    IntervalNumber d4l1z(p4.dfilter[3] * p1.dfilter[2]);
    IntervalNumber d4l2x(p4.dfilter[3] * p2.dfilter[0]);
    IntervalNumber d4l2y(p4.dfilter[3] * p2.dfilter[1]);
    IntervalNumber d4l2z(p4.dfilter[3] * p2.dfilter[2]);
    IntervalNumber d4l3x(p4.dfilter[3] * p3.dfilter[0]);
    IntervalNumber d4l3y(p4.dfilter[3] * p3.dfilter[1]);
    IntervalNumber d4l3z(p4.dfilter[3] * p3.dfilter[2]);
    IntervalNumber p1p4x(d4l1x - d1p4x);
    IntervalNumber p1p4y(d4l1y - d1p4y);
    IntervalNumber p1p4z(d4l1z - d1p4z);
    IntervalNumber p2p4x(d4l2x - d2p4x);
    IntervalNumber p2p4y(d4l2y - d2p4y);
    IntervalNumber p2p4z(d4l2z - d2p4z);
    IntervalNumber p3p4x(d4l3x - d3p4x);
    IntervalNumber p3p4y(d4l3y - d3p4y);
    IntervalNumber p3p4z(d4l3z - d3p4z);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);
    return m012.sign();
}

inline int orient3d_LTTT_exact(
    const ImplicitPointLPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const ImplicitPointTPI& p4
) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, l3x, l3y, l3z, l4x, l4y, l4z, d1, d2, d3, d4;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, l3x_len, l3y_len, l3z_len, l4x_len, l4y_len, l4z_len,
        d1_len, d2_len, d3_len, d4_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    p3.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    p4.get_exact_lambda(l4x, l4x_len, l4y, l4y_len, l4z, l4z_len, d4, d4_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0 || d3.data()[d3_len - 1] == 0.0 ||
        d4.data()[d4_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    std::vector<double> d1p4x(static_cast<uint32_t>((d1_len * l4x_len) << 1));
    int d1p4x_len = product_expansion_zeroelim(d1_len, d1.data(), l4x_len, l4x.data(), d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>((d1_len * l4y_len) << 1));
    int d1p4y_len = product_expansion_zeroelim(d1_len, d1.data(), l4y_len, l4y.data(), d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>((d1_len * l4z_len) << 1));
    int d1p4z_len = product_expansion_zeroelim(d1_len, d1.data(), l4z_len, l4z.data(), d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>((d2_len * l4x_len) << 1));
    int d2p4x_len = product_expansion_zeroelim(d2_len, d2.data(), l4x_len, l4x.data(), d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>((d2_len * l4y_len) << 1));
    int d2p4y_len = product_expansion_zeroelim(d2_len, d2.data(), l4y_len, l4y.data(), d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>((d2_len * l4z_len) << 1));
    int d2p4z_len = product_expansion_zeroelim(d2_len, d2.data(), l4z_len, l4z.data(), d2p4z.data());
    std::vector<double> d3p4x(static_cast<uint32_t>((d3_len * l4x_len) << 1));
    int d3p4x_len = product_expansion_zeroelim(d3_len, d3.data(), l4x_len, l4x.data(), d3p4x.data());
    std::vector<double> d3p4y(static_cast<uint32_t>((d3_len * l4y_len) << 1));
    int d3p4y_len = product_expansion_zeroelim(d3_len, d3.data(), l4y_len, l4y.data(), d3p4y.data());
    std::vector<double> d3p4z(static_cast<uint32_t>((d3_len * l4z_len) << 1));
    int d3p4z_len = product_expansion_zeroelim(d3_len, d3.data(), l4z_len, l4z.data(), d3p4z.data());
    std::vector<double> d4l1x(static_cast<uint32_t>((d4_len * l1x_len) << 1));
    int d4l1x_len = product_expansion_zeroelim(d4_len, d4.data(), l1x_len, l1x.data(), d4l1x.data());
    std::vector<double> d4l1y(static_cast<uint32_t>((d4_len * l1y_len) << 1));
    int d4l1y_len = product_expansion_zeroelim(d4_len, d4.data(), l1y_len, l1y.data(), d4l1y.data());
    std::vector<double> d4l1z(static_cast<uint32_t>((d4_len * l1z_len) << 1));
    int d4l1z_len = product_expansion_zeroelim(d4_len, d4.data(), l1z_len, l1z.data(), d4l1z.data());
    std::vector<double> d4l2x(static_cast<uint32_t>((d4_len * l2x_len) << 1));
    int d4l2x_len = product_expansion_zeroelim(d4_len, d4.data(), l2x_len, l2x.data(), d4l2x.data());
    std::vector<double> d4l2y(static_cast<uint32_t>((d4_len * l2y_len) << 1));
    int d4l2y_len = product_expansion_zeroelim(d4_len, d4.data(), l2y_len, l2y.data(), d4l2y.data());
    std::vector<double> d4l2z(static_cast<uint32_t>((d4_len * l2z_len) << 1));
    int d4l2z_len = product_expansion_zeroelim(d4_len, d4.data(), l2z_len, l2z.data(), d4l2z.data());
    std::vector<double> d4l3x(static_cast<uint32_t>((d4_len * l3x_len) << 1));
    int d4l3x_len = product_expansion_zeroelim(d4_len, d4.data(), l3x_len, l3x.data(), d4l3x.data());
    std::vector<double> d4l3y(static_cast<uint32_t>((d4_len * l3y_len) << 1));
    int d4l3y_len = product_expansion_zeroelim(d4_len, d4.data(), l3y_len, l3y.data(), d4l3y.data());
    std::vector<double> d4l3z(static_cast<uint32_t>((d4_len * l3z_len) << 1));
    int d4l3z_len = product_expansion_zeroelim(d4_len, d4.data(), l3z_len, l3z.data(), d4l3z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(d4l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(d4l1x_len, d4l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(d4l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(d4l1y_len, d4l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(d4l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(d4l1z_len, d4l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(d4l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(d4l2x_len, d4l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(d4l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(d4l2y_len, d4l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(d4l2z_len + d2p4z_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(d4l2z_len, d4l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    std::vector<double> p3p4x(static_cast<uint32_t>(d4l3x_len + d3p4x_len));
    int p3p4x_len = fast_expansion_diff_zeroelim(d4l3x_len, d4l3x.data(), d3p4x_len, d3p4x.data(), p3p4x.data());
    std::vector<double> p3p4y(static_cast<uint32_t>(d4l3y_len + d3p4y_len));
    int p3p4y_len = fast_expansion_diff_zeroelim(d4l3y_len, d4l3y.data(), d3p4y_len, d3p4y.data(), p3p4y.data());
    std::vector<double> p3p4z(static_cast<uint32_t>(d4l3z_len + d3p4z_len));
    int p3p4z_len = fast_expansion_diff_zeroelim(d4l3z_len, d4l3z.data(), d3p4z_len, d3p4z.data(), p3p4z.data());
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>((m01_len * p3p4z_len) << 1));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), p3p4z_len, p3p4z.data(), mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>((m02_len * p3p4y_len) << 1));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), p3p4y_len, p3p4y.data(), mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>((m12_len * p3p4x_len) << 1));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), p3p4x_len, p3p4x.data(), mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient3d_LTTT(
    const ImplicitPointLPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const ImplicitPointTPI& p4
) {
    int ret = orient3d_LTTT_interval(p1, p2, p3, p4);
    if (ret != 0) {
        return ret;
    }
    return orient3d_LTTT_exact(p1, p2, p3, p4);
}

inline int orient3d_TTTT_interval(
    const ImplicitPointTPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const ImplicitPointTPI& p4
) {
    if (!p1.get_interval_lambda() || !p2.get_interval_lambda() || !p3.get_interval_lambda() ||
        !p4.get_interval_lambda()) {
        return IPSign::ZERO;
    }

    IntervalNumber d1p4x(p1.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d1p4y(p1.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d1p4z(p1.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d2p4x(p2.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d2p4y(p2.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d2p4z(p2.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d3p4x(p3.dfilter[3] * p4.dfilter[0]);
    IntervalNumber d3p4y(p3.dfilter[3] * p4.dfilter[1]);
    IntervalNumber d3p4z(p3.dfilter[3] * p4.dfilter[2]);
    IntervalNumber d4l1x(p4.dfilter[3] * p1.dfilter[0]);
    IntervalNumber d4l1y(p4.dfilter[3] * p1.dfilter[1]);
    IntervalNumber d4l1z(p4.dfilter[3] * p1.dfilter[2]);
    IntervalNumber d4l2x(p4.dfilter[3] * p2.dfilter[0]);
    IntervalNumber d4l2y(p4.dfilter[3] * p2.dfilter[1]);
    IntervalNumber d4l2z(p4.dfilter[3] * p2.dfilter[2]);
    IntervalNumber d4l3x(p4.dfilter[3] * p3.dfilter[0]);
    IntervalNumber d4l3y(p4.dfilter[3] * p3.dfilter[1]);
    IntervalNumber d4l3z(p4.dfilter[3] * p3.dfilter[2]);
    IntervalNumber p1p4x(d4l1x - d1p4x);
    IntervalNumber p1p4y(d4l1y - d1p4y);
    IntervalNumber p1p4z(d4l1z - d1p4z);
    IntervalNumber p2p4x(d4l2x - d2p4x);
    IntervalNumber p2p4y(d4l2y - d2p4y);
    IntervalNumber p2p4z(d4l2z - d2p4z);
    IntervalNumber p3p4x(d4l3x - d3p4x);
    IntervalNumber p3p4y(d4l3y - d3p4y);
    IntervalNumber p3p4z(d4l3z - d3p4z);
    IntervalNumber tmc_a(p1p4x * p2p4y);
    IntervalNumber tmc_b(p1p4y * p2p4x);
    IntervalNumber m01(tmc_a - tmc_b);
    IntervalNumber tmi_a(p1p4x * p2p4z);
    IntervalNumber tmi_b(p1p4z * p2p4x);
    IntervalNumber m02(tmi_a - tmi_b);
    IntervalNumber tma_a(p1p4y * p2p4z);
    IntervalNumber tma_b(p1p4z * p2p4y);
    IntervalNumber m12(tma_a - tma_b);
    IntervalNumber mt1(m01 * p3p4z);
    IntervalNumber mt2(m02 * p3p4y);
    IntervalNumber mt3(m12 * p3p4x);
    IntervalNumber mtt(mt2 - mt1);
    IntervalNumber m012(mtt - mt3);
    return m012.sign();
}

inline int orient3d_TTTT_exact(
    const ImplicitPointTPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const ImplicitPointTPI& p4
) {
    std::vector<double> l1x, l1y, l1z, l2x, l2y, l2z, l3x, l3y, l3z, l4x, l4y, l4z, d1, d2, d3, d4;
    int l1x_len, l1y_len, l1z_len, l2x_len, l2y_len, l2z_len, l3x_len, l3y_len, l3z_len, l4x_len, l4y_len, l4z_len,
        d1_len, d2_len, d3_len, d4_len;
    p1.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    p2.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    p3.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    p4.get_exact_lambda(l4x, l4x_len, l4y, l4y_len, l4z, l4z_len, d4, d4_len);
    if (d1.data()[d1_len - 1] == 0.0 || d2.data()[d2_len - 1] == 0.0 || d3.data()[d3_len - 1] == 0.0 ||
        d4.data()[d4_len - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }

    std::vector<double> d1p4x(static_cast<uint32_t>((d1_len * l4x_len) << 1));
    int d1p4x_len = product_expansion_zeroelim(d1_len, d1.data(), l4x_len, l4x.data(), d1p4x.data());
    std::vector<double> d1p4y(static_cast<uint32_t>((d1_len * l4y_len) << 1));
    int d1p4y_len = product_expansion_zeroelim(d1_len, d1.data(), l4y_len, l4y.data(), d1p4y.data());
    std::vector<double> d1p4z(static_cast<uint32_t>((d1_len * l4z_len) << 1));
    int d1p4z_len = product_expansion_zeroelim(d1_len, d1.data(), l4z_len, l4z.data(), d1p4z.data());
    std::vector<double> d2p4x(static_cast<uint32_t>((d2_len * l4x_len) << 1));
    int d2p4x_len = product_expansion_zeroelim(d2_len, d2.data(), l4x_len, l4x.data(), d2p4x.data());
    std::vector<double> d2p4y(static_cast<uint32_t>((d2_len * l4y_len) << 1));
    int d2p4y_len = product_expansion_zeroelim(d2_len, d2.data(), l4y_len, l4y.data(), d2p4y.data());
    std::vector<double> d2p4z(static_cast<uint32_t>((d2_len * l4z_len) << 1));
    int d2p4z_len = product_expansion_zeroelim(d2_len, d2.data(), l4z_len, l4z.data(), d2p4z.data());
    std::vector<double> d3p4x(static_cast<uint32_t>((d3_len * l4x_len) << 1));
    int d3p4x_len = product_expansion_zeroelim(d3_len, d3.data(), l4x_len, l4x.data(), d3p4x.data());
    std::vector<double> d3p4y(static_cast<uint32_t>((d3_len * l4y_len) << 1));
    int d3p4y_len = product_expansion_zeroelim(d3_len, d3.data(), l4y_len, l4y.data(), d3p4y.data());
    std::vector<double> d3p4z(static_cast<uint32_t>((d3_len * l4z_len) << 1));
    int d3p4z_len = product_expansion_zeroelim(d3_len, d3.data(), l4z_len, l4z.data(), d3p4z.data());
    std::vector<double> d4l1x(static_cast<uint32_t>((d4_len * l1x_len) << 1));
    int d4l1x_len = product_expansion_zeroelim(d4_len, d4.data(), l1x_len, l1x.data(), d4l1x.data());
    std::vector<double> d4l1y(static_cast<uint32_t>((d4_len * l1y_len) << 1));
    int d4l1y_len = product_expansion_zeroelim(d4_len, d4.data(), l1y_len, l1y.data(), d4l1y.data());
    std::vector<double> d4l1z(static_cast<uint32_t>((d4_len * l1z_len) << 1));
    int d4l1z_len = product_expansion_zeroelim(d4_len, d4.data(), l1z_len, l1z.data(), d4l1z.data());
    std::vector<double> d4l2x(static_cast<uint32_t>((d4_len * l2x_len) << 1));
    int d4l2x_len = product_expansion_zeroelim(d4_len, d4.data(), l2x_len, l2x.data(), d4l2x.data());
    std::vector<double> d4l2y(static_cast<uint32_t>((d4_len * l2y_len) << 1));
    int d4l2y_len = product_expansion_zeroelim(d4_len, d4.data(), l2y_len, l2y.data(), d4l2y.data());
    std::vector<double> d4l2z(static_cast<uint32_t>((d4_len * l2z_len) << 1));
    int d4l2z_len = product_expansion_zeroelim(d4_len, d4.data(), l2z_len, l2z.data(), d4l2z.data());
    std::vector<double> d4l3x(static_cast<uint32_t>((d4_len * l3x_len) << 1));
    int d4l3x_len = product_expansion_zeroelim(d4_len, d4.data(), l3x_len, l3x.data(), d4l3x.data());
    std::vector<double> d4l3y(static_cast<uint32_t>((d4_len * l3y_len) << 1));
    int d4l3y_len = product_expansion_zeroelim(d4_len, d4.data(), l3y_len, l3y.data(), d4l3y.data());
    std::vector<double> d4l3z(static_cast<uint32_t>((d4_len * l3z_len) << 1));
    int d4l3z_len = product_expansion_zeroelim(d4_len, d4.data(), l3z_len, l3z.data(), d4l3z.data());
    std::vector<double> p1p4x(static_cast<uint32_t>(d4l1x_len + d1p4x_len));
    int p1p4x_len = fast_expansion_diff_zeroelim(d4l1x_len, d4l1x.data(), d1p4x_len, d1p4x.data(), p1p4x.data());
    std::vector<double> p1p4y(static_cast<uint32_t>(d4l1y_len + d1p4y_len));
    int p1p4y_len = fast_expansion_diff_zeroelim(d4l1y_len, d4l1y.data(), d1p4y_len, d1p4y.data(), p1p4y.data());
    std::vector<double> p1p4z(static_cast<uint32_t>(d4l1z_len + d1p4z_len));
    int p1p4z_len = fast_expansion_diff_zeroelim(d4l1z_len, d4l1z.data(), d1p4z_len, d1p4z.data(), p1p4z.data());
    std::vector<double> p2p4x(static_cast<uint32_t>(d4l2x_len + d2p4x_len));
    int p2p4x_len = fast_expansion_diff_zeroelim(d4l2x_len, d4l2x.data(), d2p4x_len, d2p4x.data(), p2p4x.data());
    std::vector<double> p2p4y(static_cast<uint32_t>(d4l2y_len + d2p4y_len));
    int p2p4y_len = fast_expansion_diff_zeroelim(d4l2y_len, d4l2y.data(), d2p4y_len, d2p4y.data(), p2p4y.data());
    std::vector<double> p2p4z(static_cast<uint32_t>(d4l2z_len + d2p4z_len));
    int p2p4z_len = fast_expansion_diff_zeroelim(d4l2z_len, d4l2z.data(), d2p4z_len, d2p4z.data(), p2p4z.data());
    std::vector<double> p3p4x(static_cast<uint32_t>(d4l3x_len + d3p4x_len));
    int p3p4x_len = fast_expansion_diff_zeroelim(d4l3x_len, d4l3x.data(), d3p4x_len, d3p4x.data(), p3p4x.data());
    std::vector<double> p3p4y(static_cast<uint32_t>(d4l3y_len + d3p4y_len));
    int p3p4y_len = fast_expansion_diff_zeroelim(d4l3y_len, d4l3y.data(), d3p4y_len, d3p4y.data(), p3p4y.data());
    std::vector<double> p3p4z(static_cast<uint32_t>(d4l3z_len + d3p4z_len));
    int p3p4z_len = fast_expansion_diff_zeroelim(d4l3z_len, d4l3z.data(), d3p4z_len, d3p4z.data(), p3p4z.data());
    std::vector<double> tmc_a(static_cast<uint32_t>((p1p4x_len * p2p4y_len) << 1));
    int tmc_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4y_len, p2p4y.data(), tmc_a.data());
    std::vector<double> tmc_b(static_cast<uint32_t>((p1p4y_len * p2p4x_len) << 1));
    int tmc_b_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4x_len, p2p4x.data(), tmc_b.data());
    std::vector<double> m01(static_cast<uint32_t>(tmc_a_len + tmc_b_len));
    int m01_len = fast_expansion_diff_zeroelim(tmc_a_len, tmc_a.data(), tmc_b_len, tmc_b.data(), m01.data());
    std::vector<double> tmi_a(static_cast<uint32_t>((p1p4x_len * p2p4z_len) << 1));
    int tmi_a_len = product_expansion_zeroelim(p1p4x_len, p1p4x.data(), p2p4z_len, p2p4z.data(), tmi_a.data());
    std::vector<double> tmi_b(static_cast<uint32_t>((p1p4z_len * p2p4x_len) << 1));
    int tmi_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4x_len, p2p4x.data(), tmi_b.data());
    std::vector<double> m02(static_cast<uint32_t>(tmi_a_len + tmi_b_len));
    int m02_len = fast_expansion_diff_zeroelim(tmi_a_len, tmi_a.data(), tmi_b_len, tmi_b.data(), m02.data());
    std::vector<double> tma_a(static_cast<uint32_t>((p1p4y_len * p2p4z_len) << 1));
    int tma_a_len = product_expansion_zeroelim(p1p4y_len, p1p4y.data(), p2p4z_len, p2p4z.data(), tma_a.data());
    std::vector<double> tma_b(static_cast<uint32_t>((p1p4z_len * p2p4y_len) << 1));
    int tma_b_len = product_expansion_zeroelim(p1p4z_len, p1p4z.data(), p2p4y_len, p2p4y.data(), tma_b.data());
    std::vector<double> m12(static_cast<uint32_t>(tma_a_len + tma_b_len));
    int m12_len = fast_expansion_diff_zeroelim(tma_a_len, tma_a.data(), tma_b_len, tma_b.data(), m12.data());
    std::vector<double> mt1(static_cast<uint32_t>((m01_len * p3p4z_len) << 1));
    int mt1_len = product_expansion_zeroelim(m01_len, m01.data(), p3p4z_len, p3p4z.data(), mt1.data());
    std::vector<double> mt2(static_cast<uint32_t>((m02_len * p3p4y_len) << 1));
    int mt2_len = product_expansion_zeroelim(m02_len, m02.data(), p3p4y_len, p3p4y.data(), mt2.data());
    std::vector<double> mt3(static_cast<uint32_t>((m12_len * p3p4x_len) << 1));
    int mt3_len = product_expansion_zeroelim(m12_len, m12.data(), p3p4x_len, p3p4x.data(), mt3.data());
    std::vector<double> mtt(static_cast<uint32_t>(mt2_len + mt1_len));
    int mtt_len = fast_expansion_diff_zeroelim(mt2_len, mt2.data(), mt1_len, mt1.data(), mtt.data());
    std::vector<double> m012(static_cast<uint32_t>(mtt_len + mt3_len));
    int m012_len = fast_expansion_diff_zeroelim(mtt_len, mtt.data(), mt3_len, mt3.data(), m012.data());
    if (m012.data()[m012_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (m012.data()[m012_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient3d_TTTT(
    const ImplicitPointTPI& p1, const ImplicitPointTPI& p2, const ImplicitPointTPI& p3, const ImplicitPointTPI& p4
) {
    int ret = orient3d_TTTT_interval(p1, p2, p3, p4);
    if (ret != 0) {
        return ret;
    }
    return orient3d_TTTT_exact(p1, p2, p3, p4);
}

int GenericPoint3D::orient3d(
    const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c, const GenericPoint3D& d
) {
    const int val = a.get_type() * 27 + b.get_type() * 9 + c.get_type() * 3 + d.get_type();
    switch (val) {
    case 0: {
        const double ret =
            ::orient3d(a.to_explicit().ptr(), b.to_explicit().ptr(), c.to_explicit().ptr(), d.to_explicit().ptr());
        return (ret > 0.0) - (ret < 0.0);
    }
    case 1: // EEEL
        return -orient3d_LEEE(d.to_lpi(), a.to_explicit().ptr(), c.to_explicit().ptr(), b.to_explicit().ptr());
    case 2: // EEET
        return -orient3d_TEEE(d.to_tpi(), a.to_explicit().ptr(), c.to_explicit().ptr(), b.to_explicit().ptr());
    case 3: // EELE
        return -orient3d_LEEE(c.to_lpi(), d.to_explicit().ptr(), a.to_explicit().ptr(), b.to_explicit().ptr());
    case 4: // EELL
        return -orient3d_LLEE(c.to_lpi(), d.to_lpi(), a.to_explicit().ptr(), b.to_explicit().ptr());
    case 5: // EELT
        return -orient3d_LTEE(c.to_lpi(), d.to_tpi(), a.to_explicit().ptr(), b.to_explicit().ptr());
    case 6: // EETE
        return -orient3d_TEEE(c.to_tpi(), d.to_explicit().ptr(), a.to_explicit().ptr(), b.to_explicit().ptr());
    case 7: // EETL
        return -orient3d_LTEE(d.to_lpi(), c.to_tpi(), b.to_explicit().ptr(), a.to_explicit().ptr());
    case 8: // EETT
        return -orient3d_TTEE(c.to_tpi(), d.to_tpi(), a.to_explicit().ptr(), b.to_explicit().ptr());
    case 9: // ELEE
        return -orient3d_LEEE(b.to_lpi(), c.to_explicit().ptr(), a.to_explicit().ptr(), d.to_explicit().ptr());
    case 10: // ELEL
        return -orient3d_LLEE(d.to_lpi(), b.to_lpi(), a.to_explicit().ptr(), c.to_explicit().ptr());
    case 11: // ELET
        return -orient3d_LTEE(b.to_lpi(), d.to_tpi(), c.to_explicit().ptr(), a.to_explicit().ptr());
    case 12: // ELLE
        return -orient3d_LLEE(b.to_lpi(), c.to_lpi(), a.to_explicit().ptr(), d.to_explicit().ptr());
    case 13: // ELLL
        return -orient3d_LLLE(b.to_lpi(), d.to_lpi(), c.to_lpi(), a.to_explicit().ptr());
    case 14: // ELLT
        return -orient3d_LLTE(c.to_lpi(), b.to_lpi(), d.to_tpi(), a.to_explicit().ptr());
    case 15: // ELTE
        return -orient3d_LTEE(b.to_lpi(), c.to_tpi(), a.to_explicit().ptr(), d.to_explicit().ptr());
    case 16: // ELTL
        return -orient3d_LLTE(b.to_lpi(), d.to_lpi(), c.to_tpi(), a.to_explicit().ptr());
    case 17: // ELTT
        return -orient3d_LTTE(b.to_lpi(), d.to_tpi(), c.to_tpi(), a.to_explicit().ptr());
    case 18: // ETEE
        return -orient3d_TEEE(b.to_tpi(), c.to_explicit().ptr(), a.to_explicit().ptr(), d.to_explicit().ptr());
    case 19: // ETEL
        return -orient3d_LTEE(d.to_lpi(), b.to_tpi(), a.to_explicit().ptr(), c.to_explicit().ptr());
    case 20: // ETET
        return -orient3d_TTEE(d.to_tpi(), b.to_tpi(), a.to_explicit().ptr(), c.to_explicit().ptr());
    case 21: // ETLE
        return -orient3d_LTEE(c.to_lpi(), b.to_tpi(), d.to_explicit().ptr(), a.to_explicit().ptr());
    case 22: // ETLL
        return -orient3d_LLTE(d.to_lpi(), c.to_lpi(), b.to_tpi(), a.to_explicit().ptr());
    case 23: // ETLT
        return -orient3d_LTTE(c.to_lpi(), b.to_tpi(), d.to_tpi(), a.to_explicit().ptr());
    case 24: // ETTE
        return -orient3d_TTEE(b.to_tpi(), c.to_tpi(), a.to_explicit().ptr(), d.to_explicit().ptr());
    case 25: // ETTL
        return -orient3d_LTTE(d.to_lpi(), c.to_tpi(), b.to_tpi(), a.to_explicit().ptr());
    case 26: // ETTT
        return -orient3d_TTTE(b.to_tpi(), d.to_tpi(), c.to_tpi(), a.to_explicit().ptr());
    case 27: // LEEE
        return -orient3d_LEEE(a.to_lpi(), b.to_explicit().ptr(), c.to_explicit().ptr(), d.to_explicit().ptr());
    case 28: // LEEL
        return -orient3d_LLEE(d.to_lpi(), a.to_lpi(), c.to_explicit().ptr(), b.to_explicit().ptr());
    case 29: // LEET
        return -orient3d_LTEE(a.to_lpi(), d.to_tpi(), b.to_explicit().ptr(), c.to_explicit().ptr());
    case 30: // LELE
        return -orient3d_LLEE(a.to_lpi(), c.to_lpi(), d.to_explicit().ptr(), b.to_explicit().ptr());
    case 31: // LELL
        return -orient3d_LLLE(a.to_lpi(), c.to_lpi(), d.to_lpi(), b.to_explicit().ptr());
    case 32: // LELT
        return -orient3d_LLTE(a.to_lpi(), c.to_lpi(), d.to_tpi(), b.to_explicit().ptr());
    case 33: // LETE
        return -orient3d_LTEE(a.to_lpi(), c.to_tpi(), d.to_explicit().ptr(), b.to_explicit().ptr());
    case 34: // LETL
        return -orient3d_LLTE(d.to_lpi(), a.to_lpi(), c.to_tpi(), b.to_explicit().ptr());
    case 35: // LETT
        return -orient3d_LTTE(a.to_lpi(), c.to_tpi(), d.to_tpi(), b.to_explicit().ptr());
    case 36: // LLEE
        return -orient3d_LLEE(a.to_lpi(), b.to_lpi(), c.to_explicit().ptr(), d.to_explicit().ptr());
    case 37: // LLEL
        return -orient3d_LLLE(d.to_lpi(), b.to_lpi(), a.to_lpi(), c.to_explicit().ptr());
    case 38: // LLET
        return -orient3d_LLTE(b.to_lpi(), a.to_lpi(), d.to_tpi(), c.to_explicit().ptr());
    case 39: // LLLE
        return -orient3d_LLLE(a.to_lpi(), b.to_lpi(), c.to_lpi(), d.to_explicit().ptr());
    case 40: // LLLL
        return -orient3d_LLLL(a.to_lpi(), b.to_lpi(), c.to_lpi(), d.to_lpi());
    case 41: // LLLT
        return -orient3d_LLLT(a.to_lpi(), b.to_lpi(), c.to_lpi(), d.to_tpi());
    case 42: // LLTE
        return -orient3d_LLTE(a.to_lpi(), b.to_lpi(), c.to_tpi(), d.to_explicit().ptr());
    case 43: // LLTL
        return -orient3d_LLLT(b.to_lpi(), a.to_lpi(), d.to_lpi(), c.to_tpi());
    case 44: // LLTT
        return -orient3d_LLTT(a.to_lpi(), b.to_lpi(), c.to_tpi(), d.to_tpi());
    case 45: // LTEE
        return -orient3d_LTEE(a.to_lpi(), b.to_tpi(), c.to_explicit().ptr(), d.to_explicit().ptr());
    case 46: // LTEL
        return -orient3d_LLTE(a.to_lpi(), d.to_lpi(), b.to_tpi(), c.to_explicit().ptr());
    case 47: // LTET
        return -orient3d_LTTE(a.to_lpi(), d.to_tpi(), b.to_tpi(), c.to_explicit().ptr());
    case 48: // LTLE
        return -orient3d_LLTE(c.to_lpi(), a.to_lpi(), b.to_tpi(), d.to_explicit().ptr());
    case 49: // LTLL
        return -orient3d_LLLT(c.to_lpi(), d.to_lpi(), a.to_lpi(), b.to_tpi());
    case 50: // LTLT
        return -orient3d_LLTT(c.to_lpi(), a.to_lpi(), b.to_tpi(), d.to_tpi());
    case 51: // LTTE
        return -orient3d_LTTE(d.to_lpi(), b.to_tpi(), c.to_tpi(), d.to_explicit().ptr());
    case 52: // LTTL
        return -orient3d_LLTT(a.to_lpi(), d.to_lpi(), b.to_tpi(), c.to_tpi());
    case 53: // LTTT
        return -orient3d_LTTT(a.to_lpi(), b.to_tpi(), c.to_tpi(), d.to_tpi());
    case 54: // TEEE
        return -orient3d_TEEE(a.to_tpi(), b.to_explicit().ptr(), c.to_explicit().ptr(), d.to_explicit().ptr());
    case 55: // TEEL
        return -orient3d_LTEE(d.to_lpi(), a.to_tpi(), c.to_explicit().ptr(), b.to_explicit().ptr());
    case 56: // TEET
        return -orient3d_TTEE(d.to_tpi(), a.to_tpi(), c.to_explicit().ptr(), b.to_explicit().ptr());
    case 57: // TELE
        return -orient3d_LTEE(c.to_lpi(), a.to_tpi(), b.to_explicit().ptr(), d.to_explicit().ptr());
    case 58: // TELL
        return -orient3d_LLTE(c.to_lpi(), d.to_lpi(), a.to_tpi(), b.to_explicit().ptr());
    case 59: // TELT
        return -orient3d_LTTE(c.to_lpi(), d.to_tpi(), a.to_tpi(), b.to_explicit().ptr());
    case 60: // TETE
        return -orient3d_TTEE(a.to_tpi(), c.to_tpi(), b.to_explicit().ptr(), d.to_explicit().ptr());
    case 61: // TETL
        return -orient3d_LTTE(c.to_lpi(), d.to_tpi(), a.to_tpi(), b.to_explicit().ptr());
    case 62: // TETT
        return -orient3d_TTTE(a.to_tpi(), c.to_tpi(), d.to_tpi(), b.to_explicit().ptr());
    case 63: // TLEE
        return -orient3d_LTEE(b.to_lpi(), a.to_tpi(), d.to_explicit().ptr(), c.to_explicit().ptr());
    case 64: // TLEL
        return -orient3d_LLTE(d.to_lpi(), b.to_lpi(), a.to_tpi(), c.to_explicit().ptr());
    case 65: // TLET
        return -orient3d_LTTE(b.to_lpi(), a.to_tpi(), d.to_tpi(), c.to_explicit().ptr());
    case 66: // TLLE
        return -orient3d_LLTE(b.to_lpi(), c.to_lpi(), a.to_tpi(), d.to_explicit().ptr());
    case 67: // TLLL
        return -orient3d_LLLT(d.to_lpi(), c.to_lpi(), b.to_lpi(), a.to_tpi());
    case 68: // TLLT
        return -orient3d_LLTT(b.to_lpi(), c.to_lpi(), a.to_tpi(), d.to_tpi());
    case 69: // TLTE
        return -orient3d_LTTE(b.to_lpi(), c.to_tpi(), a.to_tpi(), d.to_explicit().ptr());
    case 70: // TLTL
        return -orient3d_LLTT(b.to_lpi(), d.to_lpi(), c.to_tpi(), a.to_tpi());
    case 71: // TLTT
        return -orient3d_LTTT(b.to_lpi(), a.to_tpi(), d.to_tpi(), c.to_tpi());
    case 72: // TTEE
        return -orient3d_TTEE(a.to_tpi(), b.to_tpi(), c.to_explicit().ptr(), d.to_explicit().ptr());
    case 73: // TTEL
        return -orient3d_LTTE(d.to_lpi(), b.to_tpi(), a.to_tpi(), c.to_explicit().ptr());
    case 74: // TTET
        return -orient3d_TTTE(d.to_tpi(), b.to_tpi(), a.to_tpi(), c.to_explicit().ptr());
    case 75: // TTLE
        return -orient3d_LTTE(c.to_lpi(), a.to_tpi(), b.to_tpi(), d.to_explicit().ptr());
    case 76: // TTLL
        return -orient3d_LLTT(c.to_lpi(), d.to_lpi(), a.to_tpi(), b.to_tpi());
    case 77: // TTLT
        return -orient3d_LTTT(c.to_lpi(), d.to_tpi(), a.to_tpi(), b.to_tpi());
    case 78: // TTTE
        return -orient3d_TTTE(a.to_tpi(), b.to_tpi(), c.to_tpi(), d.to_explicit().ptr());
    case 79: // TTTL
        return -orient3d_LTTT(d.to_lpi(), c.to_tpi(), b.to_tpi(), a.to_tpi());
    case 80: // TTTT
        return -orient3d_TTTT(a.to_tpi(), b.to_tpi(), c.to_tpi(), d.to_tpi());
    default:
        return IPSign::UNDEFINED;
    }
}

inline int
orient_xy_lee_filtered(const double* l1, const double d1, const double* p2, const double* p3, double max_var) {
    double t1x = p2[1] - p3[1];
    double t1y = p3[0] - p2[0];
    double e2 = l1[0] * t1x;
    double e3 = l1[1] * t1y;
    double e = e2 + e3;
    double pr1 = p2[0] * p3[1];
    double pr2 = p2[1] * p3[0];
    double pr = pr1 - pr2;
    double dpr = d1 * pr;
    double det = dpr + e;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(p2[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p2[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(t1x)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(t1y)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= 4.752773695437811e-14;
    if (det > epsilon) return IPSign::POSITIVE;
    if (-det > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient_xy_iee_interval(
    const IntervalNumber& l1x, const IntervalNumber& l1y, const IntervalNumber& d1, const IntervalNumber& p2x,
    const IntervalNumber& p2y, const IntervalNumber& p3x, const IntervalNumber& p3y
) {
    IntervalNumber t1x(p2y - p3y);
    IntervalNumber t1y(p3x - p2x);
    IntervalNumber e2(l1x * t1x);
    IntervalNumber e3(l1y * t1y);
    IntervalNumber e(e2 + e3);
    IntervalNumber pr1(p2x * p3y);
    IntervalNumber pr2(p2y * p3x);
    IntervalNumber pr(pr1 - pr2);
    IntervalNumber dpr(d1 * pr);
    IntervalNumber det(dpr + e);
    return det.sign();
}

inline int orient_xy_iee_exact(
    const std::vector<double>& l1x, const int l1x_len, const std::vector<double>& l1y, const int l1y_len,
    const std::vector<double>& d1, const int d1_len, const double* p2, const double* p3
) {
    double t1x[2];
    two_diff(p2[1], p3[1], t1x);
    double t1y[2];
    two_diff(p3[0], p2[0], t1y);
    std::vector<double> e2(static_cast<uint32_t>(l1x_len << 2));
    int e2_len = product_expansion_zeroelim(l1x_len, l1x.data(), 2, t1x, e2.data());
    std::vector<double> e3(static_cast<uint32_t>(l1y_len << 2));
    int e3_len = product_expansion_zeroelim(l1y_len, l1y.data(), 2, t1y, e3.data());
    std::vector<double> e(static_cast<uint32_t>(e2_len + e3_len));
    int e_len = fast_expansion_sum_zeroelim(e2_len, e2.data(), e3_len, e3.data(), e.data());
    double pr1[2];
    two_prod(p2[0], p3[1], pr1);
    double pr2[2];
    two_prod(p2[1], p3[0], pr2);
    double pr[4];
    two_two_diff(pr1, pr2, pr);
    std::vector<double> dpr(static_cast<uint32_t>(d1_len << 3));
    int dpr_len = product_expansion_zeroelim(d1_len, d1.data(), 4, pr, dpr.data());
    std::vector<double> det(static_cast<uint32_t>(dpr_len + e_len));
    int det_len = fast_expansion_sum_zeroelim(dpr_len, dpr.data(), e_len, e.data(), det.data());
    if (det.data()[det_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (det.data()[det_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient_xy_lee(const ImplicitPointLPI& a, const ExplicitPoint3D& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_lee_filtered(a.ssfilter.data(), a.ssfilter[3], b.ptr(), c.ptr(), max_var)) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda()) {
        if ((ret = orient_xy_iee_interval(a.dfilter[0], a.dfilter[1], a.dfilter[3], b.x, b.y, c.x, c.y)) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1;
    int l1x_len, l1y_len, l1z_len, d1_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    return orient_xy_iee_exact(l1x, l1x_len, l1y, l1y_len, d1, d1_len, b.ptr(), c.ptr());
}

inline int orient_yz_lee(const ImplicitPointLPI& a, const ExplicitPoint3D& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_lee_filtered(a.ssfilter.data() + 1, a.ssfilter[3], b.ptr() + 1, c.ptr() + 1, max_var)) !=
            0.0)
            return ret;
    }
    if (a.get_interval_lambda()) {
        if ((ret = orient_xy_iee_interval(a.dfilter[1], a.dfilter[2], a.dfilter[3], b.y, b.z, c.y, c.z)) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1;
    int l1x_len, l1y_len, l1z_len, d1_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    return orient_xy_iee_exact(l1y, l1y_len, l1z, l1z_len, d1, d1_len, b.ptr() + 1, c.ptr() + 1);
}

inline int orient_zx_lee(const ImplicitPointLPI& a, const ExplicitPoint3D& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    const double p2[2]{b.z, b.x};
    const double p3[2]{c.z, c.x};
    if (a.get_filtered_lambda(max_var)) {
        const double l1[2]{a.ssfilter[2], a.ssfilter[0]};
        if ((ret = orient_xy_lee_filtered(l1, a.ssfilter[3], p2, p3, max_var)) != 0.0) return ret;
    }
    if (a.get_interval_lambda()) {
        if ((ret = orient_xy_iee_interval(a.dfilter[2], a.dfilter[0], a.dfilter[3], b.z, b.x, c.z, c.x)) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1;
    int l1x_len, l1y_len, l1z_len, d1_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    return orient_xy_iee_exact(l1z, l1z_len, l1x, l1x_len, d1, d1_len, p2, p3);
}

inline int
orient_xy_tee_filtered(const double* l1, const double d1, const double* p2, const double* p3, double max_var) {
    double t1x = p2[1] - p3[1];
    double t1y = p3[0] - p2[0];
    double e2 = l1[0] * t1x;
    double e3 = l1[1] * t1y;
    double e = e2 + e3;
    double pr1 = p2[0] * p3[1];
    double pr2 = p2[1] * p3[0];
    double pr = pr1 - pr2;
    double dpr = d1 * pr;
    double det = dpr + e;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(p2[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p2[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(p3[1])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(t1x)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= 9.061883188277186e-13;
    if ((_tmp_fabs = fabs(t1y)) > max_var) max_var = _tmp_fabs;
    if (det > epsilon) return IPSign::POSITIVE;
    if (-det > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient_xy_tee(const ImplicitPointTPI& a, const ExplicitPoint3D& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_tee_filtered(a.ssfilter.data(), a.ssfilter[3], b.ptr(), c.ptr(), max_var)) != 0.0)
            return ret;
    }

    if (a.get_interval_lambda()) {
        if ((ret = orient_xy_iee_interval(a.dfilter[0], a.dfilter[1], a.dfilter[3], b.x, b.y, c.x, c.y)) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1;
    int l1x_len, l1y_len, l1z_len, d1_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    return orient_xy_iee_exact(l1x, l1x_len, l1y, l1y_len, d1, d1_len, b.ptr(), c.ptr());
}

inline int orient_yz_tee(const ImplicitPointTPI& a, const ExplicitPoint3D& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_tee_filtered(a.ssfilter.data() + 1, a.ssfilter[3], b.ptr() + 1, c.ptr() + 1, max_var)) !=
            0.0)
            return ret;
    }
    if (a.get_interval_lambda()) {
        if ((ret = orient_xy_iee_interval(a.dfilter[1], a.dfilter[2], a.dfilter[3], b.y, b.z, c.y, c.z)) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1;
    int l1x_len, l1y_len, l1z_len, d1_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    return orient_xy_iee_exact(l1y, l1y_len, l1z, l1z_len, d1, d1_len, b.ptr() + 1, c.ptr() + 1);
}

inline int orient_zx_tee(const ImplicitPointTPI& a, const ExplicitPoint3D& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    const double p2[2]{b.z, b.x};
    const double p3[2]{c.z, c.x};
    int ret;
    if (a.get_filtered_lambda(max_var)) {
        const double l1[2]{a.ssfilter[2], a.ssfilter[0]};
        if ((ret = orient_xy_tee_filtered(l1, a.ssfilter[3], p2, p3, max_var)) != 0.0) return ret;
    }
    if (a.get_interval_lambda()) {
        if ((ret = orient_xy_iee_interval(a.dfilter[2], a.dfilter[0], a.dfilter[3], b.z, b.x, c.z, c.x)) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1;
    int l1x_len, l1y_len, l1z_len, d1_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    return orient_xy_iee_exact(l1z, l1z_len, l1x, l1x_len, d1, d1_len, p2, p3);
}

inline int orient_xy_lle_filtered(
    const double* l1, const double d1, const double* l2, const double d2, const double* op3, double max_var
) {
    double a = d1 * l2[0];
    double b = d2 * l1[0];
    double c = d1 * op3[1];
    double e = d1 * l2[1];
    double f = d2 * l1[1];
    double g = d1 * op3[0];
    double ab = a - b;
    double cd = c - l1[1];
    double ef = e - f;
    double gh = g - l1[0];
    double abcd = ab * cd;
    double efgh = ef * gh;
    double L = abcd - efgh;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(op3[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(op3[1])) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 1.699690735379461e-11;

    if (L > epsilon) return IPSign::POSITIVE;
    if (-L > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient_xy_iie_interval(
    const IntervalNumber& l1x, const IntervalNumber& l1y, const IntervalNumber& d1, const IntervalNumber& l2x,
    const IntervalNumber& l2y, const IntervalNumber& d2, const IntervalNumber& op3x, const IntervalNumber& op3y
) {
    IntervalNumber a(d1 * l2x);
    IntervalNumber b(d2 * l1x);
    IntervalNumber c(d1 * op3y);
    IntervalNumber e(d1 * l2y);
    IntervalNumber f(d2 * l1y);
    IntervalNumber g(d1 * op3x);
    IntervalNumber ab(a - b);
    IntervalNumber cd(c - l1y);
    IntervalNumber ef(e - f);
    IntervalNumber gh(g - l1x);
    IntervalNumber abcd(ab * cd);
    IntervalNumber efgh(ef * gh);
    IntervalNumber L(abcd - efgh);
    return L.sign();
}

inline int orient_xy_iie_exact(
    const std::vector<double>& l1x, const int l1x_len, const std::vector<double>& l1y, const int l1y_len,
    const std::vector<double>& d1, const int d1_len, const std::vector<double>& l2x, const int l2x_len,
    const std::vector<double>& l2y, const int l2y_len, const std::vector<double>& d2, const int d2_len,
    const double* op3
) {
    std::vector<double> a(static_cast<uint32_t>((d1_len * l2x_len) << 1));
    int a_len = product_expansion_zeroelim(d1_len, d1.data(), l2x_len, l2x.data(), a.data());
    std::vector<double> b(static_cast<uint32_t>((d2_len * l1x_len) << 1));
    int b_len = product_expansion_zeroelim(d2_len, d2.data(), l1x_len, l1x.data(), b.data());
    std::vector<double> c(static_cast<uint32_t>(d1_len << 1));
    int c_len = scale_expansion_zeroelim(d1_len, d1.data(), op3[1], c.data());
    std::vector<double> e(static_cast<uint32_t>((d1_len * l2y_len) << 1));
    int e_len = product_expansion_zeroelim(d1_len, d1.data(), l2y_len, l2y.data(), e.data());
    std::vector<double> f(static_cast<uint32_t>((d2_len * l1y_len) << 1));
    int f_len = product_expansion_zeroelim(d2_len, d2.data(), l1y_len, l1y.data(), f.data());
    std::vector<double> g(static_cast<uint32_t>(d1_len << 1));
    int g_len = scale_expansion_zeroelim(d1_len, d1.data(), op3[0], g.data());
    std::vector<double> ab(static_cast<uint32_t>(a_len + b_len));
    int ab_len = fast_expansion_diff_zeroelim(a_len, a.data(), b_len, b.data(), ab.data());
    std::vector<double> cd(static_cast<uint32_t>(c_len + l1y_len));
    int cd_len = fast_expansion_diff_zeroelim(c_len, c.data(), l1y_len, l1y.data(), cd.data());
    std::vector<double> ef(static_cast<uint32_t>(e_len + f_len));
    int ef_len = fast_expansion_diff_zeroelim(e_len, e.data(), f_len, f.data(), ef.data());
    std::vector<double> gh(static_cast<uint32_t>(g_len + l1x_len));
    int gh_len = fast_expansion_diff_zeroelim(g_len, g.data(), l1x_len, l1x.data(), gh.data());
    std::vector<double> abcd(static_cast<uint32_t>((ab_len * cd_len) << 1));
    int abcd_len = product_expansion_zeroelim(ab_len, ab.data(), cd_len, cd.data(), abcd.data());
    std::vector<double> efgh(static_cast<uint32_t>((ef_len * gh_len) << 1));
    int efgh_len = product_expansion_zeroelim(ef_len, ef.data(), gh_len, gh.data(), efgh.data());
    std::vector<double> L(static_cast<uint32_t>(abcd_len + efgh_len));
    int L_len = fast_expansion_diff_zeroelim(abcd_len, abcd.data(), efgh_len, efgh.data(), L.data());
    if (L.data()[L_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (L.data()[L_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient_xy_lle(const ImplicitPointLPI& a, const ImplicitPointLPI& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_lle_filtered(
                 a.ssfilter.data(), a.ssfilter[3], b.ssfilter.data(), b.ssfilter[3], c.ptr(), max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda()) {
        if ((ret = orient_xy_iie_interval(
                 a.dfilter[0], a.dfilter[1], a.dfilter[3], b.dfilter[0], b.dfilter[1], b.dfilter[3], c.x, c.y
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    return orient_xy_iie_exact(l1x, l1x_len, l1y, l1y_len, d1, d1_len, l2x, l2x_len, l2y, l2y_len, d2, d2_len, c.ptr());
}

inline int orient_yz_lle(const ImplicitPointLPI& a, const ImplicitPointLPI& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_lle_filtered(
                 a.ssfilter.data() + 1, a.ssfilter[3], b.ssfilter.data() + 1, b.ssfilter[3], c.ptr() + 1, max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda()) {
        if ((ret = orient_xy_iie_interval(
                 a.dfilter[1], a.dfilter[2], a.dfilter[3], b.dfilter[1], b.dfilter[2], b.dfilter[3], c.y, c.z
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    return orient_xy_iie_exact(
        l1y, l1y_len, l1z, l1z_len, d1, d1_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len, c.ptr() + 1
    );
}

inline int orient_zx_lle(const ImplicitPointLPI& a, const ImplicitPointLPI& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    const double p3[2]{c.z, c.x};
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var)) {
        const double l1[2]{a.ssfilter[2], a.ssfilter[0]};
        const double l2[2]{b.ssfilter[2], b.ssfilter[0]};
        if ((ret = orient_xy_lle_filtered(l1, a.ssfilter[3], l2, b.ssfilter[3], p3, max_var)) != 0.0) return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda()) {
        if ((ret = orient_xy_iie_interval(
                 a.dfilter[2], a.dfilter[0], a.dfilter[3], b.dfilter[2], b.dfilter[0], b.dfilter[3], c.z, c.x
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    return orient_xy_iie_exact(l1z, l1z_len, l1x, l1x_len, d1, d1_len, l2z, l2z_len, l2x, l2x_len, d2, d2_len, p3);
}

inline int orient_xy_lte_filtered(
    const double* l1, const double d1, const double* l2, const double d2, const double* op3, double max_var
) {
    double a = d1 * l2[0];
    double b = d2 * l1[0];
    double c = d1 * op3[1];
    double e = d1 * l2[1];
    double f = d2 * l1[1];
    double g = d1 * op3[0];
    double ab = a - b;
    double cd = c - l1[1];
    double ef = e - f;
    double gh = g - l1[0];
    double abcd = ab * cd;
    double efgh = ef * gh;
    double L = abcd - efgh;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(op3[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(op3[1])) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 2.184958117212875e-10;

    if (L > epsilon) return IPSign::POSITIVE;
    if (-L > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient_xy_lte(const ImplicitPointLPI& a, const ImplicitPointTPI& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_lte_filtered(
                 a.ssfilter.data(), a.ssfilter[3], b.ssfilter.data(), b.ssfilter[3], c.ptr(), max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda()) {
        if ((ret = orient_xy_iie_interval(
                 a.dfilter[0], a.dfilter[1], a.dfilter[3], b.dfilter[0], b.dfilter[1], b.dfilter[3], c.x, c.y
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    return orient_xy_iie_exact(l1x, l1x_len, l1y, l1y_len, d1, d1_len, l2x, l2x_len, l2y, l2y_len, d2, d2_len, c.ptr());
}

inline int orient_yz_lte(const ImplicitPointLPI& a, const ImplicitPointTPI& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_lte_filtered(
                 a.ssfilter.data() + 1, a.ssfilter[3], b.ssfilter.data() + 1, b.ssfilter[3], c.ptr() + 1, max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda()) {
        if ((ret = orient_xy_iie_interval(
                 a.dfilter[1], a.dfilter[2], a.dfilter[3], b.dfilter[1], b.dfilter[2], b.dfilter[3], c.y, c.z
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    return orient_xy_iie_exact(
        l1y, l1y_len, l1z, l1z_len, d1, d1_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len, c.ptr() + 1
    );
}

inline int orient_zx_lte(const ImplicitPointLPI& a, const ImplicitPointTPI& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    const double p3[2]{c.z, c.x};
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var)) {
        const double l1[2]{a.ssfilter[2], a.ssfilter[0]};
        const double l2[2]{b.ssfilter[2], b.ssfilter[0]};
        if ((ret = orient_xy_lte_filtered(l1, a.ssfilter[3], l2, b.ssfilter[3], p3, max_var)) != 0.0) return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda()) {
        if ((ret = orient_xy_iie_interval(
                 a.dfilter[2], a.dfilter[0], a.dfilter[3], b.dfilter[2], b.dfilter[0], b.dfilter[3], c.z, c.x
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    return orient_xy_iie_exact(l1z, l1z_len, l1x, l1x_len, d1, d1_len, l2z, l2z_len, l2x, l2x_len, d2, d2_len, p3);
}

inline int orient_xy_tte_filtered(
    const double* l1, const double d1, const double* l2, const double d2, const double* op3, double max_var
) {
    double a = d1 * l2[0];
    double b = d2 * l1[0];
    double c = d1 * op3[1];
    double e = d1 * l2[1];
    double f = d2 * l1[1];
    double g = d1 * op3[0];
    double ab = a - b;
    double cd = c - l1[1];
    double ef = e - f;
    double gh = g - l1[0];
    double abcd = ab * cd;
    double efgh = ef * gh;
    double L = abcd - efgh;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(op3[0])) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(op3[1])) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 3.307187945722514e-08;

    if (L > epsilon) return IPSign::POSITIVE;
    if (-L > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient_xy_tte(const ImplicitPointTPI& a, const ImplicitPointTPI& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_tte_filtered(
                 a.ssfilter.data(), a.ssfilter[3], b.ssfilter.data(), b.ssfilter[3], c.ptr(), max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda()) {
        if ((ret = orient_xy_iie_interval(
                 a.dfilter[0], a.dfilter[1], a.dfilter[3], b.dfilter[0], b.dfilter[1], b.dfilter[3], c.x, c.y
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    return orient_xy_iie_exact(l1x, l1x_len, l1y, l1y_len, d1, d1_len, l2x, l2x_len, l2y, l2y_len, d2, d2_len, c.ptr());
}

inline int orient_yz_tte(const ImplicitPointTPI& a, const ImplicitPointTPI& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_tte_filtered(
                 a.ssfilter.data() + 1, a.ssfilter[3], b.ssfilter.data() + 1, b.ssfilter[3], c.ptr() + 1, max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda()) {
        if ((ret = orient_xy_iie_interval(
                 a.dfilter[1], a.dfilter[2], a.dfilter[3], b.dfilter[1], b.dfilter[2], b.dfilter[3], c.y, c.z
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    return orient_xy_iie_exact(
        l1y, l1y_len, l1z, l1z_len, d1, d1_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len, c.ptr() + 1
    );
}

inline int orient_zx_tte(const ImplicitPointTPI& a, const ImplicitPointTPI& b, const ExplicitPoint3D& c) {
    double max_var = 0.0;
    int ret;
    const double p3[2]{c.z, c.x};
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var)) {
        const double l1[2]{a.ssfilter[2], a.ssfilter[0]};
        const double l2[2]{b.ssfilter[2], b.ssfilter[0]};
        if ((ret = orient_xy_tte_filtered(l1, a.ssfilter[3], l2, b.ssfilter[3], p3, max_var)) != 0.0) return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda()) {
        if ((ret = orient_xy_iie_interval(
                 a.dfilter[2], a.dfilter[0], a.dfilter[3], b.dfilter[2], b.dfilter[0], b.dfilter[3], c.z, c.x
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    return orient_xy_iie_exact(l1z, l1z_len, l1x, l1x_len, d1, d1_len, l2z, l2z_len, l2x, l2x_len, d2, d2_len, p3);
}

inline int orient_xy_lll_filtered(
    const double* l1, const double d1, const double* l2, const double d2, const double* l3, const double d3,
    double max_var
) {
    double a = d1 * l2[0];
    double b = d2 * l1[0];
    double c = d1 * l3[1];
    double d = d3 * l1[1];
    double e = d1 * l2[1];
    double f = d2 * l1[1];
    double g = d1 * l3[0];
    double h = d3 * l1[0];
    double ab = a - b;
    double cd = c - d;
    double ef = e - f;
    double gh = g - h;
    double abcd = ab * cd;
    double efgh = ef * gh;
    double L = abcd - efgh;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 1.75634284893534e-10;
    if (L > epsilon) return IPSign::POSITIVE;
    if (-L > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient_xy_iii_interval(
    const IntervalNumber& l1x, const IntervalNumber& l1y, const IntervalNumber& d1, const IntervalNumber& l2x,
    const IntervalNumber& l2y, const IntervalNumber& d2, const IntervalNumber& l3x, const IntervalNumber& l3y,
    const IntervalNumber& d3
) {
    IntervalNumber a(d1 * l2x);
    IntervalNumber b(d2 * l1x);
    IntervalNumber c(d1 * l3y);
    IntervalNumber d(d3 * l1y);
    IntervalNumber e(d1 * l2y);
    IntervalNumber f(d2 * l1y);
    IntervalNumber g(d1 * l3x);
    IntervalNumber h(d3 * l1x);
    IntervalNumber ab(a - b);
    IntervalNumber cd(c - d);
    IntervalNumber ef(e - f);
    IntervalNumber gh(g - h);
    IntervalNumber abcd(ab * cd);
    IntervalNumber efgh(ef * gh);
    IntervalNumber L(abcd - efgh);
    return L.sign();
}
inline int orient_xy_iii_exact(
    const std::vector<double>& l1x, const int l1x_len, const std::vector<double>& l1y, const int l1y_len,
    const std::vector<double>& d1, const int d1_len, const std::vector<double>& l2x, const int l2x_len,
    const std::vector<double>& l2y, const int l2y_len, const std::vector<double>& d2, const int d2_len,
    const std::vector<double>& l3x, const int l3x_len, const std::vector<double>& l3y, const int l3y_len,
    const std::vector<double>& d3, const int d3_len
) {
    std::vector<double> a(static_cast<uint32_t>((d1_len * l2x_len) << 1));
    int a_len = product_expansion_zeroelim(d1_len, d1.data(), l2x_len, l2x.data(), a.data());
    std::vector<double> b(static_cast<uint32_t>((d2_len * l1x_len) << 1));
    int b_len = product_expansion_zeroelim(d2_len, d2.data(), l1x_len, l1x.data(), b.data());
    std::vector<double> c(static_cast<uint32_t>((d1_len * l3y_len) << 1));
    int c_len = product_expansion_zeroelim(d1_len, d1.data(), l3y_len, l3y.data(), c.data());
    std::vector<double> d(static_cast<uint32_t>((d3_len * l1y_len) << 1));
    int d_len = product_expansion_zeroelim(d3_len, d3.data(), l1y_len, l1y.data(), d.data());
    std::vector<double> e(static_cast<uint32_t>((d1_len * l2y_len) << 1));
    int e_len = product_expansion_zeroelim(d1_len, d1.data(), l2y_len, l2y.data(), e.data());
    std::vector<double> f(static_cast<uint32_t>((d2_len * l1y_len) << 1));
    int f_len = product_expansion_zeroelim(d2_len, d2.data(), l1y_len, l1y.data(), f.data());
    std::vector<double> g(static_cast<uint32_t>((d1_len * l3x_len) << 1));
    int g_len = product_expansion_zeroelim(d1_len, d1.data(), l3x_len, l3x.data(), g.data());
    std::vector<double> h(static_cast<uint32_t>((d3_len * l1x_len) << 1));
    int h_len = product_expansion_zeroelim(d3_len, d3.data(), l1x_len, l1x.data(), h.data());
    std::vector<double> ab(static_cast<uint32_t>(a_len + b_len));
    int ab_len = fast_expansion_diff_zeroelim(a_len, a.data(), b_len, b.data(), ab.data());
    std::vector<double> cd(static_cast<uint32_t>(c_len + d_len));
    int cd_len = fast_expansion_diff_zeroelim(c_len, c.data(), d_len, d.data(), cd.data());
    std::vector<double> ef(static_cast<uint32_t>(e_len + f_len));
    int ef_len = fast_expansion_diff_zeroelim(e_len, e.data(), f_len, f.data(), ef.data());
    std::vector<double> gh(static_cast<uint32_t>(g_len + h_len));
    int gh_len = fast_expansion_diff_zeroelim(g_len, g.data(), h_len, h.data(), gh.data());
    std::vector<double> abcd(static_cast<uint32_t>((ab_len * cd_len) << 1));
    int abcd_len = product_expansion_zeroelim(ab_len, ab.data(), cd_len, cd.data(), abcd.data());
    std::vector<double> efgh(static_cast<uint32_t>((ef_len * gh_len) << 1));
    int efgh_len = product_expansion_zeroelim(ef_len, ef.data(), gh_len, gh.data(), efgh.data());
    std::vector<double> L(static_cast<uint32_t>(abcd_len + efgh_len));
    int L_len = fast_expansion_diff_zeroelim(abcd_len, abcd.data(), efgh_len, efgh.data(), L.data());
    if (L.data()[L_len - 1] > 0) {
        return IPSign::POSITIVE;
    } else if (L.data()[L_len - 1] < 0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int orient_xy_lll(const ImplicitPointLPI& a, const ImplicitPointLPI& b, const ImplicitPointLPI& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var) && c.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_lll_filtered(
                 a.ssfilter.data(), a.ssfilter[3], b.ssfilter.data(), b.ssfilter[3], c.ssfilter.data(), c.ssfilter[3],
                 max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[0], a.dfilter[1], a.dfilter[3], b.dfilter[0], b.dfilter[1], b.dfilter[3], c.dfilter[0],
                 c.dfilter[1], c.dfilter[3]
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1x, l1x_len, l1y, l1y_len, d1, d1_len, l2x, l2x_len, l2y, l2y_len, d2, d2_len, l3x, l3x_len, l3y, l3y_len, d3,
        d3_len
    );
}

inline int orient_yz_lll(const ImplicitPointLPI& a, const ImplicitPointLPI& b, const ImplicitPointLPI& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var) && c.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_lll_filtered(
                 a.ssfilter.data() + 1, a.ssfilter[3], b.ssfilter.data() + 1, b.ssfilter[3], c.ssfilter.data() + 1,
                 c.ssfilter[3], max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[1], a.dfilter[2], a.dfilter[3], b.dfilter[1], b.dfilter[2], b.dfilter[3], c.dfilter[1],
                 c.dfilter[2], c.dfilter[3]
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1y, l1y_len, l1z, l1z_len, d1, d1_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len, l3y, l3y_len, l3z, l3z_len, d3,
        d3_len
    );
}

inline int orient_zx_lll(const ImplicitPointLPI& a, const ImplicitPointLPI& b, const ImplicitPointLPI& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var) && c.get_filtered_lambda(max_var)) {
        const double l1[2]{a.ssfilter[2], a.ssfilter[0]};
        const double l2[2]{b.ssfilter[2], b.ssfilter[0]};
        const double l3[2]{c.ssfilter[2], c.ssfilter[0]};
        if ((ret = orient_xy_lll_filtered(l1, a.ssfilter[3], l2, b.ssfilter[3], l3, c.ssfilter[3], max_var)) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[2], a.dfilter[0], a.dfilter[3], b.dfilter[2], b.dfilter[0], b.dfilter[3], c.dfilter[2],
                 c.dfilter[0], c.dfilter[3]
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1z, l1z_len, l1x, l1x_len, d1, d1_len, l2z, l2z_len, l2x, l2x_len, d2, d2_len, l3z, l3z_len, l3x, l3x_len, d3,
        d3_len
    );
}

inline int orient_xy_llt_filtered(
    const double* l1, const double d1, const double* l2, const double d2, const double* l3, const double d3,
    double max_var
) {
    double a = d1 * l2[0];
    double b = d2 * l1[0];
    double c = d1 * l3[1];
    double d = d3 * l1[1];
    double e = d1 * l2[1];
    double f = d2 * l1[1];
    double g = d1 * l3[0];
    double h = d3 * l1[0];
    double ab = a - b;
    double cd = c - d;
    double ef = e - f;
    double gh = g - h;
    double abcd = ab * cd;
    double efgh = ef * gh;
    double L = abcd - efgh;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= 2.144556754402072e-09;
    if (L > epsilon) return IPSign::POSITIVE;
    if (-L > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient_xy_llt(const ImplicitPointLPI& a, const ImplicitPointLPI& b, const ImplicitPointTPI& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var) && c.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_llt_filtered(
                 a.ssfilter.data(), a.ssfilter[3], b.ssfilter.data(), b.ssfilter[3], c.ssfilter.data(), c.ssfilter[3],
                 max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[0], a.dfilter[1], a.dfilter[3], b.dfilter[0], b.dfilter[1], b.dfilter[3], c.dfilter[0],
                 c.dfilter[1], c.dfilter[3]
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1x, l1x_len, l1y, l1y_len, d1, d1_len, l2x, l2x_len, l2y, l2y_len, d2, d2_len, l3x, l3x_len, l3y, l3y_len, d3,
        d3_len
    );
}

inline int orient_yz_llt(const ImplicitPointLPI& a, const ImplicitPointLPI& b, const ImplicitPointTPI& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var) && c.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_llt_filtered(
                 a.ssfilter.data() + 1, a.ssfilter[3], b.ssfilter.data() + 1, b.ssfilter[3], c.ssfilter.data() + 1,
                 c.ssfilter[3], max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[1], a.dfilter[2], a.dfilter[3], b.dfilter[1], b.dfilter[2], b.dfilter[3], c.dfilter[1],
                 c.dfilter[2], c.dfilter[3]
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1y, l1y_len, l1z, l1z_len, d1, d1_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len, l3y, l3y_len, l3z, l3z_len, d3,
        d3_len
    );
}

inline int orient_zx_llt(const ImplicitPointLPI& a, const ImplicitPointLPI& b, const ImplicitPointTPI& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var) && c.get_filtered_lambda(max_var)) {
        const double l1[2]{a.ssfilter[2], a.ssfilter[0]};
        const double l2[2]{b.ssfilter[2], b.ssfilter[0]};
        const double l3[2]{c.ssfilter[2], c.ssfilter[0]};
        if ((ret = orient_xy_llt_filtered(l1, a.ssfilter[3], l2, b.ssfilter[3], l3, c.ssfilter[3], max_var)) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[2], a.dfilter[0], a.dfilter[3], b.dfilter[2], b.dfilter[0], b.dfilter[3], c.dfilter[2],
                 c.dfilter[0], c.dfilter[3]
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1z, l1z_len, l1x, l1x_len, d1, d1_len, l2z, l2z_len, l2x, l2x_len, d2, d2_len, l3z, l3z_len, l3x, l3x_len, d3,
        d3_len
    );
}

inline int orient_xy_ltt_filtered(
    const double* l1, const double d1, const double* l2, const double d2, const double* l3, const double d3,
    double max_var
) {
    double a = d1 * l2[0];
    double b = d2 * l1[0];
    double c = d1 * l3[1];
    double d = d3 * l1[1];
    double e = d1 * l2[1];
    double f = d2 * l1[1];
    double g = d1 * l3[0];
    double h = d3 * l1[0];
    double ab = a - b;
    double cd = c - d;
    double ef = e - f;
    double gh = g - h;
    double abcd = ab * cd;
    double efgh = ef * gh;
    double L = abcd - efgh;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 2.535681042914479e-08;
    if (L > epsilon) return IPSign::POSITIVE;
    if (-L > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int orient_xy_ltt(const ImplicitPointLPI& a, const ImplicitPointTPI& b, const ImplicitPointTPI& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var) && c.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_ltt_filtered(
                 a.ssfilter.data(), a.ssfilter[3], b.ssfilter.data(), b.ssfilter[3], c.ssfilter.data(), c.ssfilter[3],
                 max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[0], a.dfilter[1], a.dfilter[3], b.dfilter[0], b.dfilter[1], b.dfilter[3], c.dfilter[0],
                 c.dfilter[1], c.dfilter[3]
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1x, l1x_len, l1y, l1y_len, d1, d1_len, l2x, l2x_len, l2y, l2y_len, d2, d2_len, l3x, l3x_len, l3y, l3y_len, d3,
        d3_len
    );
}

inline int orient_yz_ltt(const ImplicitPointLPI& a, const ImplicitPointTPI& b, const ImplicitPointTPI& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var) && c.get_filtered_lambda(max_var)) {
        if ((ret = orient_xy_ltt_filtered(
                 a.ssfilter.data() + 1, a.ssfilter[3], b.ssfilter.data() + 1, b.ssfilter[3], c.ssfilter.data() + 1,
                 c.ssfilter[3], max_var
             )) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[1], a.dfilter[2], a.dfilter[3], b.dfilter[1], b.dfilter[2], b.dfilter[3], c.dfilter[1],
                 c.dfilter[2], c.dfilter[3]
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1y, l1y_len, l1z, l1z_len, d1, d1_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len, l3y, l3y_len, l3z, l3z_len, d3,
        d3_len
    );
}

inline int orient_zx_ltt(const ImplicitPointLPI& a, const ImplicitPointTPI& b, const ImplicitPointTPI& c) {
    double max_var = 0.0;
    int ret;
    if (a.get_filtered_lambda(max_var) && b.get_filtered_lambda(max_var) && c.get_filtered_lambda(max_var)) {
        const double l1[2]{a.ssfilter[2], a.ssfilter[0]};
        const double l2[2]{b.ssfilter[2], b.ssfilter[0]};
        const double l3[2]{c.ssfilter[2], c.ssfilter[0]};
        if ((ret = orient_xy_ltt_filtered(l1, a.ssfilter[3], l2, b.ssfilter[3], l3, c.ssfilter[3], max_var)) != 0.0)
            return ret;
    }
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[2], a.dfilter[0], a.dfilter[3], b.dfilter[2], b.dfilter[0], b.dfilter[3], c.dfilter[2],
                 c.dfilter[0], c.dfilter[3]
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1z, l1z_len, l1x, l1x_len, d1, d1_len, l2z, l2z_len, l2x, l2x_len, d2, d2_len, l3z, l3z_len, l3x, l3x_len, d3,
        d3_len
    );
}

inline int orient_xy_ttt(const ImplicitPointTPI& a, const ImplicitPointTPI& b, const ImplicitPointTPI& c) {
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        int ret;
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[0], a.dfilter[1], a.dfilter[3], b.dfilter[0], b.dfilter[1], b.dfilter[3], c.dfilter[0],
                 c.dfilter[1], c.dfilter[3]
             )) != 0.0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1x, l1x_len, l1y, l1y_len, d1, d1_len, l2x, l2x_len, l2y, l2y_len, d2, d2_len, l3x, l3x_len, l3y, l3y_len, d3,
        d3_len
    );
}

inline int orient_yz_ttt(const ImplicitPointTPI& a, const ImplicitPointTPI& b, const ImplicitPointTPI& c) {
    if (a.get_interval_lambda() && b.get_interval_lambda() && c.get_interval_lambda()) {
        int ret;
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[1], a.dfilter[2], a.dfilter[3], b.dfilter[1], b.dfilter[2], b.dfilter[3], c.dfilter[1],
                 c.dfilter[2], c.dfilter[3]
             )) != 0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1y, l1y_len, l1z, l1z_len, d1, d1_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len, l3y, l3y_len, l3z, l3z_len, d3,
        d3_len
    );
}

inline int orient_zx_ttt(const ImplicitPointTPI& a, const ImplicitPointTPI& b, const ImplicitPointTPI& c) {
    if (a.get_interval_lambda() && !b.get_interval_lambda() && c.get_interval_lambda()) {
        int ret;
        if ((ret = orient_xy_iii_interval(
                 a.dfilter[2], a.dfilter[0], a.dfilter[3], b.dfilter[2], b.dfilter[0], b.dfilter[3], c.dfilter[2],
                 c.dfilter[0], c.dfilter[3]
             )) != 0)
            return ret;
    }

    std::vector<double> l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
    int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
    a.get_exact_lambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
    b.get_exact_lambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
    c.get_exact_lambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
    return orient_xy_iii_exact(
        l1z, l1z_len, l1x, l1x_len, d1, d1_len, l2z, l2z_len, l2x, l2x_len, d2, d2_len, l3z, l3z_len, l3x, l3x_len, d3,
        d3_len
    );
}

int GenericPoint3D::orient_xy(const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c) {
    const int val = a.get_type() * 9 + b.get_type() * 3 + c.get_type();
    switch (val) {
    case 0: {
        const double ret = ::orient2d(a.to_explicit().ptr(), b.to_explicit().ptr(), c.to_explicit().ptr());
        return (ret > 0) - (ret < 0);
    }
    case 1: // EEL
        return orient_xy_lee(c.to_lpi(), a.to_explicit(), b.to_explicit());
    case 2: // EET
        return orient_xy_tee(c.to_tpi(), a.to_explicit(), b.to_explicit());
    case 3: // ELE
        return orient_xy_lee(b.to_lpi(), c.to_explicit(), a.to_explicit());
    case 4: // ELL
        return orient_xy_lle(b.to_lpi(), c.to_lpi(), a.to_explicit());
    case 5: // ELT
        return orient_xy_lte(b.to_lpi(), c.to_tpi(), a.to_explicit());
    case 6: // ETE
        return orient_xy_tee(b.to_tpi(), c.to_explicit(), a.to_explicit());
    case 7: // ETL
        return -orient_xy_lte(c.to_lpi(), b.to_tpi(), a.to_explicit());
    case 8: // ETT
        return orient_xy_tte(b.to_tpi(), c.to_tpi(), a.to_explicit());
    case 9: // LEE
        return orient_xy_lee(a.to_lpi(), b.to_explicit(), c.to_explicit());
    case 10: // LEL
        return orient_xy_lle(c.to_lpi(), a.to_lpi(), b.to_explicit());
    case 11: // LET
        return -orient_xy_lte(a.to_lpi(), c.to_tpi(), b.to_explicit());
    case 12: // LLE
        return orient_xy_lle(a.to_lpi(), b.to_lpi(), c.to_explicit());
    case 13: // LLL
        return orient_xy_lll(a.to_lpi(), b.to_lpi(), c.to_lpi());
    case 14: // LLT
        return orient_xy_llt(a.to_lpi(), b.to_lpi(), c.to_tpi());
    case 15: // LTE
        return orient_xy_lte(a.to_lpi(), b.to_tpi(), c.to_explicit());
    case 16: // LTL
        return orient_xy_llt(c.to_lpi(), a.to_lpi(), b.to_tpi());
    case 17: // LTT
        return orient_xy_ltt(a.to_lpi(), b.to_tpi(), c.to_tpi());
    case 18: // TEE
        return orient_xy_tee(a.to_tpi(), b.to_explicit(), c.to_explicit());
    case 19: // TEL
        return orient_xy_lte(c.to_lpi(), a.to_tpi(), b.to_explicit());
    case 20: // TET
        return orient_xy_tte(c.to_tpi(), a.to_tpi(), b.to_explicit());
    case 21: // TLE
        return -orient_xy_lte(b.to_lpi(), a.to_tpi(), c.to_explicit());
    case 22: // TLL
        return orient_xy_llt(b.to_lpi(), c.to_lpi(), a.to_tpi());
    case 23: // TLT
        return orient_xy_ltt(b.to_lpi(), c.to_tpi(), a.to_tpi());
    case 24: // TTE
        return orient_xy_tte(a.to_tpi(), b.to_tpi(), c.to_explicit());
    case 25: // TLT
        return orient_xy_ltt(c.to_lpi(), a.to_tpi(), b.to_tpi());
    case 26: // TTT
        return orient_xy_ttt(a.to_tpi(), b.to_tpi(), c.to_tpi());
    default:
        return IPSign::UNDEFINED;
    }
}

int GenericPoint3D::orient_yz(const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c) {
    const int val = a.get_type() * 9 + b.get_type() * 3 + c.get_type();
    switch (val) {
    case 0: {
        const double ret = ::orient2d(a.to_explicit().ptr() + 1, b.to_explicit().ptr() + 1, c.to_explicit().ptr() + 1);
        return (ret > 0) - (ret < 0);
    }
    case 1: // EEL
        return orient_yz_lee(c.to_lpi(), a.to_explicit(), b.to_explicit());
    case 2: // EET
        return orient_yz_tee(c.to_tpi(), a.to_explicit(), b.to_explicit());
    case 3: // ELE
        return orient_yz_lee(b.to_lpi(), c.to_explicit(), a.to_explicit());
    case 4: // ELL
        return orient_yz_lle(b.to_lpi(), c.to_lpi(), a.to_explicit());
    case 5: // ELT
        return orient_yz_lte(b.to_lpi(), c.to_tpi(), a.to_explicit());
    case 6: // ETE
        return orient_yz_tee(b.to_tpi(), c.to_explicit(), a.to_explicit());
    case 7: // ETL
        return -orient_yz_lte(c.to_lpi(), b.to_tpi(), a.to_explicit());
    case 8: // ETT
        return orient_yz_tte(b.to_tpi(), c.to_tpi(), a.to_explicit());
    case 9: // LEE
        return orient_yz_lee(a.to_lpi(), b.to_explicit(), c.to_explicit());
    case 10: // LEL
        return orient_yz_lle(c.to_lpi(), a.to_lpi(), b.to_explicit());
    case 11: // LET
        return -orient_yz_lte(a.to_lpi(), c.to_tpi(), b.to_explicit());
    case 12: // LLE
        return orient_yz_lle(a.to_lpi(), b.to_lpi(), c.to_explicit());
    case 13: // LLL
        return orient_yz_lll(a.to_lpi(), b.to_lpi(), c.to_lpi());
    case 14: // LLT
        return orient_yz_llt(a.to_lpi(), b.to_lpi(), c.to_tpi());
    case 15: // LTE
        return orient_yz_lte(a.to_lpi(), b.to_tpi(), c.to_explicit());
    case 16: // LTL
        return orient_yz_llt(c.to_lpi(), a.to_lpi(), b.to_tpi());
    case 17: // LTT
        return orient_yz_ltt(a.to_lpi(), b.to_tpi(), c.to_tpi());
    case 18: // TEE
        return orient_yz_tee(a.to_tpi(), b.to_explicit(), c.to_explicit());
    case 19: // TEL
        return orient_yz_lte(c.to_lpi(), a.to_tpi(), b.to_explicit());
    case 20: // TET
        return orient_yz_tte(c.to_tpi(), a.to_tpi(), b.to_explicit());
    case 21: // TLE
        return -orient_yz_lte(b.to_lpi(), a.to_tpi(), c.to_explicit());
    case 22: // TLL
        return orient_yz_llt(b.to_lpi(), c.to_lpi(), a.to_tpi());
    case 23: // TLT
        return orient_yz_ltt(b.to_lpi(), c.to_tpi(), a.to_tpi());
    case 24: // TTE
        return orient_yz_tte(a.to_tpi(), b.to_tpi(), c.to_explicit());
    case 25: // TLT
        return orient_yz_ltt(c.to_lpi(), a.to_tpi(), b.to_tpi());
    case 26: // TTT
        return orient_yz_ttt(a.to_tpi(), b.to_tpi(), c.to_tpi());
    default:
        return IPSign::UNDEFINED;
    }
}

int GenericPoint3D::orient_zx(const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c) {
    const int val = a.get_type() * 9 + b.get_type() * 3 + c.get_type();
    switch (val) {
    case 0: {
        const ExplicitPoint3D& pa = a.to_explicit();
        const ExplicitPoint3D& pb = b.to_explicit();
        const ExplicitPoint3D& pc = c.to_explicit();
        const double aptr[2] = {pa.z, pa.x};
        const double bptr[2] = {pb.z, pb.x};
        const double cptr[2] = {pc.z, pc.x};
        const double ret = ::orient2d(aptr, bptr, cptr);
        return (ret > 0) - (ret < 0);
    }
    case 1: // EEL
        return orient_zx_lee(c.to_lpi(), a.to_explicit(), b.to_explicit());
    case 2: // EET
        return orient_zx_tee(c.to_tpi(), a.to_explicit(), b.to_explicit());
    case 3: // ELE
        return orient_zx_lee(b.to_lpi(), c.to_explicit(), a.to_explicit());
    case 4: // ELL
        return orient_zx_lle(b.to_lpi(), c.to_lpi(), a.to_explicit());
    case 5: // ELT
        return orient_zx_lte(b.to_lpi(), c.to_tpi(), a.to_explicit());
    case 6: // ETE
        return orient_zx_tee(b.to_tpi(), c.to_explicit(), a.to_explicit());
    case 7: // ETL
        return -orient_zx_lte(c.to_lpi(), b.to_tpi(), a.to_explicit());
    case 8: // ETT
        return orient_zx_tte(b.to_tpi(), c.to_tpi(), a.to_explicit());
    case 9: // LEE
        return orient_zx_lee(a.to_lpi(), b.to_explicit(), c.to_explicit());
    case 10: // LEL
        return orient_zx_lle(c.to_lpi(), a.to_lpi(), b.to_explicit());
    case 11: // LET
        return -orient_zx_lte(a.to_lpi(), c.to_tpi(), b.to_explicit());
    case 12: // LLE
        return orient_zx_lle(a.to_lpi(), b.to_lpi(), c.to_explicit());
    case 13: // LLL
        return orient_zx_lll(a.to_lpi(), b.to_lpi(), c.to_lpi());
    case 14: // LLT
        return orient_zx_llt(a.to_lpi(), b.to_lpi(), c.to_tpi());
    case 15: // LTE
        return orient_zx_lte(a.to_lpi(), b.to_tpi(), c.to_explicit());
    case 16: // LTL
        return orient_zx_llt(c.to_lpi(), a.to_lpi(), b.to_tpi());
    case 17: // LTT
        return orient_zx_ltt(a.to_lpi(), b.to_tpi(), c.to_tpi());
    case 18: // TEE
        return orient_zx_tee(a.to_tpi(), b.to_explicit(), c.to_explicit());
    case 19: // TEL
        return orient_zx_lte(c.to_lpi(), a.to_tpi(), b.to_explicit());
    case 20: // TET
        return orient_zx_tte(c.to_tpi(), a.to_tpi(), b.to_explicit());
    case 21: // TLE
        return -orient_zx_lte(b.to_lpi(), a.to_tpi(), c.to_explicit());
    case 22: // TLL
        return orient_zx_llt(b.to_lpi(), c.to_lpi(), a.to_tpi());
    case 23: // TLT
        return orient_zx_ltt(b.to_lpi(), c.to_tpi(), a.to_tpi());
    case 24: // TTE
        return orient_zx_tte(a.to_tpi(), b.to_tpi(), c.to_explicit());
    case 25: // TLT
        return orient_zx_ltt(c.to_lpi(), a.to_tpi(), b.to_tpi());
    case 26: // TTT
        return orient_zx_ttt(a.to_tpi(), b.to_tpi(), c.to_tpi());
    default:
        return IPSign::UNDEFINED;
    }
}

inline int less_than_le_filtered(const double& l1x, const double& d1, double max_var, const double& bx) {
    double dbx = bx * d1;
    double kx = l1x - dbx;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(bx)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= 1.932297637868842e-14;
    if (kx > epsilon) return IPSign::POSITIVE;
    if (-kx > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int less_than_xe_interval(const IntervalNumber& l1x, const IntervalNumber& d1, const IntervalNumber& bx) {
    IntervalNumber dbx(bx * d1);
    IntervalNumber kx(l1x - dbx);
    if (!kx.signIsReliable()) return IPSign::ZERO;
    return kx.sign();
}

inline int less_than_xe_exact(
    const std::vector<double>& l1x, const int l1x_len, const std::vector<double>& d1, const int d1_len, const double bx
) {
    std::vector<double> dbx(128);
    int dbx_len = scale_expansion_zeroelim(d1_len, d1.data(), bx, dbx.data());
    std::vector<double> kx(128);
    int kx_len = fast_expansion_diff_zeroelim(l1x_len, l1x.data(), dbx_len, dbx.data(), kx.data());

    const double return_value = kx.data()[kx_len - 1];
    if (return_value > 0) return IPSign::POSITIVE;
    if (return_value < 0) return IPSign::NEGATIVE;
    if (return_value == 0.0) return IPSign::ZERO;
    return IPSign::UNDEFINED;
}

inline int less_than_on_x_le(const ImplicitPointLPI& a, const double& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv)) {
        if ((ret = less_than_le_filtered(a.ssfilter[0], a.ssfilter[3], mv, b)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda()) {
        if ((ret = less_than_xe_interval(a.dfilter[0], a.dfilter[3], b)) != 0) {
            return ret;
        }
    }
    std::vector<double> lx, ly, lz, d;
    int lxl, lyl, lzl, dl;
    a.get_exact_lambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
    if (d.data()[dl - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xe_exact(lx, lxl, d, dl, b);
}

inline int less_than_on_y_le(const ImplicitPointLPI& a, const double& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv)) {
        if ((ret = less_than_le_filtered(a.ssfilter[1], a.ssfilter[3], mv, b)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda()) {
        if ((ret = less_than_xe_interval(a.dfilter[1], a.dfilter[3], b)) != 0) {
            return ret;
        }
    }
    std::vector<double> lx, ly, lz, d;
    int lxl, lyl, lzl, dl;
    a.get_exact_lambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
    if (d.data()[dl - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xe_exact(ly, lyl, d, dl, b);
}

inline int less_than_on_z_le(const ImplicitPointLPI& a, const double& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv)) {
        if ((ret = less_than_le_filtered(a.ssfilter[2], a.ssfilter[3], mv, b)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda()) {
        if ((ret = less_than_xe_interval(a.dfilter[2], a.dfilter[3], b)) != 0) {
            return ret;
        }
    }
    std::vector<double> lx, ly, lz, d;
    int lxl, lyl, lzl, dl;
    a.get_exact_lambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
    if (d.data()[dl - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xe_exact(lz, lzl, d, dl, b);
}

inline int less_than_te_filtered(const double& l1x, const double& d1, double max_var, const double& bx) {
    double dbx = bx * d1;
    double kx = l1x - dbx;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(bx)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 3.980270973924514e-13;
    if (kx > epsilon) return IPSign::POSITIVE;
    if (-kx > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int less_than_on_x_te(const ImplicitPointTPI& a, const double& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv)) {
        if ((ret = less_than_te_filtered(a.ssfilter[0], a.ssfilter[3], mv, b)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda()) {
        if ((ret = less_than_xe_interval(a.dfilter[0], a.dfilter[3], b)) != 0) {
            return ret;
        }
    }
    std::vector<double> lx, ly, lz, d;
    int lxl, lyl, lzl, dl;
    a.get_exact_lambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
    if (d.data()[dl - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xe_exact(lx, lxl, d, dl, b);
}

inline int less_than_on_y_te(const ImplicitPointTPI& a, const double& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv)) {
        if ((ret = less_than_te_filtered(a.ssfilter[1], a.ssfilter[3], mv, b)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda()) {
        if ((ret = less_than_xe_interval(a.dfilter[1], a.dfilter[3], b)) != 0) {
            return ret;
        }
    }
    std::vector<double> lx, ly, lz, d;
    int lxl, lyl, lzl, dl;
    a.get_exact_lambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
    if (d.data()[dl - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xe_exact(ly, lyl, d, dl, b);
}

inline int less_than_on_z_te(const ImplicitPointTPI& a, const double& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv)) {
        if ((ret = less_than_te_filtered(a.ssfilter[2], a.ssfilter[3], mv, b)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda()) {
        if ((ret = less_than_xe_interval(a.dfilter[2], a.dfilter[3], b)) != 0) {
            return ret;
        }
    }
    std::vector<double> lx, ly, lz, d;
    int lxl, lyl, lzl, dl;
    a.get_exact_lambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
    if (d.data()[dl - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xe_exact(lz, lzl, d, dl, b);
}

inline int
less_than_ll_filtered(const double& l1x, const double& d1, const double& l2x, const double& d2, double max_var) {
    double k1 = d2 * l1x;
    double k2 = d1 * l2x;
    double kx = k1 - k2;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 2.922887626377607e-13;
    if (kx > epsilon) return IPSign::POSITIVE;
    if (-kx > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int less_than_xx_interval(
    const IntervalNumber& l1x, const IntervalNumber& d1, const IntervalNumber& l2x, const IntervalNumber& d2
) {
    IntervalNumber k1(d2 * l1x);
    IntervalNumber k2(d1 * l2x);
    IntervalNumber kx(k1 - k2);
    if (!kx.signIsReliable()) return IPSign::ZERO;
    return kx.sign();
}

inline int less_than_xx_exact(
    const std::vector<double>& l1x, const int l1x_len, const std::vector<double>& d1, const int d1_len,
    const std::vector<double>& l2x, const int l2x_len, const std::vector<double>& d2, const int d2_len
) {
    std::vector<double> k1(128);
    int k1_len = product_expansion_zeroelim(d2_len, d2.data(), l1x_len, l1x.data(), k1.data());
    std::vector<double> k2(128);
    int k2_len = product_expansion_zeroelim(d1_len, d1.data(), l2x_len, l2x.data(), k2.data());
    std::vector<double> kx(128);
    int kx_len = fast_expansion_diff_zeroelim(k1_len, k1.data(), k2_len, k2.data(), kx.data());
    const double return_value = kx.data()[kx_len - 1];
    if (return_value > 0.0) {
        return IPSign::POSITIVE;
    } else if (return_value < 0.0) {
        return IPSign::NEGATIVE;
    } else {
        return IPSign::ZERO;
    }
}

inline int less_than_on_x_ll(const ImplicitPointLPI& a, const ImplicitPointLPI& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv) && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_ll_filtered(a.ssfilter[0], a.ssfilter[3], b.ssfilter[0], b.ssfilter[3], mv)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda() && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_xx_interval(a.dfilter[0], a.dfilter[3], b.dfilter[0], b.dfilter[3])) != 0) {
            return ret;
        }
    }
    std::vector<double> l1x, l1y, l1z, d1;
    int l1xl, l1yl, l1zl, d1l;
    std::vector<double> l2x, l2y, l2z, d2;
    int l2xl, l2yl, l2zl, d2l;
    a.get_exact_lambda(l1x, l1xl, l1y, l1yl, l1z, l1zl, d1, d1l);
    b.get_exact_lambda(l2x, l2xl, l2y, l2yl, l2z, l2zl, d2, d2l);
    if (d1.data()[d1l - 1] == 0.0 || d2.data()[d2l - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xx_exact(l1x, l1xl, d1, d1l, l2x, l2xl, d2, d2l);
}

inline int less_than_on_y_ll(const ImplicitPointLPI& a, const ImplicitPointLPI& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv) && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_ll_filtered(a.ssfilter[1], a.ssfilter[3], b.ssfilter[1], b.ssfilter[3], mv)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda() && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_xx_interval(a.dfilter[1], a.dfilter[3], b.dfilter[1], b.dfilter[3])) != 0) {
            return ret;
        }
    }
    std::vector<double> l1x, l1y, l1z, d1;
    int l1xl, l1yl, l1zl, d1l;
    std::vector<double> l2x, l2y, l2z, d2;
    int l2xl, l2yl, l2zl, d2l;
    a.get_exact_lambda(l1x, l1xl, l1y, l1yl, l1z, l1zl, d1, d1l);
    b.get_exact_lambda(l2x, l2xl, l2y, l2yl, l2z, l2zl, d2, d2l);
    if (d1.data()[d1l - 1] == 0.0 || d2.data()[d2l - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xx_exact(l1y, l1yl, d1, d1l, l2y, l2yl, d2, d2l);
}

inline int less_than_on_z_ll(const ImplicitPointLPI& a, const ImplicitPointLPI& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv) && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_ll_filtered(a.ssfilter[2], a.ssfilter[3], b.ssfilter[2], b.ssfilter[3], mv)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda() && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_xx_interval(a.dfilter[2], a.dfilter[3], b.dfilter[2], b.dfilter[3])) != 0) {
            return ret;
        }
    }
    std::vector<double> l1x, l1y, l1z, d1;
    int l1xl, l1yl, l1zl, d1l;
    std::vector<double> l2x, l2y, l2z, d2;
    int l2xl, l2yl, l2zl, d2l;
    a.get_exact_lambda(l1x, l1xl, l1y, l1yl, l1z, l1zl, d1, d1l);
    b.get_exact_lambda(l2x, l2xl, l2y, l2yl, l2z, l2zl, d2, d2l);
    if (d1.data()[d1l - 1] == 0.0 || d2.data()[d2l - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xx_exact(l1z, l1zl, d1, d1l, l2z, l2zl, d2, d2l);
}

inline int
less_than_lt_filtered(const double& l1x, const double& d1, const double& l2x, const double& d2, double max_var) {
    double k1 = d2 * l1x;
    double k2 = d1 * l2x;
    double kx = k1 - k2;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 4.321380059346694e-12;
    if (kx > epsilon) return IPSign::POSITIVE;
    if (-kx > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int less_than_on_x_lt(const ImplicitPointLPI& a, const ImplicitPointTPI& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv) && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_lt_filtered(a.ssfilter[0], a.ssfilter[3], b.ssfilter[0], b.ssfilter[3], mv)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda() && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_xx_interval(a.dfilter[0], a.dfilter[3], b.dfilter[0], b.dfilter[3])) != 0) {
            return ret;
        }
    }
    std::vector<double> l1x, l1y, l1z, d1;
    int l1xl, l1yl, l1zl, d1l;
    std::vector<double> l2x, l2y, l2z, d2;
    int l2xl, l2yl, l2zl, d2l;
    a.get_exact_lambda(l1x, l1xl, l1y, l1yl, l1z, l1zl, d1, d1l);
    b.get_exact_lambda(l2x, l2xl, l2y, l2yl, l2z, l2zl, d2, d2l);
    if (d1.data()[d1l - 1] == 0.0 || d2.data()[d2l - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xx_exact(l1x, l1xl, d1, d1l, l2x, l2xl, d2, d2l);
}

inline int less_than_on_y_lt(const ImplicitPointLPI& a, const ImplicitPointTPI& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv) && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_lt_filtered(a.ssfilter[1], a.ssfilter[3], b.ssfilter[1], b.ssfilter[3], mv)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda() && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_xx_interval(a.dfilter[1], a.dfilter[3], b.dfilter[1], b.dfilter[3])) != 0) {
            return ret;
        }
    }
    std::vector<double> l1x, l1y, l1z, d1;
    int l1xl, l1yl, l1zl, d1l;
    std::vector<double> l2x, l2y, l2z, d2;
    int l2xl, l2yl, l2zl, d2l;
    a.get_exact_lambda(l1x, l1xl, l1y, l1yl, l1z, l1zl, d1, d1l);
    b.get_exact_lambda(l2x, l2xl, l2y, l2yl, l2z, l2zl, d2, d2l);
    if (d1.data()[d1l - 1] == 0.0 || d2.data()[d2l - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xx_exact(l1y, l1yl, d1, d1l, l2y, l2yl, d2, d2l);
}

inline int less_than_on_z_lt(const ImplicitPointLPI& a, const ImplicitPointTPI& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv) && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_lt_filtered(a.ssfilter[2], a.ssfilter[3], b.ssfilter[2], b.ssfilter[3], mv)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda() && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_xx_interval(a.dfilter[2], a.dfilter[3], b.dfilter[2], b.dfilter[3])) != 0) {
            return ret;
        }
    }
    std::vector<double> l1x, l1y, l1z, d1;
    int l1xl, l1yl, l1zl, d1l;
    std::vector<double> l2x, l2y, l2z, d2;
    int l2xl, l2yl, l2zl, d2l;
    a.get_exact_lambda(l1x, l1xl, l1y, l1yl, l1z, l1zl, d1, d1l);
    b.get_exact_lambda(l2x, l2xl, l2y, l2yl, l2z, l2zl, d2, d2l);
    if (d1.data()[d1l - 1] == 0.0 || d2.data()[d2l - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xx_exact(l1z, l1zl, d1, d1l, l2z, l2zl, d2, d2l);
}

inline int
less_than_tt_filtered(const double& l1x, const double& d1, const double& l2x, const double& d2, double max_var) {
    double k1 = d2 * l1x;
    double k2 = d1 * l2x;
    double kx = k1 - k2;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= max_var;
    epsilon *= 4.321380059346694e-12;
    if (kx > epsilon) return IPSign::POSITIVE;
    if (-kx > epsilon) return IPSign::NEGATIVE;
    return IPSign::ZERO;
}

inline int less_than_on_x_tt(const ImplicitPointTPI& a, const ImplicitPointTPI& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv) && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_tt_filtered(a.ssfilter[0], a.ssfilter[3], b.ssfilter[0], b.ssfilter[3], mv)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda() && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_xx_interval(a.dfilter[0], a.dfilter[3], b.dfilter[0], b.dfilter[3])) != 0) {
            return ret;
        }
    }
    std::vector<double> l1x, l1y, l1z, d1;
    int l1xl, l1yl, l1zl, d1l;
    std::vector<double> l2x, l2y, l2z, d2;
    int l2xl, l2yl, l2zl, d2l;
    a.get_exact_lambda(l1x, l1xl, l1y, l1yl, l1z, l1zl, d1, d1l);
    b.get_exact_lambda(l2x, l2xl, l2y, l2yl, l2z, l2zl, d2, d2l);
    if (d1.data()[d1l - 1] == 0.0 || d2.data()[d2l - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xx_exact(l1x, l1xl, d1, d1l, l2x, l2xl, d2, d2l);
}

inline int less_than_on_y_tt(const ImplicitPointTPI& a, const ImplicitPointTPI& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv) && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_tt_filtered(a.ssfilter[1], a.ssfilter[3], b.ssfilter[1], b.ssfilter[3], mv)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda() && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_xx_interval(a.dfilter[1], a.dfilter[3], b.dfilter[1], b.dfilter[3])) != 0) {
            return ret;
        }
    }
    std::vector<double> l1x, l1y, l1z, d1;
    int l1xl, l1yl, l1zl, d1l;
    std::vector<double> l2x, l2y, l2z, d2;
    int l2xl, l2yl, l2zl, d2l;
    a.get_exact_lambda(l1x, l1xl, l1y, l1yl, l1z, l1zl, d1, d1l);
    b.get_exact_lambda(l2x, l2xl, l2y, l2yl, l2z, l2zl, d2, d2l);
    if (d1.data()[d1l - 1] == 0.0 || d2.data()[d2l - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xx_exact(l1y, l1yl, d1, d1l, l2y, l2yl, d2, d2l);
}

inline int less_than_on_z_tt(const ImplicitPointTPI& a, const ImplicitPointTPI& b) {
    double mv = 0.0;
    int ret;
    if (a.get_filtered_lambda(mv) && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_tt_filtered(a.ssfilter[2], a.ssfilter[3], b.ssfilter[2], b.ssfilter[3], mv)) != 0) {
            return ret;
        }
    }
    if (a.get_interval_lambda() && b.get_filtered_lambda(mv)) {
        if ((ret = less_than_xx_interval(a.dfilter[2], a.dfilter[3], b.dfilter[2], b.dfilter[3])) != 0) {
            return ret;
        }
    }
    std::vector<double> l1x, l1y, l1z, d1;
    int l1xl, l1yl, l1zl, d1l;
    std::vector<double> l2x, l2y, l2z, d2;
    int l2xl, l2yl, l2zl, d2l;
    a.get_exact_lambda(l1x, l1xl, l1y, l1yl, l1z, l1zl, d1, d1l);
    b.get_exact_lambda(l2x, l2xl, l2y, l2yl, l2z, l2zl, d2, d2l);
    if (d1.data()[d1l - 1] == 0.0 || d2.data()[d2l - 1] == 0.0) {
        return IPSign::UNDEFINED;
    }
    return less_than_xx_exact(l1z, l1zl, d1, d1l, l2z, l2zl, d2, d2l);
}


int GenericPoint3D::less_than_on_x(const GenericPoint3D& a, const GenericPoint3D& b) {
    int type = a.get_type() * 3 + b.get_type();
    switch (type) {
    case 0:
        return (a.to_explicit().x > b.to_explicit().x) - (a.to_explicit().x < b.to_explicit().x);
    case 1: // EL
        return -less_than_on_x_le(b.to_lpi(), a.to_explicit().x);
    case 2: // ET
        return -less_than_on_x_te(b.to_tpi(), a.to_explicit().x);
    case 3: // LE
        return less_than_on_x_le(a.to_lpi(), b.to_explicit().x);
    case 4: // LL
        return less_than_on_x_ll(a.to_lpi(), b.to_lpi());
    case 5: // LT
        return less_than_on_x_lt(a.to_lpi(), b.to_tpi());
    case 6: // TE
        return less_than_on_x_te(a.to_tpi(), b.to_explicit().x);
    case 7: // TL
        return -less_than_on_x_lt(b.to_lpi(), a.to_tpi());
    case 8: // TT
        return less_than_on_x_tt(a.to_tpi(), b.to_tpi());
    default:
        return IPSign::UNDEFINED;
    }
}

int GenericPoint3D::less_than_on_y(const GenericPoint3D& a, const GenericPoint3D& b) {
    int type = a.get_type() * 3 + b.get_type();
    switch (type) {
    case 0:
        return (a.to_explicit().y > b.to_explicit().y) - (a.to_explicit().y < b.to_explicit().y);
    case 1: // EL
        return -less_than_on_y_le(b.to_lpi(), a.to_explicit().y);
    case 2: // ET
        return -less_than_on_y_te(b.to_tpi(), a.to_explicit().y);
    case 3: // LE
        return less_than_on_y_le(a.to_lpi(), b.to_explicit().y);
    case 4: // LL
        return less_than_on_y_ll(a.to_lpi(), b.to_lpi());
    case 5: // LT
        return less_than_on_y_lt(a.to_lpi(), b.to_tpi());
    case 6: // TE
        return less_than_on_y_te(a.to_tpi(), b.to_explicit().y);
    case 7: // TL
        return -less_than_on_y_lt(b.to_lpi(), a.to_tpi());
    case 8: // TT
        return less_than_on_y_tt(a.to_tpi(), b.to_tpi());
    default:
        return IPSign::UNDEFINED;
    }
}

int GenericPoint3D::less_than_on_z(const GenericPoint3D& a, const GenericPoint3D& b) {
    int type = a.get_type() * 3 + b.get_type();
    switch (type) {
    case 0:
        return (a.to_explicit().z > b.to_explicit().z) - (a.to_explicit().z < b.to_explicit().z);
    case 1: // EL
        return -less_than_on_z_le(b.to_lpi(), a.to_explicit().z);
    case 2: // ET
        return -less_than_on_z_te(b.to_tpi(), a.to_explicit().z);
    case 3: // LE
        return less_than_on_z_le(a.to_lpi(), b.to_explicit().z);
    case 4: // LL
        return less_than_on_z_ll(a.to_lpi(), b.to_lpi());
    case 5: // LT
        return less_than_on_z_lt(a.to_lpi(), b.to_tpi());
    case 6: // TE
        return less_than_on_z_te(a.to_tpi(), b.to_explicit().z);
    case 7: // TL
        return -less_than_on_z_lt(b.to_lpi(), a.to_tpi());
    case 8: // TT
        return less_than_on_z_tt(a.to_tpi(), b.to_tpi());
    default:
        return IPSign::UNDEFINED;
    }
}

int GenericPoint3D::orient2d(
    const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c, const int n_max
) {
    if (n_max == 0)
        return orient_yz(a, b, c);
    else if (n_max == 1)
        return orient_zx(a, b, c);
    else
        return orient_xy(a, b, c);
}

inline int max_component_at_triangle_normal_filtered(const double* ov1, const double* ov2, const double* ov3) {
    double v3x = ov3[0] - ov2[0];
    double v3y = ov3[1] - ov2[1];
    double v3z = ov3[2] - ov2[2];
    double v2x = ov2[0] - ov1[0];
    double v2y = ov2[1] - ov1[1];
    double v2z = ov2[2] - ov1[2];
    double nvx1 = v2y * v3z;
    double nvx2 = v2z * v3y;
    double nvx = nvx1 - nvx2;
    double nvy1 = v3x * v2z;
    double nvy2 = v3z * v2x;
    double nvy = nvy1 - nvy2;
    double nvz1 = v2x * v3y;
    double nvz2 = v2y * v3x;
    double nvz = nvz1 - nvz2;

    double _tmp_fabs, max_var = 0;
    if ((_tmp_fabs = fabs(v3x)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(v3y)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(v3z)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(v2x)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(v2y)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = fabs(v2z)) > max_var) max_var = _tmp_fabs;
    double epsilon = 8.88395e-016 * max_var * max_var;

    double nvxc = fabs(nvx);
    double nvyc = fabs(nvy);
    double nvzc = fabs(nvz);
    double nv = nvxc;
    if (nvyc > nv) nv = nvyc;
    if (nvzc > nv) nv = nvzc;

    if (nv > epsilon) {
        if (nv == nvx) return 0;
        if (nv == nvy) return 1;
        if (nv == nvz) return 2;
    }
    return -1;
}
inline int max_component_at_triangle_normal_exact(const double* ov1, const double* ov2, const double* ov3) {
    double v3x[2];
    two_diff(ov3[0], ov2[0], v3x);
    double v3y[2];
    two_diff(ov3[1], ov2[1], v3y);
    double v3z[2];
    two_diff(ov3[2], ov2[2], v3z);
    double v2x[2];
    two_diff(ov2[0], ov1[0], v2x);
    double v2y[2];
    two_diff(ov2[1], ov1[1], v2y);
    double v2z[2];
    two_diff(ov2[2], ov1[2], v2z);
    double nvx1[8];
    two_two_prod(v2y, v3z, nvx1);
    double nvx2[8];
    two_two_prod(v2z, v3y, nvx2);
    double nvx[16];
    int nvx_len = fast_expansion_diff_zeroelim(8, nvx1, 8, nvx2, nvx);
    double nvy1[8];
    two_two_prod(v3x, v2z, nvy1);
    double nvy2[8];
    two_two_prod(v3z, v2x, nvy2);
    double nvy[16];
    int nvy_len = fast_expansion_diff_zeroelim(8, nvy1, 8, nvy2, nvy);
    double nvz1[8];
    two_two_prod(v2x, v3y, nvz1);
    double nvz2[8];
    two_two_prod(v2y, v3x, nvz2);
    double nvz[16];
    int nvz_len = fast_expansion_diff_zeroelim(8, nvz1, 8, nvz2, nvz);

    double nvxc = fabs(nvx[nvx_len - 1]);
    double nvyc = fabs(nvy[nvy_len - 1]);
    double nvzc = fabs(nvz[nvz_len - 1]);
    double nv = nvxc;
    if (nvyc > nv) nv = nvyc;
    if (nvzc > nv) return 2;
    if (nv == nvxc) return 0;
    return 1;
}

int GenericPoint3D::max_component_at_triangle_normal(const double* v1, const double* v2, const double* v3) {
    int ret;
    if ((ret = max_component_at_triangle_normal_filtered(v1, v2, v3)) >= 0) return ret;
    return max_component_at_triangle_normal_exact(v1, v2, v3);
}

int GenericPoint3D::sign_orient3d(const double* v1, const double* v2, const double* v3, const double* v4) {
    const double ret = ::orient3d(v1, v2, v3, v4);
    return (ret > 0.0) - (ret < 0.0);
}

bool GenericPoint3D::inner_segment_cross_inner_triangle(
    const double* u1, const double* u2, const double* v1, const double* v2, const double* v3
) {
    double bound;
    bound = std::fmin(u1[0], u2[0]); // min(u1,u2) alogng x-axis
    if (v1[0] <= bound && v2[0] <= bound && v3[0] <= bound) return false;
    bound = std::fmax(u1[0], u2[0]); // max(u1,u2) alogng x-axis
    if (v1[0] >= bound && v2[0] >= bound && v3[0] >= bound) return false;
    bound = std::fmin(u1[1], u2[1]); // min(u1,u2) alogng y-axis
    if (v1[1] <= bound && v2[1] <= bound && v3[1] <= bound) return false;
    bound = std::fmax(u1[1], u2[1]); // max(u1,u2) alogng y-axis
    if (v1[1] >= bound && v2[1] >= bound && v3[1] >= bound) return false;
    bound = std::fmin(u1[2], u2[2]); // min(u1,u2) alogng z-axis
    if (v1[2] <= bound && v2[2] <= bound && v3[2] <= bound) return false;
    bound = std::fmax(u1[2], u2[2]); // max(u1,u2) alogng z-axis
    if (v1[2] >= bound && v2[2] >= bound && v3[2] >= bound) return false;

    const int orient_u1_tri = sign_orient3d(u1, v1, v2, v3);
    const int orient_u2_tri = sign_orient3d(u2, v1, v2, v3);

    // Check if triangle vertices and at least one of the segment endpoints are coplanar:
    // in this case there is no proper intersection.
    if (orient_u1_tri == 0 || orient_u2_tri == 0) return false;

    // Endpoints of one segment cannot stay both in one of the same half-space defined by the triangle.
    if (orient_u1_tri == orient_u2_tri) return false;

    // Since now, endpoints are one abouve and one below the triangle-plane.

    // Intersection between segment and triangle sides are not proper.
    // Check also if segment intersect the triangle-plane outside the triangle.
    const int orient_u_v1v2 = sign_orient3d(u1, u2, v1, v2);
    const int orient_u_v2v3 = sign_orient3d(u1, u2, v2, v3);

    if (orient_u_v1v2 == 0 || orient_u_v2v3 == 0) return false;
    if (orient_u_v1v2 != orient_u_v2v3) return false;

    const int orient_u_v3v1 = sign_orient3d(u1, u2, v3, v1);

    if (orient_u_v3v1 == 0) return false;
    if (orient_u_v3v1 != orient_u_v2v3) return false;

    // Finally, we have a proper intersection.
    return true;
}

bool GenericPoint3D::inner_segment_cross_triangle(
    const double* u1, const double* u2, const double* v1, const double* v2, const double* v3
) {
    return point_in_inner_segment(v1, u1, u2) || point_in_inner_segment(v2, u1, u2) ||
           point_in_inner_segment(v3, u1, u2) || inner_segments_cross(v2, v3, u1, u2) ||
           inner_segments_cross(v3, v1, u1, u2) || inner_segments_cross(v1, v2, u1, u2) ||
           inner_segment_cross_inner_triangle(u1, u2, v1, v2, v3);
}

bool GenericPoint3D::mis_alignment(const double* p, const double* q, const double* r) {
    // Projection on (x,y)-plane
    if (::orient2d(p, q, r) != .0) return true;

    // Projection on (y,z)-plane
    if (::orient2d(p + 1, q + 1, r + 1) != .0) return true;

    // Projection on (x,z)-plane
    const double pxz[] = {p[0], p[2]};
    const double qxz[] = {q[0], q[2]};
    const double rxz[] = {r[0], r[2]};
    return ::orient2d(pxz, qxz, rxz) != .0;
}

bool GenericPoint3D::mis_alignment(const GenericPoint3D& p, const GenericPoint3D& q, const GenericPoint3D& r) {
    return orient_xy(p, q, r) || orient_yz(p, q, r) || orient_zx(p, q, r);
}

bool GenericPoint3D::same_half_plane(const double* p, const double* q, const double* v1, const double* v2) {

    // Projection on (x,y)-plane
    if (GenericPoint2D::sign_orient2d(p, v1, v2) != GenericPoint2D::sign_orient2d(q, v1, v2)) return false;

    // Projection on (y,z)-plane
    if (GenericPoint2D::sign_orient2d(p + 1, v1 + 1, v2 + 1) != GenericPoint2D::sign_orient2d(q + 1, v1 + 1, v2 + 1))
        return 0;

    // Projection on (x,z)-plane
    const double pxz[] = {p[0], p[2]};
    const double qxz[] = {q[0], q[2]};
    const double v1xz[] = {v1[0], v1[2]};
    const double v2xz[] = {v2[0], v2[2]};
    return (GenericPoint2D::sign_orient2d(pxz, v1xz, v2xz) == GenericPoint2D::sign_orient2d(qxz, v1xz, v2xz));
}

bool GenericPoint3D::point_in_inner_segment(const double* p, const double* v1, const double* v2) {
    return (
        mis_alignment(p, v1, v2) == 0 &&
        ((v1[0] < v2[0] && v1[0] < p[0] && p[0] < v2[0]) || (v1[0] > v2[0] && v1[0] > p[0] && p[0] > v2[0]) ||
         (v1[1] < v2[1] && v1[1] < p[1] && p[1] < v2[1]) || (v1[1] > v2[1] && v1[1] > p[1] && p[1] > v2[1]) ||
         (v1[2] < v2[2] && v1[2] < p[2] && p[2] < v2[2]) || (v1[2] > v2[2] && v1[2] > p[2] && p[2] > v2[2]))
    );
}

bool GenericPoint3D::point_in_inner_segment(
    const GenericPoint3D& p, const GenericPoint3D& v1, const GenericPoint3D& v2, int xyz
) {
    int lt2, lt3;
    if (xyz == 0) {
        if (orient_yz(p, v1, v2)) return false;
        lt2 = less_than_on_y(v1, p);
        lt3 = less_than_on_y(p, v2);
        if (lt2) return (lt2 == lt3);
        lt2 = less_than_on_z(v1, p);
        lt3 = less_than_on_z(p, v2);
    } else if (xyz == 1) {
        if (orient_zx(p, v1, v2)) return false;
        lt2 = less_than_on_x(v1, p);
        lt3 = less_than_on_x(p, v2);
        if (lt2) return (lt2 == lt3);
        lt2 = less_than_on_z(v1, p);
        lt3 = less_than_on_z(p, v2);
    } else {
        if (orient_xy(p, v1, v2)) return false;
        lt2 = less_than_on_x(v1, p);
        lt3 = less_than_on_x(p, v2);
        if (lt2) return (lt2 == lt3);
        lt2 = less_than_on_y(v1, p);
        lt3 = less_than_on_y(p, v2);
    }
    return (lt2 && (lt2 == lt3));
}

bool GenericPoint3D::point_in_segment(const double* p, const double* v1, const double* v2) {
    return same_point(p, v1) || same_point(p, v2) || point_in_inner_segment(p, v1, v2);
}

bool GenericPoint3D::point_in_inner_triangle(const double* p, const double* v1, const double* v2, const double* v3) {
    int o1, o2, oo2, oo4, oo6;

    // Projection on (x,y)-plane -> p VS v1
    o1 = GenericPoint2D::sign_orient2d(p, v2, v3);
    o2 = GenericPoint2D::sign_orient2d(v1, v2, v3);
    oo2 = o2;
    if (o1 != o2) return false;

    // Projection on (y,z)-plane -> p VS v1
    o1 = GenericPoint2D::sign_orient2d(p + 1, v2 + 1, v3 + 1);
    o2 = GenericPoint2D::sign_orient2d(v1 + 1, v2 + 1, v3 + 1);
    oo4 = o2;
    if (o1 != o2) return false;

    // Projection on (x,z)-plane -> p VS v1
    const double pxz[] = {p[0], p[2]};
    const double v1xz[] = {v1[0], v1[2]};
    const double v2xz[] = {v2[0], v2[2]};
    const double v3xz[] = {v3[0], v3[2]};
    o1 = GenericPoint2D::sign_orient2d(pxz, v2xz, v3xz);
    o2 = GenericPoint2D::sign_orient2d(v1xz, v2xz, v3xz);
    oo6 = o2;
    if (o1 != o2) return false;

    // Projection on (x,y)-plane -> p VS v2
    o1 = GenericPoint2D::sign_orient2d(p, v3, v1);
    o2 = oo2;
    if (o1 != o2) return false;

    // Projection on (y,z)-plane -> p VS v2
    o1 = GenericPoint2D::sign_orient2d(p + 1, v3 + 1, v1 + 1);
    o2 = oo4;
    if (o1 != o2) return false;

    // Projection on (x,z)-plane -> p VS v2
    o1 = GenericPoint2D::sign_orient2d(pxz, v3xz, v1xz);
    o2 = oo6;
    if (o1 != o2) return false;

    // Projection on (x,y)-plane -> p VS v3
    o1 = GenericPoint2D::sign_orient2d(p, v1, v2);
    o2 = oo2;
    if (o1 != o2) return false;

    // Projection on (y,z)-plane -> p VS v3
    o1 = GenericPoint2D::sign_orient2d(p + 1, v1 + 1, v2 + 1);
    o2 = oo4;
    if (o1 != o2) return false;

    // Projection on (x,z)-plane -> p VS v3
    o1 = GenericPoint2D::sign_orient2d(pxz, v1xz, v2xz);
    o2 = oo6;
    if (o1 != o2) return false;

    return true;
}
bool GenericPoint3D::point_in_triangle(const double* p, const double* v1, const double* v2, const double* v3) {
    return GenericPoint3D::point_in_segment(p, v1, v2) || GenericPoint3D::point_in_segment(p, v2, v3) ||
           GenericPoint3D::point_in_segment(p, v3, v1) || GenericPoint3D::point_in_inner_triangle(p, v1, v2, v3);
}

bool GenericPoint3D::point_in_triangle(
    const GenericPoint3D& p, const GenericPoint3D& v1, const GenericPoint3D& v2, const GenericPoint3D& v3, const int xyz
) {
    int o1, o2, o3;
    if (xyz == 2) {
        o1 = GenericPoint3D::orient_xy(p, v1, v2);
        o2 = GenericPoint3D::orient_xy(p, v2, v3);
        o3 = GenericPoint3D::orient_xy(p, v3, v1);
    } else if (xyz == 0) {
        o1 = GenericPoint3D::orient_yz(p, v1, v2);
        o2 = GenericPoint3D::orient_yz(p, v2, v3);
        o3 = GenericPoint3D::orient_yz(p, v3, v1);
    } else {
        o1 = GenericPoint3D::orient_zx(p, v1, v2);
        o2 = GenericPoint3D::orient_zx(p, v2, v3);
        o3 = GenericPoint3D::orient_zx(p, v3, v1);
    }
    return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
}

bool GenericPoint3D::inner_segments_cross(const double* u1, const double* u2, const double* v1, const double* v2) {
    // The 4 endpoints must be coplanar
    if (::orient3d(u1, u2, v1, v2) != 0.) return false;

    // Endpoints of one segment cannot stay either on the same side of the other one.
    if (same_half_plane(u1, u2, v1, v2) || same_half_plane(v1, v2, u1, u2)) return false;

    // Each segment endpoint cannot be aligned with the other segment.
    if (!mis_alignment(u1, v1, v2)) return false;
    if (!mis_alignment(u2, v1, v2)) return false;
    if (!mis_alignment(v1, u1, u2)) return false;
    if (!mis_alignment(v2, u1, u2)) return false;

    // If the segment projected on one coordinate plane cross -> segmant cross.
    // Projection on (x,y)-plane
    if (::orient2d(u1, u2, v1) != 0.) return true;
    if (::orient2d(v1, v2, u2) != 0.) return true;
    // Projection on (y,z)-plane
    if (::orient2d(u1 + 1, u2 + 1, v1 + 1) != 0.) return true;
    if (::orient2d(v1 + 1, v2 + 1, u2 + 1) != 0.) return true;
    // Projection on (z,x)-plane
    const double u1xz[] = {u1[0], u1[2]};
    const double u2xz[] = {u2[0], u2[2]};
    const double v1xz[] = {v1[0], v1[2]};
    const double v2xz[] = {v2[0], v2[2]};
    if (::orient2d(u1xz, u2xz, v1xz) != 0.) return true;
    if (::orient2d(v1xz, v2xz, u2xz) != 0.) return true;

    return false;
}

bool GenericPoint3D::inner_segments_cross(
    const GenericPoint3D& u1, const GenericPoint3D& u2, const GenericPoint3D& v1, const GenericPoint3D& v2,
    const int xyz
) {
    int o11, o12, o21, o22;

    if (xyz == 2) {
        o11 = orient_xy(v1, u1, u2);
        o12 = orient_xy(v2, u2, u1);
        o21 = orient_xy(u1, v1, v2);
        o22 = orient_xy(u2, v2, v1);
    } else if (xyz == 0) {
        o11 = orient_yz(v1, u1, u2);
        o12 = orient_yz(v2, u2, u1);
        o21 = orient_yz(u1, v1, v2);
        o22 = orient_yz(u2, v2, v1);
    } else {
        o11 = orient_zx(v1, u1, u2);
        o12 = orient_zx(v2, u2, u1);
        o21 = orient_zx(u1, v1, v2);
        o22 = orient_zx(u2, v2, v1);
    }

    return ((o11 || o21 || o12 || o22) && o11 == o12 && o21 == o22);
}

inline bool lambda3d_LPI_filtered(
    const double* p, const double* q, const double* r, const double* s, const double* t, double* filter
) {
    double a11 = p[0] - q[0];
    double a12 = p[1] - q[1];
    double a13 = p[2] - q[2];
    double a21 = s[0] - r[0];
    double a22 = s[1] - r[1];
    double a23 = s[2] - r[2];
    double a31 = t[0] - r[0];
    double a32 = t[1] - r[1];
    double a33 = t[2] - r[2];
    double tv1 = a22 * a33;
    double tv2 = a23 * a32;
    double a2233 = tv1 - tv2;
    double tv3 = a21 * a33;
    double tv4 = a23 * a31;
    double a2133 = tv3 - tv4;
    double tv5 = a21 * a32;
    double tv6 = a22 * a31;
    double a2132 = tv5 - tv6;
    double tv7 = a11 * a2233;
    double tv8 = a12 * a2133;
    double tv9 = a13 * a2132;
    double tt1 = tv7 - tv8;
    filter[3] = tt1 + tv9;
    double px_rx = p[0] - r[0];
    double py_ry = p[1] - r[1];
    double pz_rz = p[2] - r[2];
    double tt2 = py_ry * a2133;
    double tt3 = px_rx * a2233;
    double tt4 = pz_rz * a2132;
    double tt5 = tt3 + tt4;
    double n = tt5 - tt2;
    double ax = a11 * n;
    double ay = a12 * n;
    double az = a13 * n;
    double dpx = filter[3] * p[0];
    double dpy = filter[3] * p[1];
    double dpz = filter[3] * p[2];
    filter[0] = dpx - ax;
    filter[1] = dpy - ay;
    filter[2] = dpz - az;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(p[0])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(p[1])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(p[2])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(a11)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(a12)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(a13)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(a21)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(a22)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(a23)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(a31)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(a32)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(a33)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(px_rx)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(py_ry)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(pz_rz)) > filter[4]) filter[4] = _tmp_fabs;
    double lambda_d_eps = filter[4];
    lambda_d_eps *= lambda_d_eps;
    lambda_d_eps *= filter[4];
    lambda_d_eps *= 4.884981308350689e-15;

    return ((filter[3] > lambda_d_eps || filter[3] < -lambda_d_eps));
}

bool ImplicitPointLPI::get_filtered_lambda(double& mv) const {
    if (need_filtered_lambda()) {
        if (!lambda3d_LPI_filtered(ip.ptr(), iq.ptr(), ir.ptr(), is.ptr(), it.ptr(), ssfilter.data())) {
            ssfilter[3] = 0.0;
        }
        if (ssfilter[3] < 0) {
            ssfilter[0] = -ssfilter[0];
            ssfilter[1] = -ssfilter[1];
            ssfilter[2] = -ssfilter[2];
            ssfilter[3] = -ssfilter[3];
        }
    }
    if (mv < ssfilter[4]) mv = ssfilter[4];
    return ssfilter[3] != 0.0;
}

inline bool lambda3d_LPI_interval(
    const IntervalNumber& px, const IntervalNumber& py, const IntervalNumber& pz, const IntervalNumber& qx,
    const IntervalNumber& qy, const IntervalNumber& qz, const IntervalNumber& rx, const IntervalNumber& ry,
    const IntervalNumber& rz, const IntervalNumber& sx, const IntervalNumber& sy, const IntervalNumber& sz,
    const IntervalNumber& tx, const IntervalNumber& ty, const IntervalNumber& tz, IntervalNumber* dfilter
) {
    IntervalNumber a11(px - qx);
    IntervalNumber a12(py - qy);
    IntervalNumber a13(pz - qz);
    IntervalNumber a21(sx - rx);
    IntervalNumber a22(sy - ry);
    IntervalNumber a23(sz - rz);
    IntervalNumber a31(tx - rx);
    IntervalNumber a32(ty - ry);
    IntervalNumber a33(tz - rz);
    IntervalNumber tv1(a22 * a33);
    IntervalNumber tv2(a23 * a32);
    IntervalNumber a2233(tv1 - tv2);
    IntervalNumber tv3(a21 * a33);
    IntervalNumber tv4(a23 * a31);
    IntervalNumber a2133(tv3 - tv4);
    IntervalNumber tv5(a21 * a32);
    IntervalNumber tv6(a22 * a31);
    IntervalNumber a2132(tv5 - tv6);
    IntervalNumber tv7(a11 * a2233);
    IntervalNumber tv8(a12 * a2133);
    IntervalNumber tv9(a13 * a2132);
    IntervalNumber tt1(tv7 - tv8);
    dfilter[3] = tt1 + tv9;
    IntervalNumber px_rx(px - rx);
    IntervalNumber py_ry(py - ry);
    IntervalNumber pz_rz(pz - rz);
    IntervalNumber tt2(py_ry * a2133);
    IntervalNumber tt3(px_rx * a2233);
    IntervalNumber tt4(pz_rz * a2132);
    IntervalNumber tt5(tt3 + tt4);
    IntervalNumber n(tt5 - tt2);
    IntervalNumber ax(a11 * n);
    IntervalNumber ay(a12 * n);
    IntervalNumber az(a13 * n);
    IntervalNumber dpx(dfilter[3] * px);
    IntervalNumber dpy(dfilter[3] * py);
    IntervalNumber dpz(dfilter[3] * pz);
    dfilter[0] = dpx - ax;
    dfilter[1] = dpy - ay;
    dfilter[2] = dpz - az;
    return dfilter[3].signIsReliable();
}

bool ImplicitPointLPI::get_interval_lambda() const {
    if (need_interval_lambda()) {
        lambda3d_LPI_interval(
            ip.x, ip.y, ip.z, iq.x, iq.y, iq.z, ir.x, ir.y, ir.z, is.x, is.y, is.z, it.x, it.y, it.z, dfilter.data()
        );
        if (dfilter[3].isNegative()) {
            dfilter[0].invert();
            dfilter[1].invert();
            dfilter[2].invert();
            dfilter[3].invert();
        }
    }
    return dfilter[3].signIsReliable();
}

inline void lambda3d_LPI_exact(
    const double* p, const double* q, const double* r, const double* s, const double* t, std::vector<double>& lambda_x,
    int& lambda_x_len, std::vector<double>& lambda_y, int& lambda_y_len, std::vector<double>& lambda_z,
    int& lambda_z_len, std::vector<double>& lambda_d, int& lambda_d_len
) {
    double a11[2];
    two_diff(p[0], q[0], a11);
    double a12[2];
    two_diff(p[1], q[1], a12);
    double a13[2];
    two_diff(p[2], q[2], a13);
    double a21[2];
    two_diff(s[0], r[0], a21);
    double a22[2];
    two_diff(s[1], r[1], a22);
    double a23[2];
    two_diff(s[2], r[2], a23);
    double a31[2];
    two_diff(t[0], r[0], a31);
    double a32[2];
    two_diff(t[1], r[1], a32);
    double a33[2];
    two_diff(t[2], r[2], a33);
    double tv1[8];
    int tv1_len = product_expansion_zeroelim(2, a22, 2, a33, tv1);
    double tv2[8];
    int tv2_len = product_expansion_zeroelim(2, a23, 2, a32, tv2);
    double a2233[16];
    int a2233_len = fast_expansion_diff_zeroelim(tv1_len, tv1, tv2_len, tv2, a2233);
    double tv3[8];
    int tv3_len = product_expansion_zeroelim(2, a21, 2, a33, tv3);
    double tv4[8];
    int tv4_len = product_expansion_zeroelim(2, a23, 2, a31, tv4);
    double a2133[16];
    int a2133_len = fast_expansion_diff_zeroelim(tv3_len, tv3, tv4_len, tv4, a2133);
    double tv5[8];
    int tv5_len = product_expansion_zeroelim(2, a21, 2, a32, tv5);
    double tv6[8];
    int tv6_len = product_expansion_zeroelim(2, a22, 2, a31, tv6);
    double a2132[16];
    int a2132_len = fast_expansion_diff_zeroelim(tv5_len, tv5, tv6_len, tv6, a2132);
    double tv7[64];
    int tv7_len = product_expansion_zeroelim(2, a11, a2233_len, a2233, tv7);
    double tv8[64];
    int tv8_len = product_expansion_zeroelim(2, a12, a2133_len, a2133, tv8);
    double tv9[64];
    int tv9_len = product_expansion_zeroelim(2, a13, a2132_len, a2132, tv9);
    double tt1[128];
    int tt1_len = fast_expansion_diff_zeroelim(tv7_len, tv7, tv8_len, tv8, tt1);
    lambda_d.resize(static_cast<uint32_t>(tt1_len + tv9_len));
    lambda_d_len = fast_expansion_sum_zeroelim(tt1_len, tt1, tv9_len, tv9, lambda_d.data());
    double px_rx[2];
    two_diff(p[0], r[0], px_rx);
    double py_ry[2];
    two_diff(p[1], r[1], py_ry);
    double pz_rz[2];
    two_diff(p[2], r[2], pz_rz);
    double tt2[64];
    int tt2_len = product_expansion_zeroelim(2, py_ry, a2133_len, a2133, tt2);
    double tt3[64];
    int tt3_len = product_expansion_zeroelim(2, px_rx, a2233_len, a2233, tt3);
    double tt4[64];
    int tt4_len = product_expansion_zeroelim(2, pz_rz, a2132_len, a2132, tt4);
    double tt5[128];
    int tt5_len = fast_expansion_sum_zeroelim(tt3_len, tt3, tt4_len, tt4, tt5);
    std::vector<double> n(static_cast<uint32_t>(tt5_len + tt2_len));
    int n_len = fast_expansion_diff_zeroelim(tt5_len, tt5, tt2_len, tt2, n.data());
    std::vector<double> ax(static_cast<uint32_t>(n_len << 2));
    int ax_len = product_expansion_zeroelim(2, a11, n_len, n.data(), ax.data());
    std::vector<double> ay(static_cast<uint32_t>(n_len << 2));
    int ay_len = product_expansion_zeroelim(2, a12, n_len, n.data(), ay.data());
    std::vector<double> az(static_cast<uint32_t>(n_len << 2));
    int az_len = product_expansion_zeroelim(2, a13, n_len, n.data(), az.data());
    std::vector<double> dpx(static_cast<uint32_t>(lambda_d_len << 1));
    int dpx_len = scale_expansion_zeroelim(lambda_d_len, lambda_d.data(), p[0], dpx.data());
    std::vector<double> dpy(static_cast<uint32_t>(lambda_d_len << 1));
    int dpy_len = scale_expansion_zeroelim(lambda_d_len, lambda_d.data(), p[1], dpy.data());
    std::vector<double> dpz(static_cast<uint32_t>(lambda_d_len << 1));
    int dpz_len = scale_expansion_zeroelim(lambda_d_len, lambda_d.data(), p[2], dpz.data());
    lambda_x.resize(static_cast<uint32_t>(dpx_len + ax_len));
    lambda_x_len = fast_expansion_diff_zeroelim(dpx_len, dpx.data(), ax_len, ax.data(), lambda_x.data());
    lambda_y.resize(static_cast<uint32_t>(dpy_len + ay_len));
    lambda_y_len = fast_expansion_diff_zeroelim(dpy_len, dpy.data(), ay_len, ay.data(), lambda_y.data());
    lambda_z.resize(static_cast<uint32_t>(dpz_len + az_len));
    lambda_z_len = fast_expansion_diff_zeroelim(dpz_len, dpz.data(), az_len, az.data(), lambda_z.data());
}

// Keeps lambda/d pairs as close to one as possible to avoid under/overflows
inline void normalizeLambda3D(double* lx, int lxl, double* ly, int lyl, double* lz, int lzl, double* d, int dl) {
    double maxd = -std::numeric_limits<double>::max(), maxsd = 0.0, ad, aad;
    if ((aad = fabs((ad = estimate(lxl, lx)))) > maxd) {
        maxd = aad;
        maxsd = ad;
    }
    if ((aad = fabs((ad = estimate(lyl, ly)))) > maxd) {
        maxd = aad;
        maxsd = ad;
    }
    if ((aad = fabs((ad = estimate(lzl, lz)))) > maxd) {
        maxd = aad;
        maxsd = ad;
    }
    if ((aad = fabs((ad = estimate(dl, d)))) > maxd) {
        maxd = aad;
        maxsd = ad;
    }

    int e;
    std::frexp(maxsd, &e);
    const double m = std::ldexp(2, -e);
    exact_scale(lxl, lx, m);
    exact_scale(lyl, ly, m);
    exact_scale(lzl, lz, m);
    exact_scale(dl, d, m);
}

void ImplicitPointLPI::get_exact_lambda(
    std::vector<double>& lx, int& lxl, std::vector<double>& ly, int& lyl, std::vector<double>& lz, int& lzl,
    std::vector<double>& d, int& dl
) const {
    lambda3d_LPI_exact(ip.ptr(), iq.ptr(), ir.ptr(), is.ptr(), it.ptr(), lx, lxl, ly, lyl, lz, lzl, d, dl);
    if (d.data()[dl - 1] < 0) {
        invert_expansion(lxl, lx.data());
        invert_expansion(lyl, ly.data());
        invert_expansion(lzl, lz.data());
        invert_expansion(dl, d.data());
    }
    normalizeLambda3D(lx.data(), lxl, ly.data(), lyl, lz.data(), lzl, d.data(), dl);
}

void ImplicitPointLPI::to_double(double* result) {
    std::vector<double> lx, ly, lz, d;
    int lxl, lyl, lzl, dl;
    get_exact_lambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
    const double lambda_x = estimate(lxl, lx.data());
    const double lambda_y = estimate(lyl, ly.data());
    const double lambda_z = estimate(lzl, lz.data());
    const double lambda_d = estimate(dl, d.data());
    result[0] = lambda_x / lambda_d;
    result[1] = lambda_y / lambda_d;
    result[2] = lambda_z / lambda_d;
}

void ImplicitPointLPI::to_double_approx(double* result) {
    double mv = 0.0;
    double lx, ly, lz, d;
    if (get_filtered_lambda(mv)) {
        lx = ssfilter[0];
        ly = ssfilter[1];
        lz = ssfilter[2];
        d = ssfilter[3];
    } else if (get_interval_lambda()) {
        lx = dfilter[0].max() + dfilter[0].min();
        ly = dfilter[1].max() + dfilter[1].min();
        lz = dfilter[2].max() + dfilter[2].min();
        d = dfilter[3].max() + dfilter[3].min();
    } else {
        to_double(result);
        return;
    }
    result[0] = lx / d;
    result[1] = ly / d;
    result[2] = lz / d;
}

inline bool lambda3d_TPI_filtered(
    const double* ov1, const double* ov2, const double* ov3, const double* ow1, const double* ow2, const double* ow3,
    const double* ou1, const double* ou2, const double* ou3, double* filter
) {
    double v3x = ov3[0] - ov2[0];
    double v3y = ov3[1] - ov2[1];
    double v3z = ov3[2] - ov2[2];
    double v2x = ov2[0] - ov1[0];
    double v2y = ov2[1] - ov1[1];
    double v2z = ov2[2] - ov1[2];
    double w3x = ow3[0] - ow2[0];
    double w3y = ow3[1] - ow2[1];
    double w3z = ow3[2] - ow2[2];
    double w2x = ow2[0] - ow1[0];
    double w2y = ow2[1] - ow1[1];
    double w2z = ow2[2] - ow1[2];
    double u3x = ou3[0] - ou2[0];
    double u3y = ou3[1] - ou2[1];
    double u3z = ou3[2] - ou2[2];
    double u2x = ou2[0] - ou1[0];
    double u2y = ou2[1] - ou1[1];
    double u2z = ou2[2] - ou1[2];
    double nvx1 = v2y * v3z;
    double nvx2 = v2z * v3y;
    double nvx = nvx1 - nvx2;
    double nvy1 = v3x * v2z;
    double nvy2 = v3z * v2x;
    double nvy = nvy1 - nvy2;
    double nvz1 = v2x * v3y;
    double nvz2 = v2y * v3x;
    double nvz = nvz1 - nvz2;
    double nwx1 = w2y * w3z;
    double nwx2 = w2z * w3y;
    double nwx = nwx1 - nwx2;
    double nwy1 = w3x * w2z;
    double nwy2 = w3z * w2x;
    double nwy = nwy1 - nwy2;
    double nwz1 = w2x * w3y;
    double nwz2 = w2y * w3x;
    double nwz = nwz1 - nwz2;
    double nux1 = u2y * u3z;
    double nux2 = u2z * u3y;
    double nux = nux1 - nux2;
    double nuy1 = u3x * u2z;
    double nuy2 = u3z * u2x;
    double nuy = nuy1 - nuy2;
    double nuz1 = u2x * u3y;
    double nuz2 = u2y * u3x;
    double nuz = nuz1 - nuz2;
    double nwyuz1 = nwy * nuz;
    double nwyuz2 = nwz * nuy;
    double nwyuz = nwyuz1 - nwyuz2;
    double nwxuz1 = nwx * nuz;
    double nwxuz2 = nwz * nux;
    double nwxuz = nwxuz1 - nwxuz2;
    double nwxuy1 = nwx * nuy;
    double nwxuy2 = nwy * nux;
    double nwxuy = nwxuy1 - nwxuy2;
    double nvyuz1 = nvy * nuz;
    double nvyuz2 = nvz * nuy;
    double nvyuz = nvyuz1 - nvyuz2;
    double nvxuz1 = nvx * nuz;
    double nvxuz2 = nvz * nux;
    double nvxuz = nvxuz1 - nvxuz2;
    double nvxuy1 = nvx * nuy;
    double nvxuy2 = nvy * nux;
    double nvxuy = nvxuy1 - nvxuy2;
    double nvywz1 = nvy * nwz;
    double nvywz2 = nvz * nwy;
    double nvywz = nvywz1 - nvywz2;
    double nvxwz1 = nvx * nwz;
    double nvxwz2 = nvz * nwx;
    double nvxwz = nvxwz1 - nvxwz2;
    double nvxwy1 = nvx * nwy;
    double nvxwy2 = nvy * nwx;
    double nvxwy = nvxwy1 - nvxwy2;
    double p1a = nvx * ov1[0];
    double p1b = nvy * ov1[1];
    double p1c = nvz * ov1[2];
    double p1ab = p1a + p1b;
    double p1 = p1ab + p1c;
    double p2a = nwx * ow1[0];
    double p2b = nwy * ow1[1];
    double p2c = nwz * ow1[2];
    double p2ab = p2a + p2b;
    double p2 = p2ab + p2c;
    double p3a = nux * ou1[0];
    double p3b = nuy * ou1[1];
    double p3c = nuz * ou1[2];
    double p3ab = p3a + p3b;
    double p3 = p3ab + p3c;
    double lxa = p1 * nwyuz;
    double lxb = p3 * nvywz;
    double lxc = p2 * nvyuz;
    double lxab = lxa + lxb;
    filter[0] = lxab - lxc;
    double lya = p2 * nvxuz;
    double lyb = p3 * nvxwz;
    double lyc = p1 * nwxuz;
    double lybc = lyc + lyb;
    filter[1] = lya - lybc;
    double lza = p3 * nvxwy;
    double lzb = p1 * nwxuy;
    double lzc = p2 * nvxuy;
    double lzab = lza + lzb;
    filter[2] = lzab - lzc;
    double da = nvx * nwyuz;
    double db = nvz * nwxuy;
    double dc = nvy * nwxuz;
    double dab = da + db;
    filter[3] = dab - dc;

    double _tmp_fabs;
    if ((_tmp_fabs = fabs(ov1[0])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ov1[1])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ov1[2])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ow1[0])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ow1[1])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ow1[2])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ou1[0])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ou1[1])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(ou1[2])) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(v3x)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(v3y)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(v3z)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(v2x)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(v2y)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(v2z)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(w3x)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(w3y)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(w3z)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(w2x)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(w2y)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(w2z)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(u3x)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(u3y)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(u3z)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(u2x)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(u2y)) > filter[4]) filter[4] = _tmp_fabs;
    if ((_tmp_fabs = fabs(u2z)) > filter[4]) filter[4] = _tmp_fabs;
    double lambda_d_eps = filter[4];
    lambda_d_eps *= lambda_d_eps;
    lambda_d_eps *= lambda_d_eps;
    lambda_d_eps *= filter[4];
    lambda_d_eps *= filter[4];
    lambda_d_eps *= 8.704148513061234e-14;

    return ((filter[3] > lambda_d_eps || filter[3] < -lambda_d_eps));
}

inline bool lambda3d_TPI_interval(
    const IntervalNumber& ov1x, const IntervalNumber& ov1y, const IntervalNumber& ov1z, const IntervalNumber& ov2x,
    const IntervalNumber& ov2y, const IntervalNumber& ov2z, const IntervalNumber& ov3x, const IntervalNumber& ov3y,
    const IntervalNumber& ov3z, const IntervalNumber& ow1x, const IntervalNumber& ow1y, const IntervalNumber& ow1z,
    const IntervalNumber& ow2x, const IntervalNumber& ow2y, const IntervalNumber& ow2z, const IntervalNumber& ow3x,
    const IntervalNumber& ow3y, const IntervalNumber& ow3z, const IntervalNumber& ou1x, const IntervalNumber& ou1y,
    const IntervalNumber& ou1z, const IntervalNumber& ou2x, const IntervalNumber& ou2y, const IntervalNumber& ou2z,
    const IntervalNumber& ou3x, const IntervalNumber& ou3y, const IntervalNumber& ou3z, IntervalNumber* filters
) {
    IntervalNumber v3x(ov3x - ov2x);
    IntervalNumber v3y(ov3y - ov2y);
    IntervalNumber v3z(ov3z - ov2z);
    IntervalNumber v2x(ov2x - ov1x);
    IntervalNumber v2y(ov2y - ov1y);
    IntervalNumber v2z(ov2z - ov1z);
    IntervalNumber w3x(ow3x - ow2x);
    IntervalNumber w3y(ow3y - ow2y);
    IntervalNumber w3z(ow3z - ow2z);
    IntervalNumber w2x(ow2x - ow1x);
    IntervalNumber w2y(ow2y - ow1y);
    IntervalNumber w2z(ow2z - ow1z);
    IntervalNumber u3x(ou3x - ou2x);
    IntervalNumber u3y(ou3y - ou2y);
    IntervalNumber u3z(ou3z - ou2z);
    IntervalNumber u2x(ou2x - ou1x);
    IntervalNumber u2y(ou2y - ou1y);
    IntervalNumber u2z(ou2z - ou1z);
    IntervalNumber nvx1(v2y * v3z);
    IntervalNumber nvx2(v2z * v3y);
    IntervalNumber nvx(nvx1 - nvx2);
    IntervalNumber nvy1(v3x * v2z);
    IntervalNumber nvy2(v3z * v2x);
    IntervalNumber nvy(nvy1 - nvy2);
    IntervalNumber nvz1(v2x * v3y);
    IntervalNumber nvz2(v2y * v3x);
    IntervalNumber nvz(nvz1 - nvz2);
    IntervalNumber nwx1(w2y * w3z);
    IntervalNumber nwx2(w2z * w3y);
    IntervalNumber nwx(nwx1 - nwx2);
    IntervalNumber nwy1(w3x * w2z);
    IntervalNumber nwy2(w3z * w2x);
    IntervalNumber nwy(nwy1 - nwy2);
    IntervalNumber nwz1(w2x * w3y);
    IntervalNumber nwz2(w2y * w3x);
    IntervalNumber nwz(nwz1 - nwz2);
    IntervalNumber nux1(u2y * u3z);
    IntervalNumber nux2(u2z * u3y);
    IntervalNumber nux(nux1 - nux2);
    IntervalNumber nuy1(u3x * u2z);
    IntervalNumber nuy2(u3z * u2x);
    IntervalNumber nuy(nuy1 - nuy2);
    IntervalNumber nuz1(u2x * u3y);
    IntervalNumber nuz2(u2y * u3x);
    IntervalNumber nuz(nuz1 - nuz2);
    IntervalNumber nwyuz1(nwy * nuz);
    IntervalNumber nwyuz2(nwz * nuy);
    IntervalNumber nwyuz(nwyuz1 - nwyuz2);
    IntervalNumber nwxuz1(nwx * nuz);
    IntervalNumber nwxuz2(nwz * nux);
    IntervalNumber nwxuz(nwxuz1 - nwxuz2);
    IntervalNumber nwxuy1(nwx * nuy);
    IntervalNumber nwxuy2(nwy * nux);
    IntervalNumber nwxuy(nwxuy1 - nwxuy2);
    IntervalNumber nvyuz1(nvy * nuz);
    IntervalNumber nvyuz2(nvz * nuy);
    IntervalNumber nvyuz(nvyuz1 - nvyuz2);
    IntervalNumber nvxuz1(nvx * nuz);
    IntervalNumber nvxuz2(nvz * nux);
    IntervalNumber nvxuz(nvxuz1 - nvxuz2);
    IntervalNumber nvxuy1(nvx * nuy);
    IntervalNumber nvxuy2(nvy * nux);
    IntervalNumber nvxuy(nvxuy1 - nvxuy2);
    IntervalNumber nvywz1(nvy * nwz);
    IntervalNumber nvywz2(nvz * nwy);
    IntervalNumber nvywz(nvywz1 - nvywz2);
    IntervalNumber nvxwz1(nvx * nwz);
    IntervalNumber nvxwz2(nvz * nwx);
    IntervalNumber nvxwz(nvxwz1 - nvxwz2);
    IntervalNumber nvxwy1(nvx * nwy);
    IntervalNumber nvxwy2(nvy * nwx);
    IntervalNumber nvxwy(nvxwy1 - nvxwy2);
    IntervalNumber p1a(nvx * ov1x);
    IntervalNumber p1b(nvy * ov1y);
    IntervalNumber p1c(nvz * ov1z);
    IntervalNumber p1ab(p1a + p1b);
    IntervalNumber p1(p1ab + p1c);
    IntervalNumber p2a(nwx * ow1x);
    IntervalNumber p2b(nwy * ow1y);
    IntervalNumber p2c(nwz * ow1z);
    IntervalNumber p2ab(p2a + p2b);
    IntervalNumber p2(p2ab + p2c);
    IntervalNumber p3a(nux * ou1x);
    IntervalNumber p3b(nuy * ou1y);
    IntervalNumber p3c(nuz * ou1z);
    IntervalNumber p3ab(p3a + p3b);
    IntervalNumber p3(p3ab + p3c);
    IntervalNumber lxa(p1 * nwyuz);
    IntervalNumber lxb(p3 * nvywz);
    IntervalNumber lxc(p2 * nvyuz);
    IntervalNumber lxab(lxa + lxb);
    filters[0] = lxab - lxc;
    IntervalNumber lya(p2 * nvxuz);
    IntervalNumber lyb(p3 * nvxwz);
    IntervalNumber lyc(p1 * nwxuz);
    IntervalNumber lybc(lyc + lyb);
    filters[1] = lya - lybc;
    IntervalNumber lza(p3 * nvxwy);
    IntervalNumber lzb(p1 * nwxuy);
    IntervalNumber lzc(p2 * nvxuy);
    IntervalNumber lzab(lza + lzb);
    filters[2] = lzab - lzc;
    IntervalNumber da(nvx * nwyuz);
    IntervalNumber db(nvz * nwxuy);
    IntervalNumber dc(nvy * nwxuz);
    IntervalNumber dab(da + db);
    filters[3] = dab - dc;
    return filters[3].signIsReliable();
}

inline void lambda3d_TPI_exact(
    const double* ov1, const double* ov2, const double* ov3, const double* ow1, const double* ow2, const double* ow3,
    const double* ou1, const double* ou2, const double* ou3, std::vector<double>& lambda_x, int& lambda_x_len,
    std::vector<double>& lambda_y, int& lambda_y_len, std::vector<double>& lambda_z, int& lambda_z_len,
    std::vector<double>& lambda_d, int& lambda_d_len
) {
    double v3x[2];
    two_diff(ov3[0], ov2[0], v3x);
    double v3y[2];
    two_diff(ov3[1], ov2[1], v3y);
    double v3z[2];
    two_diff(ov3[2], ov2[2], v3z);
    double v2x[2];
    two_diff(ov2[0], ov1[0], v2x);
    double v2y[2];
    two_diff(ov2[1], ov1[1], v2y);
    double v2z[2];
    two_diff(ov2[2], ov1[2], v2z);
    double w3x[2];
    two_diff(ow3[0], ow2[0], w3x);
    double w3y[2];
    two_diff(ow3[1], ow2[1], w3y);
    double w3z[2];
    two_diff(ow3[2], ow2[2], w3z);
    double w2x[2];
    two_diff(ow2[0], ow1[0], w2x);
    double w2y[2];
    two_diff(ow2[1], ow1[1], w2y);
    double w2z[2];
    two_diff(ow2[2], ow1[2], w2z);
    double u3x[2];
    two_diff(ou3[0], ou2[0], u3x);
    double u3y[2];
    two_diff(ou3[1], ou2[1], u3y);
    double u3z[2];
    two_diff(ou3[2], ou2[2], u3z);
    double u2x[2];
    two_diff(ou2[0], ou1[0], u2x);
    double u2y[2];
    two_diff(ou2[1], ou1[1], u2y);
    double u2z[2];
    two_diff(ou2[2], ou1[2], u2z);
    double nvx1[8];
    int nvx1_len = product_expansion_zeroelim(2, v2y, 2, v3z, nvx1);
    double nvx2[8];
    int nvx2_len = product_expansion_zeroelim(2, v2z, 2, v3y, nvx2);
    double nvx[16];
    int nvx_len = fast_expansion_diff_zeroelim(nvx1_len, nvx1, nvx2_len, nvx2, nvx);
    double nvy1[8];
    int nvy1_len = product_expansion_zeroelim(2, v3x, 2, v2z, nvy1);
    double nvy2[8];
    int nvy2_len = product_expansion_zeroelim(2, v3z, 2, v2x, nvy2);
    double nvy[16];
    int nvy_len = fast_expansion_diff_zeroelim(nvy1_len, nvy1, nvy2_len, nvy2, nvy);
    double nvz1[8];
    int nvz1_len = product_expansion_zeroelim(2, v2x, 2, v3y, nvz1);
    double nvz2[8];
    int nvz2_len = product_expansion_zeroelim(2, v2y, 2, v3x, nvz2);
    double nvz[16];
    int nvz_len = fast_expansion_diff_zeroelim(nvz1_len, nvz1, nvz2_len, nvz2, nvz);
    double nwx1[8];
    int nwx1_len = product_expansion_zeroelim(2, w2y, 2, w3z, nwx1);
    double nwx2[8];
    int nwx2_len = product_expansion_zeroelim(2, w2z, 2, w3y, nwx2);
    double nwx[16];
    int nwx_len = fast_expansion_diff_zeroelim(nwx1_len, nwx1, nwx2_len, nwx2, nwx);
    double nwy1[8];
    int nwy1_len = product_expansion_zeroelim(2, w3x, 2, w2z, nwy1);
    double nwy2[8];
    int nwy2_len = product_expansion_zeroelim(2, w3z, 2, w2x, nwy2);
    double nwy[16];
    int nwy_len = fast_expansion_diff_zeroelim(nwy1_len, nwy1, nwy2_len, nwy2, nwy);
    double nwz1[8];
    int nwz1_len = product_expansion_zeroelim(2, w2x, 2, w3y, nwz1);
    double nwz2[8];
    int nwz2_len = product_expansion_zeroelim(2, w2y, 2, w3x, nwz2);
    double nwz[16];
    int nwz_len = fast_expansion_diff_zeroelim(nwz1_len, nwz1, nwz2_len, nwz2, nwz);
    double nux1[8];
    int nux1_len = product_expansion_zeroelim(2, u2y, 2, u3z, nux1);
    double nux2[8];
    int nux2_len = product_expansion_zeroelim(2, u2z, 2, u3y, nux2);
    double nux[16];
    int nux_len = fast_expansion_diff_zeroelim(nux1_len, nux1, nux2_len, nux2, nux);
    double nuy1[8];
    int nuy1_len = product_expansion_zeroelim(2, u3x, 2, u2z, nuy1);
    double nuy2[8];
    int nuy2_len = product_expansion_zeroelim(2, u3z, 2, u2x, nuy2);
    double nuy[16];
    int nuy_len = fast_expansion_diff_zeroelim(nuy1_len, nuy1, nuy2_len, nuy2, nuy);
    double nuz1[8];
    int nuz1_len = product_expansion_zeroelim(2, u2x, 2, u3y, nuz1);
    double nuz2[8];
    int nuz2_len = product_expansion_zeroelim(2, u2y, 2, u3x, nuz2);
    double nuz[16];
    int nuz_len = fast_expansion_diff_zeroelim(nuz1_len, nuz1, nuz2_len, nuz2, nuz);
    std::vector<double> nwyuz1(static_cast<uint32_t>((nwy_len * nuz_len) << 1));
    int nwyuz1_len = product_expansion_zeroelim(nwy_len, nwy, nuz_len, nuz, nwyuz1.data());
    std::vector<double> nwyuz2(static_cast<uint32_t>((nwz_len * nuy_len) << 1));
    int nwyuz2_len = product_expansion_zeroelim(nwz_len, nwz, nuy_len, nuy, nwyuz2.data());
    std::vector<double> nwyuz(static_cast<uint32_t>(nwyuz1_len + nwyuz2_len));
    int nwyuz_len = fast_expansion_diff_zeroelim(nwyuz1_len, nwyuz1.data(), nwyuz2_len, nwyuz2.data(), nwyuz.data());
    std::vector<double> nwxuz1(static_cast<uint32_t>((nwx_len * nuz_len) << 1));
    int nwxuz1_len = product_expansion_zeroelim(nwx_len, nwx, nuz_len, nuz, nwxuz1.data());
    std::vector<double> nwxuz2(static_cast<uint32_t>((nwz_len * nux_len) << 1));
    int nwxuz2_len = product_expansion_zeroelim(nwz_len, nwz, nux_len, nux, nwxuz2.data());
    std::vector<double> nwxuz(static_cast<uint32_t>(nwxuz1_len + nwxuz2_len));
    int nwxuz_len = fast_expansion_diff_zeroelim(nwxuz1_len, nwxuz1.data(), nwxuz2_len, nwxuz2.data(), nwxuz.data());
    std::vector<double> nwxuy1(static_cast<uint32_t>((nwx_len * nuy_len) << 1));
    int nwxuy1_len = product_expansion_zeroelim(nwx_len, nwx, nuy_len, nuy, nwxuy1.data());
    std::vector<double> nwxuy2(static_cast<uint32_t>((nwy_len * nux_len) << 1));
    int nwxuy2_len = product_expansion_zeroelim(nwy_len, nwy, nux_len, nux, nwxuy2.data());
    std::vector<double> nwxuy(static_cast<uint32_t>(nwxuy1_len + nwxuy2_len));
    int nwxuy_len = fast_expansion_diff_zeroelim(nwxuy1_len, nwxuy1.data(), nwxuy2_len, nwxuy2.data(), nwxuy.data());
    std::vector<double> nvyuz1(static_cast<uint32_t>((nvy_len * nuz_len) << 1));
    int nvyuz1_len = product_expansion_zeroelim(nvy_len, nvy, nuz_len, nuz, nvyuz1.data());
    std::vector<double> nvyuz2(static_cast<uint32_t>((nvz_len * nuy_len) << 1));
    int nvyuz2_len = product_expansion_zeroelim(nvz_len, nvz, nuy_len, nuy, nvyuz2.data());
    std::vector<double> nvyuz(static_cast<uint32_t>(nvyuz1_len + nvyuz2_len));
    int nvyuz_len = fast_expansion_diff_zeroelim(nvyuz1_len, nvyuz1.data(), nvyuz2_len, nvyuz2.data(), nvyuz.data());
    std::vector<double> nvxuz1(static_cast<uint32_t>((nvx_len * nuz_len) << 1));
    int nvxuz1_len = product_expansion_zeroelim(nvx_len, nvx, nuz_len, nuz, nvxuz1.data());
    std::vector<double> nvxuz2(static_cast<uint32_t>((nvz_len * nux_len) << 1));
    int nvxuz2_len = product_expansion_zeroelim(nvz_len, nvz, nux_len, nux, nvxuz2.data());
    std::vector<double> nvxuz(static_cast<uint32_t>(nvxuz1_len + nvxuz2_len));
    int nvxuz_len = fast_expansion_diff_zeroelim(nvxuz1_len, nvxuz1.data(), nvxuz2_len, nvxuz2.data(), nvxuz.data());
    std::vector<double> nvxuy1(static_cast<uint32_t>((nvx_len * nuy_len) << 1));
    int nvxuy1_len = product_expansion_zeroelim(nvx_len, nvx, nuy_len, nuy, nvxuy1.data());
    std::vector<double> nvxuy2(static_cast<uint32_t>((nvy_len * nux_len) << 1));
    int nvxuy2_len = product_expansion_zeroelim(nvy_len, nvy, nux_len, nux, nvxuy2.data());
    std::vector<double> nvxuy(static_cast<uint32_t>(nvxuy1_len + nvxuy2_len));
    int nvxuy_len = fast_expansion_diff_zeroelim(nvxuy1_len, nvxuy1.data(), nvxuy2_len, nvxuy2.data(), nvxuy.data());
    std::vector<double> nvywz1(static_cast<uint32_t>((nvy_len * nwz_len) << 1));
    int nvywz1_len = product_expansion_zeroelim(nvy_len, nvy, nwz_len, nwz, nvywz1.data());
    std::vector<double> nvywz2(static_cast<uint32_t>((nvz_len * nwy_len) << 1));
    int nvywz2_len = product_expansion_zeroelim(nvz_len, nvz, nwy_len, nwy, nvywz2.data());
    std::vector<double> nvywz(static_cast<uint32_t>(nvywz1_len + nvywz2_len));
    int nvywz_len = fast_expansion_diff_zeroelim(nvywz1_len, nvywz1.data(), nvywz2_len, nvywz2.data(), nvywz.data());
    std::vector<double> nvxwz1(static_cast<uint32_t>((nvx_len * nwz_len) << 1));
    int nvxwz1_len = product_expansion_zeroelim(nvx_len, nvx, nwz_len, nwz, nvxwz1.data());
    std::vector<double> nvxwz2(static_cast<uint32_t>((nvz_len * nwx_len) << 1));
    int nvxwz2_len = product_expansion_zeroelim(nvz_len, nvz, nwx_len, nwx, nvxwz2.data());
    std::vector<double> nvxwz(static_cast<uint32_t>(nvxwz1_len + nvxwz2_len));
    int nvxwz_len = fast_expansion_diff_zeroelim(nvxwz1_len, nvxwz1.data(), nvxwz2_len, nvxwz2.data(), nvxwz.data());
    std::vector<double> nvxwy1(static_cast<uint32_t>((nvx_len * nwy_len) << 1));
    int nvxwy1_len = product_expansion_zeroelim(nvx_len, nvx, nwy_len, nwy, nvxwy1.data());
    std::vector<double> nvxwy2(static_cast<uint32_t>((nvy_len * nwx_len) << 1));
    int nvxwy2_len = product_expansion_zeroelim(nvy_len, nvy, nwx_len, nwx, nvxwy2.data());
    std::vector<double> nvxwy(static_cast<uint32_t>(nvxwy1_len + nvxwy2_len));
    int nvxwy_len = fast_expansion_diff_zeroelim(nvxwy1_len, nvxwy1.data(), nvxwy2_len, nvxwy2.data(), nvxwy.data());
    std::vector<double> p1a(static_cast<uint32_t>(nvx_len << 1));
    int p1a_len = scale_expansion_zeroelim(nvx_len, nvx, ov1[0], p1a.data());
    std::vector<double> p1b(static_cast<uint32_t>(nvy_len << 1));
    int p1b_len = scale_expansion_zeroelim(nvy_len, nvy, ov1[1], p1b.data());
    std::vector<double> p1c(static_cast<uint32_t>(nvz_len << 1));
    int p1c_len = scale_expansion_zeroelim(nvz_len, nvz, ov1[2], p1c.data());
    std::vector<double> p1ab(static_cast<uint32_t>(p1a_len + p1b_len));
    int p1ab_len = fast_expansion_sum_zeroelim(p1a_len, p1a.data(), p1b_len, p1b.data(), p1ab.data());
    std::vector<double> p1(static_cast<uint32_t>(p1ab_len + p1c_len));
    int p1_len = fast_expansion_sum_zeroelim(p1ab_len, p1ab.data(), p1c_len, p1c.data(), p1.data());
    std::vector<double> p2a(static_cast<uint32_t>(nwx_len << 1));
    int p2a_len = scale_expansion_zeroelim(nwx_len, nwx, ow1[0], p2a.data());
    std::vector<double> p2b(static_cast<uint32_t>(nwy_len << 1));
    int p2b_len = scale_expansion_zeroelim(nwy_len, nwy, ow1[1], p2b.data());
    std::vector<double> p2c(static_cast<uint32_t>(nwz_len << 1));
    int p2c_len = scale_expansion_zeroelim(nwz_len, nwz, ow1[2], p2c.data());
    std::vector<double> p2ab(static_cast<uint32_t>(p2a_len + p2b_len));
    int p2ab_len = fast_expansion_sum_zeroelim(p2a_len, p2a.data(), p2b_len, p2b.data(), p2ab.data());
    std::vector<double> p2(static_cast<uint32_t>(p2ab_len + p2c_len));
    int p2_len = fast_expansion_sum_zeroelim(p2ab_len, p2ab.data(), p2c_len, p2c.data(), p2.data());
    std::vector<double> p3a(static_cast<uint32_t>(nux_len << 1));
    int p3a_len = scale_expansion_zeroelim(nux_len, nux, ou1[0], p3a.data());
    std::vector<double> p3b(static_cast<uint32_t>(nuy_len << 1));
    int p3b_len = scale_expansion_zeroelim(nuy_len, nuy, ou1[1], p3b.data());
    std::vector<double> p3c(static_cast<uint32_t>(nuz_len << 1));
    int p3c_len = scale_expansion_zeroelim(nuz_len, nuz, ou1[2], p3c.data());
    std::vector<double> p3ab(static_cast<uint32_t>(p3a_len + p3b_len));
    int p3ab_len = fast_expansion_sum_zeroelim(p3a_len, p3a.data(), p3b_len, p3b.data(), p3ab.data());
    std::vector<double> p3(static_cast<uint32_t>(p3ab_len + p3c_len));
    int p3_len = fast_expansion_sum_zeroelim(p3ab_len, p3ab.data(), p3c_len, p3c.data(), p3.data());
    std::vector<double> lxa(static_cast<uint32_t>((p1_len * nwyuz_len) << 1));
    int lxa_len = product_expansion_zeroelim(p1_len, p1.data(), nwyuz_len, nwyuz.data(), lxa.data());
    std::vector<double> lxb(static_cast<uint32_t>((p3_len * nvywz_len) << 1));
    int lxb_len = product_expansion_zeroelim(p3_len, p3.data(), nvywz_len, nvywz.data(), lxb.data());
    std::vector<double> lxc(static_cast<uint32_t>((p2_len * nvyuz_len) << 1));
    int lxc_len = product_expansion_zeroelim(p2_len, p2.data(), nvyuz_len, nvyuz.data(), lxc.data());
    std::vector<double> lxab(static_cast<uint32_t>(lxa_len + lxb_len));
    int lxab_len = fast_expansion_sum_zeroelim(lxa_len, lxa.data(), lxb_len, lxb.data(), lxab.data());
    lambda_x.resize(static_cast<uint32_t>(lxab_len + lxc_len));
    lambda_x_len = fast_expansion_diff_zeroelim(lxab_len, lxab.data(), lxc_len, lxc.data(), lambda_x.data());
    std::vector<double> lya(static_cast<uint32_t>((p2_len * nvxuz_len) << 1));
    int lya_len = product_expansion_zeroelim(p2_len, p2.data(), nvxuz_len, nvxuz.data(), lya.data());
    std::vector<double> lyb(static_cast<uint32_t>((p3_len * nvxwz_len) << 1));
    int lyb_len = product_expansion_zeroelim(p3_len, p3.data(), nvxwz_len, nvxwz.data(), lyb.data());
    std::vector<double> lyc(static_cast<uint32_t>((p1_len * nwxuz_len) << 1));
    int lyc_len = product_expansion_zeroelim(p1_len, p1.data(), nwxuz_len, nwxuz.data(), lyc.data());
    std::vector<double> lybc(static_cast<uint32_t>(lyc_len + lyb_len));
    int lybc_len = fast_expansion_sum_zeroelim(lyc_len, lyc.data(), lyb_len, lyb.data(), lybc.data());
    lambda_y.resize(static_cast<uint32_t>(lya_len + lybc_len));
    lambda_y_len = fast_expansion_diff_zeroelim(lya_len, lya.data(), lybc_len, lybc.data(), lambda_y.data());
    std::vector<double> lza(static_cast<uint32_t>((p3_len * nvxwy_len) << 1));
    int lza_len = product_expansion_zeroelim(p3_len, p3.data(), nvxwy_len, nvxwy.data(), lza.data());
    std::vector<double> lzb(static_cast<uint32_t>((p1_len * nwxuy_len) << 1));
    int lzb_len = product_expansion_zeroelim(p1_len, p1.data(), nwxuy_len, nwxuy.data(), lzb.data());
    std::vector<double> lzc(static_cast<uint32_t>((p2_len * nvxuy_len) << 1));
    int lzc_len = product_expansion_zeroelim(p2_len, p2.data(), nvxuy_len, nvxuy.data(), lzc.data());
    std::vector<double> lzab(static_cast<uint32_t>(lza_len + lzb_len));
    int lzab_len = fast_expansion_sum_zeroelim(lza_len, lza.data(), lzb_len, lzb.data(), lzab.data());
    lambda_z.resize(static_cast<uint32_t>(lzab_len + lzc_len));
    lambda_z_len = fast_expansion_diff_zeroelim(lzab_len, lzab.data(), lzc_len, lzc.data(), lambda_z.data());
    std::vector<double> da(static_cast<uint32_t>((nvx_len * nwyuz_len) << 1));
    int da_len = product_expansion_zeroelim(nvx_len, nvx, nwyuz_len, nwyuz.data(), da.data());
    std::vector<double> db(static_cast<uint32_t>((nvz_len * nwxuy_len) << 1));
    int db_len = product_expansion_zeroelim(nvz_len, nvz, nwxuy_len, nwxuy.data(), db.data());
    std::vector<double> dc(static_cast<uint32_t>((nvy_len * nwxuz_len) << 1));
    int dc_len = product_expansion_zeroelim(nvy_len, nvy, nwxuz_len, nwxuz.data(), dc.data());
    std::vector<double> dab(static_cast<uint32_t>(da_len + db_len));
    int dab_len = fast_expansion_sum_zeroelim(da_len, da.data(), db_len, db.data(), dab.data());
    lambda_d.resize(static_cast<uint32_t>(dab_len + dc_len));
    lambda_d_len = fast_expansion_diff_zeroelim(dab_len, dab.data(), dc_len, dc.data(), lambda_d.data());
}

bool ImplicitPointTPI::get_filtered_lambda(double& mv) const {
    if (need_filtered_lambda()) {
        if (!lambda3d_TPI_filtered(
                iv1.ptr(), iv2.ptr(), iv3.ptr(), iw1.ptr(), iw2.ptr(), iw3.ptr(), iu1.ptr(), iu2.ptr(), iu3.ptr(),
                ssfilter.data()
            )) {
            ssfilter[3] = 0.0;
        }

        if (ssfilter[3] < 0) {
            ssfilter[0] = -ssfilter[0];
            ssfilter[1] = -ssfilter[1];
            ssfilter[2] = -ssfilter[2];
            ssfilter[3] = -ssfilter[3];
        }
    }
    if (mv < ssfilter[4]) {
        mv = ssfilter[4];
    }
    return ssfilter[3] != 0.0;
}

bool ImplicitPointTPI::get_interval_lambda() const {
    if (need_interval_lambda()) {
        lambda3d_TPI_interval(
            iv1.x, iv1.y, iv1.z, iv2.x, iv2.y, iv2.z, iv3.x, iv3.y, iv3.z, iw1.x, iw1.y, iw1.z, iw2.x, iw2.y, iw2.z,
            iw3.x, iw3.y, iw3.z, iu1.x, iu1.y, iu1.z, iu2.x, iu2.y, iu2.z, iu3.x, iu3.y, iu3.z, dfilter.data()
        );
        if (dfilter[3].isNegative()) {
            dfilter[0].invert();
            dfilter[1].invert();
            dfilter[2].invert();
            dfilter[3].invert();
        }
    }
    return dfilter[3].signIsReliable();
}

void ImplicitPointTPI::get_exact_lambda(
    std::vector<double>& lx, int& lxl, std::vector<double>& ly, int& lyl, std::vector<double>& lz, int& lzl,
    std::vector<double>& d, int& dl
) const {
    lambda3d_TPI_exact(
        iv1.ptr(), iv2.ptr(), iv3.ptr(), iw1.ptr(), iw2.ptr(), iw3.ptr(), iu1.ptr(), iu2.ptr(), iu3.ptr(), lx, lxl, ly,
        lyl, lz, lzl, d, dl
    );
    if (d.data()[dl - 1] < 0) {
        invert_expansion(lxl, lx.data());
        invert_expansion(lyl, ly.data());
        invert_expansion(lzl, lz.data());
        invert_expansion(dl, d.data());
    }
    normalizeLambda3D(lx.data(), lxl, ly.data(), lyl, lz.data(), lzl, d.data(), dl);
}

void ImplicitPointTPI::to_double(double* result) {
    std::vector<double> lx, ly, lz, d;
    int lxl, lyl, lzl, dl;
    get_exact_lambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
    const double lambda_x = estimate(lxl, lx.data());
    const double lambda_y = estimate(lyl, ly.data());
    const double lambda_z = estimate(lzl, lz.data());
    const double lambda_d = estimate(dl, d.data());
    result[0] = lambda_x / lambda_d;
    result[1] = lambda_y / lambda_d;
    result[2] = lambda_z / lambda_d;
}

void ImplicitPointTPI::to_double_approx(double* result) {
    double mv = 0.0;
    double lx, ly, lz, d;
    if (get_filtered_lambda(mv)) {
        lx = ssfilter[0];
        ly = ssfilter[1];
        lz = ssfilter[2];
        d = ssfilter[3];
    } else if (get_interval_lambda()) {
        lx = dfilter[0].max() + dfilter[0].min();
        ly = dfilter[1].max() + dfilter[1].min();
        lz = dfilter[2].max() + dfilter[2].min();
        d = dfilter[3].max() + dfilter[3].min();
    } else {
        to_double_approx(result);
        return;
    }
    result[0] = lx / d;
    result[1] = ly / d;
    result[2] = lz / d;
}
