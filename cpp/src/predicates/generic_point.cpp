#include <predicates/generic_point.h>
#include <predicates/predicates.h>

#include <vector>

inline int orient2d_EEE(const explicitPoint2D& p1, const explicitPoint2D& p2, const explicitPoint2D& p3) {
    const double ret =
        orient2d(const_cast<double*>(p1.ptr()), const_cast<double*>(p2.ptr()), const_cast<double*>(p3.ptr()));
    return (ret > 0) - (ret < 0);
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
    IntervalNumber l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
    if (!p1.getIntervalLambda() || !p2.getIntervalLambda() || !p3.getIntervalLambda()) return IPSign::ZERO;

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

int genericPoint::orient2D(const genericPoint& a, const genericPoint& b, const genericPoint& c) {
    switch (a.isExplicit2D() << 2 | b.isExplicit2D() << 1 | c.isExplicit2D()) {
    case 0:
        return orient2d_SSS(a.toSSI(), b.toSSI(), c.toSSI());
    case 7:
        return orient2d_EEE(a.toExplicit2D(), b.toExplicit2D(), c.toExplicit2D());
    }
    return 0;
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
