#pragma once
#include "interval_number.h"

#include <array>
#include <vector>

enum IPSign { ZERO = 0, POSITIVE = 1, NEGATIVE = -1, UNDEFINED = 2 };
enum class Point2DType {
    EXPLICIT2D,
    SSI // segment-segment intersection
};

enum class Point3DType {
    EXPLICIT3D = 0,
    LPI = 1, // line-plane intersection
    TPI = 2  // three plane intersection
};

class GenericPoint2D {
  protected:
    Point2DType type;

  public:
    GenericPoint2D(const Point2DType t) : type(t) {}

    bool isExplicit() const { return (type == Point2DType::EXPLICIT2D); }

    const class explicitPoint2D& toExplicit2D() const { return reinterpret_cast<const explicitPoint2D&>(*this); }
    const class implicitPoint2D_SSI& toSSI() const { return reinterpret_cast<const implicitPoint2D_SSI&>(*this); }

    static int sign_orient2d(const double* p, const double* q, const double* r);
    static int orient2D(const GenericPoint2D& a, const GenericPoint2D& b, const GenericPoint2D& c);
};


class explicitPoint2D : public GenericPoint2D {
    double x{0.0}, y{0.0};

  public:
    explicitPoint2D() : GenericPoint2D(Point2DType::EXPLICIT2D) {}
    explicitPoint2D(double _x, double _y) : GenericPoint2D(Point2DType::EXPLICIT2D), x(_x), y(_y) {}
    explicitPoint2D(const explicitPoint2D& b) : GenericPoint2D(Point2DType::EXPLICIT2D), x(b.x), y(b.y) {}

    void operator=(const explicitPoint2D& b) {
        type = Point2DType::EXPLICIT2D;
        x = b.x;
        y = b.y;
    }
    void set(const double a, const double b) {
        x = a;
        y = b;
    }

    double X() const { return x; }
    double Y() const { return y; }

    const double* ptr() const { return &x; }
};

class implicitPoint2D_SSI : public GenericPoint2D {
    const explicitPoint2D &l1_1, &l1_2, &l2_1, &l2_2;

  public:
    // 0: lambda x; 1: lambda y; 2: demininator; 3: max value
    mutable std::array<double, 4> ssfilter;        // semi-static filter
    mutable std::array<IntervalNumber, 3> dfilter; // dynamic filter

    implicitPoint2D_SSI(
        const explicitPoint2D& l11, const explicitPoint2D& l12, const explicitPoint2D& l21, const explicitPoint2D& l22
    )
        : GenericPoint2D(Point2DType::SSI), l1_1(l11), l1_2(l12), l2_1(l21),
          l2_2(l22), ssfilter{
                         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN(), 0.0}} {}

    inline const explicitPoint2D& L1_1() const { return l1_1; }
    inline const explicitPoint2D& L1_2() const { return l1_2; }
    inline const explicitPoint2D& L2_1() const { return l2_1; }
    inline const explicitPoint2D& L2_2() const { return l2_2; }

    // Calculates an explicit approximation of the implicit point.
    // Returns false if point is undefined
    bool approxExplicit(explicitPoint2D&) const;

    // Same as above, but the approximation is as precise as possible.
    // Slightly slower.
    bool apapExplicit(explicitPoint2D&) const;

    bool getFilteredLambda(double& mv) const;
    bool getIntervalLambda() const;
    void getExactLambda(double* lx, int& lxl, double* ly, int& lyl, double* d, int& dl) const;

  private:
    inline bool needsFilteredLambda() const { return (ssfilter[2] != ssfilter[2]); }
    inline bool needsIntervalLambda() const { return (dfilter[2].isNAN()); }
};

class GenericPoint3D {
  protected:
    Point3DType type;

  public:
    GenericPoint3D(const Point3DType t) : type(t) {}
    virtual ~GenericPoint3D() {}
        
    int get_type() const { return static_cast<int>(type); }
    bool is_explicit() const { return type == Point3DType::EXPLICIT3D; }
    bool is_lpi() const { return type == Point3DType::LPI; }
    bool is_tpi() const { return type == Point3DType::TPI; }

    const class ExplicitPoint3D& to_explicit() const { return reinterpret_cast<const ExplicitPoint3D&>(*this); }
    const class ImplicitPointLPI& to_lpi() const { return reinterpret_cast<const ImplicitPointLPI&>(*this); }
    const class ImplicitPointTPI& to_tpi() const { return reinterpret_cast<const ImplicitPointTPI&>(*this); }
    
    virtual void to_double(double*) = 0;
    virtual void to_double_approx(double*) = 0;

    static int orient3d(const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c, const GenericPoint3D& d);

    static int orient_xy(const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c);
    static int orient_yz(const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c);
    static int orient_zx(const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c);
    
    static int orient2d(const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c, const int n_max);

    static int max_component_at_triangle_normal(const double* v1, const double* v2, const double* v3);

    static int sign_orient3d(const double* v1, const double* v2, const double* v3, const double* v4);

    // segment and triangle properly intersects, i.e.
    // intersection occurs in both segment and triangle interior.
    static bool inner_segment_cross_inner_triangle(const double* u1, const double* u2, const double* v1, const double* v2, const double* v3);

    static bool inner_segment_cross_triangle(const double* u1, const double* u2, const double* v1, const double* v2, const double* v3);

    // return true when three points are not aligned
    static bool mis_alignment(const double* p, const double* q, const double* r);

    // return true when points p and q lies both on the same side of the straight line passing through v1 and v2.
    // Note. points and segment must be coplanar
    static bool same_half_plane(const double* p, const double* q, const double* v1, const double* v2);

    // return true when segments properly intersect.
    static bool inner_segments_cross(const double* u1, const double* u2, const double* v1, const double* v2);

    static bool point_in_inner_segment(const double* p, const double* v1, const double* v2);

    static bool point_in_segment(const double* p, const double* v1, const double* v2);

    static bool point_in_inner_triangle(const double* p, const double* v1, const double* v2, const double* v3);

    static bool point_in_triangle(const double* p, const double* v1, const double* v2, const double* v3);
};

class ExplicitPoint3D : public GenericPoint3D {
  public:
    double x, y, z;
    inline ExplicitPoint3D() : GenericPoint3D(Point3DType::EXPLICIT3D), x{0.0}, y{0.0}, z{0.0} {}
    inline ExplicitPoint3D(double _x, double _y, double _z)
        : GenericPoint3D(Point3DType::EXPLICIT3D), x(_x), y(_y), z(_z) {}
    inline ExplicitPoint3D(const ExplicitPoint3D& b) : GenericPoint3D(Point3DType::EXPLICIT3D), x(b.x), y(b.y), z(b.z) {}
    const double* ptr() const { return &x; }
        
    void to_double(double* result) override {
        result[0] = x;
        result[1] = y;
        result[2] = z;
    }
    virtual void to_double_approx(double* result) override {
        result[0] = x;
        result[1] = y;
        result[2] = z;
    }
};

class ImplicitPointLPI : public GenericPoint3D {
    const ExplicitPoint3D &ip, &iq;      // The line
    const ExplicitPoint3D &ir, &is, &it; // The plane

  public:
    mutable std::array<double, 5> ssfilter;        // semi-static filter
    mutable std::array<IntervalNumber, 4> dfilter; // dynamic filter

    ImplicitPointLPI(
        const ExplicitPoint3D& _p, const ExplicitPoint3D& _q, const ExplicitPoint3D& _r, const ExplicitPoint3D& _s,
        const ExplicitPoint3D& _t
    )
        : GenericPoint3D(Point3DType::LPI), ip(_p), iq(_q), ir(_r), is(_s),
          it(_t), ssfilter{
                      {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 0.0}} {}

    bool get_filtered_lambda(double& mv) const;
    bool get_interval_lambda() const;
    void get_exact_lambda(
        std::vector<double>& lx, int& lxl, std::vector<double>& ly, int& lyl, std::vector<double>& lz, int& lzl,
        std::vector<double>& d, int& dl
    ) const;
    void to_double(double*) override;
    virtual void to_double_approx(double*) override;
        
    const ExplicitPoint3D& p() const { return ip; }
    const ExplicitPoint3D& q() const { return iq; }
    const ExplicitPoint3D& r() const { return ir; }
    const ExplicitPoint3D& s() const { return is; }
    const ExplicitPoint3D& t() const { return it; }

  private:
    bool need_filtered_lambda() const { return (ssfilter[3] != ssfilter[3]); } // TRUE if NAN
    bool need_interval_lambda() const { return (dfilter[3].isNAN()); }         // TRUE if NAN
};

class ImplicitPointTPI : public GenericPoint3D {
    const ExplicitPoint3D &iv1, &iv2, &iv3; // Plane 1
    const ExplicitPoint3D &iw1, &iw2, &iw3; // Plane 2
    const ExplicitPoint3D &iu1, &iu2, &iu3; // Plane 3

  public:
    mutable std::array<double, 5> ssfilter;        // semi-static filter
    mutable std::array<IntervalNumber, 4> dfilter; // dynamic filter
    ImplicitPointTPI(
        const ExplicitPoint3D& _v1, const ExplicitPoint3D& _v2, const ExplicitPoint3D& _v3, const ExplicitPoint3D& _w1,
        const ExplicitPoint3D& _w2, const ExplicitPoint3D& _w3, const ExplicitPoint3D& _u1, const ExplicitPoint3D& _u2,
        const ExplicitPoint3D& _u3
    )
        : GenericPoint3D(Point3DType::TPI), iv1(_v1), iv2(_v2), iv3(_v3), iw1(_w1), iw2(_w2), iw3(_w3), iu1(_u1),
          iu2(_u2),
          iu3(_u3), ssfilter{
                        {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 0.0}} {}

    bool get_filtered_lambda(double& mv) const;
    bool get_interval_lambda() const;
    void get_exact_lambda(
        std::vector<double>& lx, int& lxl, std::vector<double>& ly, int& lyl, std::vector<double>& lz, int& lzl,
        std::vector<double>& d, int& dl
    ) const;
    void to_double(double*) override;
    virtual void to_double_approx(double*) override;
        
    const ExplicitPoint3D& v1() const { return iv1; }
    const ExplicitPoint3D& v2() const { return iv2; }
    const ExplicitPoint3D& v3() const { return iv3; }
    const ExplicitPoint3D& w1() const { return iw1; }
    const ExplicitPoint3D& w2() const { return iw2; }
    const ExplicitPoint3D& w3() const { return iw3; }
    const ExplicitPoint3D& u1() const { return iu1; }
    const ExplicitPoint3D& u2() const { return iu2; }
    const ExplicitPoint3D& u3() const { return iu3; }

  private:
    bool need_filtered_lambda() const { return (ssfilter[3] != ssfilter[3]); } // TRUE if NAN
    bool need_interval_lambda() const { return (dfilter[3].isNAN()); }         // TRUE if NAN
};
