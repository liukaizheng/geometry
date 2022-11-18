#pragma once
#include "interval_number.h"

#include <array>

enum IPSign { ZERO = 0, POSITIVE = 1, NEGATIVE = -1, UNDEFINED = 2 };
enum class Point2DType {
    EXPLICIT2D,
    SSI // segment-segment intersection
};

enum class Point3DType {
    EXPLICIT3D = 3,
    LPI = 4, // line-plane intersection
    TPI = 5  // three plane intersection
};

class GenericPoint2D {
  protected:
    Point2DType type;

  public:
    GenericPoint2D(const Point2DType t) : type(t) {}

    bool isExplicit() const { return (type == Point2DType::EXPLICIT2D); }

    const class explicitPoint2D& toExplicit2D() const { return reinterpret_cast<const explicitPoint2D&>(*this); }
    const class implicitPoint2D_SSI& toSSI() const { return reinterpret_cast<const implicitPoint2D_SSI&>(*this); }
    const class explicitPoint3D& toExplicit3D() const { return reinterpret_cast<const explicitPoint3D&>(*this); }
    const class implicitPoint3D_LPI& toLPI() const { return reinterpret_cast<const implicitPoint3D_LPI&>(*this); }
    const class implicitPoint3D_TPI& toTPI() const { return reinterpret_cast<const implicitPoint3D_TPI&>(*this); }
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
    mutable std::array<double, 4> ssfilter;             // semi-static filter
    mutable std::array<IntervalNumber, 3> dfilter;      // dynamic filter

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

    inline bool needsFilteredLambda() const { return (ssfilter[2] != ssfilter[2]); }
    inline bool needsIntervalLambda() const { return (dfilter[2].isNAN()); }

    bool getFilteredLambda(double& mv) const;
    bool getIntervalLambda() const;
    void getExactLambda(double* lx, int& lxl, double* ly, int& lyl, double* d, int& dl) const;
};


