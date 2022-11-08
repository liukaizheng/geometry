#pragma once

#include <cmath>
#include <limits>

static inline double next(const double& a) { return std::nextafter(a, std::numeric_limits<double>::infinity()); }

class IntervalNumber {
  public:
    IntervalNumber() : inf(-NAN), sup(NAN) {}
    IntervalNumber(const double& a) : inf(-a), sup(a) {}
    IntervalNumber(const double& inf, const double& sup) : inf(inf), sup(sup) {}
    IntervalNumber& operator=(const IntervalNumber& b) {
        inf = b.inf;
        sup = b.sup;
        return *this;
    }
    IntervalNumber operator+(const IntervalNumber& b) const {
        return IntervalNumber(next(inf + b.inf), next(sup + b.sup));
    }
    IntervalNumber operator-(const IntervalNumber& b) const {
        return IntervalNumber(next(inf + b.sup), next(sup + b.inf));
    }
    IntervalNumber operator*(const IntervalNumber& b) const {
        const auto sign = [](const IntervalNumber& num) { return std::signbit(num.inf) << 1 | std::signbit(num.sup); };
        switch (sign(*this) << 2 | sign(b)) {
        case 0:
            return IntervalNumber(next(std::max(inf * b.sup, sup * b.inf)), next(std::max(inf * b.inf, sup * b.sup)));
        case 1:
            return IntervalNumber(next(sup * b.inf), next(inf * b.inf));
        case 2:
            return IntervalNumber(next(inf * b.sup), next(sup * b.sup));
        case 4:
            return IntervalNumber(next(inf * b.sup), next(inf * b.inf));
        case 5:
            return IntervalNumber(next((-sup) * b.sup), next(inf * b.inf));
        case 6:
            return IntervalNumber(next(inf * b.sup), next(sup * (-b.inf)));
        case 8:
            return IntervalNumber(next(sup * b.inf), next(sup * b.sup));
        case 9:
            return IntervalNumber(next(sup * b.inf), next((-inf) * b.sup));
        case 10:
            return IntervalNumber(next(inf * (-b.inf)), next(sup * b.sup));
        default:
            break;
        }
        return IntervalNumber(NAN);
    }

  private:
    double inf;
    double sup;
};
