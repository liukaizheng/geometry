#pragma once

#include <cmath>
#include <limits>
#include <utility>

static inline double next(const double& a) { return std::nextafter(a, std::numeric_limits<double>::infinity()); }

class IntervalNumber {
  public:
    IntervalNumber() : inf(-std::numeric_limits<double>::quiet_NaN()), sup(std::numeric_limits<double>::quiet_NaN()) {}
    IntervalNumber(const double& a) : inf(-a), sup(a) {}
    IntervalNumber(const double& _inf, const double& _sup) : inf(_inf), sup(_sup) {}
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
            return IntervalNumber(next(std::fmax(inf * b.sup, sup * b.inf)), next(std::fmax(inf * b.inf, sup * b.sup)));
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
        return IntervalNumber{};
    }

    bool isNAN() { return sup != sup; }
    bool isNegative() const { return (sup < 0.0); }
    bool isPosivite() const { return inf < 0.0; }
    bool signIsReliable() const { return inf < 0.0 || sup < 0.0; }
    int sign() const { return inf < 0.0 ? 1 : (sup < 0.0 ? -1 : 0); }
    void invert() { std::swap(inf, sup); }


  private:
    double inf;
    double sup;
};
