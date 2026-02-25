#pragma once
#ifndef _MATH_MPFR_H_
#define _MATH_MPFR_H_

#include <mpfr.h>

struct Mpfr;

inline std::ostream& operator<< (std::ostream& stream, const Mpfr& val);
inline Mpfr operator+(const Mpfr& lhs, const Mpfr& rhs);
inline Mpfr operator-(const Mpfr& lhs, const Mpfr& rhs);
inline Mpfr operator/(const Mpfr& lhs, const Mpfr& rhs);
inline Mpfr operator*(const Mpfr& lhs, const Mpfr& rhs);

struct Mpfr {
    mpfr_t v;
    void init()
    {
        mpfr_init2(v, 256);
    }
    inline void copy(const Mpfr& rhs)
    {
        mpfr_set(v, rhs.v, MPFR_RNDN);
    }
    Mpfr() {
        init();
        mpfr_set_d(v, 0, MPFR_RNDN);
    }
    Mpfr(double d) {
        init();
        mpfr_set_d(v, d, MPFR_RNDN);
    }
    Mpfr(const Mpfr& rhs)
    {
        init();
        copy(rhs);
    }
    ~Mpfr()
    {
        mpfr_clear(v);
    }
    const Mpfr& operator=(const Mpfr& rhs)
    {
        if (this != &rhs)
        {
            copy(rhs);
        }
        return *this;
    }
    inline Mpfr operator-() const
    {
        Mpfr tmp;
        mpfr_neg(tmp.v, v, MPFR_RNDN);
        return tmp;
    }
    inline const Mpfr& operator*=(const Mpfr& rhs)
    {
        *this = *this * rhs;
        return *this;
    }
    inline const Mpfr& operator/=(const Mpfr& rhs)
    {
        *this = *this / rhs;
        return *this;
    }
    inline const Mpfr& operator-=(const Mpfr& rhs)
    {
        *this = *this - rhs;
        return *this;
    }
    inline const Mpfr& operator+=(const Mpfr& rhs)
    {
        *this = *this + rhs;
        return *this;
    }
    inline bool operator<(const Mpfr& rhs) const
    {
        return mpfr_cmp(v, rhs.v) < 0;
    }
    inline bool operator>(const Mpfr& rhs) const
    {
        return mpfr_cmp(v, rhs.v) > 0;
    }
    inline bool operator<=(const Mpfr& rhs) const
    {
        return mpfr_cmp(v, rhs.v) <= 0;
    }
    inline bool operator>=(const Mpfr& rhs) const
    {
        return mpfr_cmp(v, rhs.v) >= 0;
    }
    inline Mpfr abs() const
    {
        Mpfr tmp;
        mpfr_abs(tmp.v, v, MPFR_RNDN);
        return tmp;
    }
    inline Mpfr sqrt() const
    {
        Mpfr tmp;
        mpfr_sqrt(tmp.v, v, MPFR_RNDN);
        return tmp;
    }
    static Mpfr tolerance() {
        return Mpfr(1e-64);
    }
    inline bool isZero() const {
        return mpfr_zero_p(v);
    }
    inline bool isInfinite() const {
        return mpfr_inf_p(v);
    }
    inline bool isfinite() const {
        return !isInfinite();
    }
};

inline Mpfr sqrt(const Mpfr & x)
{
    return x.sqrt();
}

inline std::ostream& operator<< (std::ostream& stream, const Mpfr& val)
{
    return stream << mpfr_get_d(val.v, MPFR_RNDN);
}
inline Mpfr operator+(const Mpfr& lhs, const Mpfr& rhs)
{
    Mpfr tmp;
    mpfr_add(tmp.v, lhs.v, rhs.v, MPFR_RNDN);
    return tmp;
}
inline Mpfr operator-(const Mpfr& lhs, const Mpfr& rhs)
{
    Mpfr tmp;
    mpfr_sub(tmp.v, lhs.v, rhs.v, MPFR_RNDN);
    return tmp;
}
inline Mpfr operator/(const Mpfr& lhs, const Mpfr& rhs)
{
    Mpfr tmp;
    mpfr_div(tmp.v, lhs.v, rhs.v, MPFR_RNDN);
    return tmp;
}
inline Mpfr operator*(const Mpfr& lhs, const Mpfr& rhs)
{
    Mpfr tmp;
    mpfr_mul(tmp.v, lhs.v, rhs.v, MPFR_RNDN);
    return tmp;
}
inline double to_double(const Mpfr& rhs)
{
    return mpfr_get_d(rhs.v, MPFR_RNDN);
}

#endif // _MATH_MPFR_H_
