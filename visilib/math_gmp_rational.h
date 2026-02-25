#pragma once
#ifndef _MATH_GMP_RATIONAL_H_
#define _MATH_GMP_RATIONAL_H_

#include <gmp.h>
struct GmpRational;

inline std::ostream& operator<< (std::ostream& stream, const GmpRational& val);
inline GmpRational operator+(const GmpRational& lhs, const GmpRational& rhs);
inline GmpRational operator-(const GmpRational& lhs, const GmpRational& rhs);
inline GmpRational operator/(const GmpRational& lhs, const GmpRational& rhs);
inline GmpRational operator*(const GmpRational& lhs, const GmpRational& rhs);

struct GmpRational {
    mpq_t v;
    void init()
    {
        mpq_init(v);
    }
    inline void copy(const GmpRational& rhs)
    {
        mpq_set(v, rhs.v);
    }
    GmpRational() {
        init();
        mpq_set_d(v, 0);
    }
    GmpRational(double d) {
        init();
        mpq_set_d(v, d);
    }
    GmpRational(const GmpRational& rhs)
    {
        init();
        copy(rhs);
    }
    ~GmpRational()
    {
        mpq_clear(v);
    }
    const GmpRational& operator=(const GmpRational& rhs)
    {
        if (this != &rhs)
        {
            copy(rhs);
        }
        return *this;
    }
    inline GmpRational operator-() const
    {
        GmpRational tmp;
        mpq_neg(tmp.v, v);
        return tmp;
    }
    inline const GmpRational& operator*=(const GmpRational& rhs)
    {
        *this = *this * rhs;
        return *this;
    }
    inline const GmpRational& operator/=(const GmpRational& rhs)
    {
        *this = *this / rhs;
        return *this;
    }
    inline const GmpRational& operator-=(const GmpRational& rhs)
    {
        *this = *this - rhs;
        return *this;
    }
    inline const GmpRational& operator+=(const GmpRational& rhs)
    {
        *this = *this + rhs;
        return *this;
    }
    inline bool operator<(const GmpRational& rhs) const
    {
        return mpq_cmp(v, rhs.v) < 0;
    }
    inline bool operator>(const GmpRational& rhs) const
    {
        return mpq_cmp(v, rhs.v) > 0;
    }
    inline bool operator<=(const GmpRational& rhs) const
    {
        return mpq_cmp(v, rhs.v) <= 0;
    }
    inline bool operator>=(const GmpRational& rhs) const
    {
        return mpq_cmp(v, rhs.v) >= 0;
    }
    inline GmpRational abs() const 
    {
        GmpRational tmp;
        mpq_abs(tmp.v, v);
        return tmp;
    }
    inline GmpRational sqrt() const 
    {
        GmpRational tmp;
        mpq_set_d(tmp.v, std::sqrt(mpq_get_d(v)));
        return tmp;
    }
    static GmpRational tolerance() {
        return GmpRational(1e-11);
    }
    inline bool isInfinite() const {
        return false;
    }
    inline bool isfinite() const {
        return !isInfinite();
    }
};

inline GmpRational sqrt(const GmpRational & x)
{
    return x.sqrt();
}

inline std::ostream& operator<< (std::ostream& stream, const GmpRational& val)
{
    return stream << mpq_get_d(val.v);
}
inline GmpRational operator+(const GmpRational& lhs, const GmpRational& rhs)
{
    GmpRational tmp;
    mpq_add(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline GmpRational operator-(const GmpRational& lhs, const GmpRational& rhs)
{
    GmpRational tmp;
    mpq_sub(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline GmpRational operator/(const GmpRational& lhs, const GmpRational& rhs)
{
    GmpRational tmp;
    mpq_div(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline GmpRational operator*(const GmpRational& lhs, const GmpRational& rhs)
{
    GmpRational tmp;
    mpq_mul(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline double to_double(const GmpRational& rhs)
{
    return mpq_get_d(rhs.v);
}

#endif _MATH_GMP_RATIONAL_H_
