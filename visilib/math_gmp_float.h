#pragma once

#include <gmp.h>
#include <cmath>
struct MathGmpFloat;

inline std::ostream& operator<< (std::ostream& stream, const MathGmpFloat& val);
inline MathGmpFloat operator+(const MathGmpFloat& lhs, const MathGmpFloat& rhs);
inline MathGmpFloat operator-(const MathGmpFloat& lhs, const MathGmpFloat& rhs);
inline MathGmpFloat operator/(const MathGmpFloat& lhs, const MathGmpFloat& rhs);
inline MathGmpFloat operator*(const MathGmpFloat& lhs, const MathGmpFloat& rhs);

struct MathGmpFloat {
    mpf_t v;
    void init()
    {
        mpf_init2(v, 256);
    }
    inline void copy(const MathGmpFloat& rhs)
    {
        mpf_set(v, rhs.v);
    }
    MathGmpFloat() {
        init();
        mpf_set_d(v, 0);
    }
    MathGmpFloat(double d) {
        init();
        mpf_set_d(v, d);
    }
    MathGmpFloat(const MathGmpFloat& rhs)
    {
        init();
        copy(rhs);
    }
    ~MathGmpFloat()
    {
        mpf_clear(v);
    }
    const MathGmpFloat& operator=(const MathGmpFloat& rhs)
    {
        if (this != &rhs)
        {
            copy(rhs);
        }
        return *this;
    }
    inline MathGmpFloat operator-() const
    {
        MathGmpFloat tmp;
        mpf_neg(tmp.v, v);
        return tmp;
    }
    inline const MathGmpFloat& operator*=(const MathGmpFloat& rhs)
    {
        *this = *this * rhs;
        return *this;
    }
    inline const MathGmpFloat& operator/=(const MathGmpFloat& rhs)
    {
        *this = *this / rhs;
        return *this;
    }
    inline const MathGmpFloat& operator-=(const MathGmpFloat& rhs)
    {
        *this = *this - rhs;
        return *this;
    }
    inline const MathGmpFloat& operator+=(const MathGmpFloat& rhs)
    {
        *this = *this + rhs;
        return *this;
    }
    inline bool operator<(const MathGmpFloat& rhs) const
    {
        return mpf_cmp(v, rhs.v) < 0;
    }
    inline bool operator>(const MathGmpFloat& rhs) const
    {
        return mpf_cmp(v, rhs.v) > 0;
    }
    inline bool operator<=(const MathGmpFloat& rhs) const
    {
        return mpf_cmp(v, rhs.v) <= 0;
    }
    inline bool operator>=(const MathGmpFloat& rhs) const
    {
        return mpf_cmp(v, rhs.v) >= 0;
    }
    inline MathGmpFloat abs() const 
    {
        MathGmpFloat tmp;
        mpf_abs(tmp.v, v);
        return tmp;
    }
    inline MathGmpFloat sqrt() const
    {
        MathGmpFloat tmp;
        mpf_sqrt(tmp.v, v);
        return tmp;
    }
    inline bool isInfinite() const {
        return false;
    }
    inline bool isfinite() const {
        return !isInfinite();
    }
    inline static mp_bitcnt_t get_default_precision()
    {
      return mpf_get_default_prec();
    }
    inline static void set_default_precision(mp_bitcnt_t prec)
    {
        mpf_set_default_prec(prec);
    }

    inline static mp_bitcnt_t digits2bits(int d)
    {
        const double LOG2_10 = 3.3219280948873624;
        return mp_bitcnt_t(std::ceil( d * LOG2_10 ));
    }

};

inline MathGmpFloat sqrt(const MathGmpFloat & x)
{
    return x.sqrt();
}

inline std::ostream& operator<< (std::ostream& stream, const MathGmpFloat& val)
{
    return stream << mpf_get_d(val.v);
}
inline MathGmpFloat operator+(const MathGmpFloat& lhs, const MathGmpFloat& rhs)
{
    MathGmpFloat tmp;
    mpf_add(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline MathGmpFloat operator-(const MathGmpFloat& lhs, const MathGmpFloat& rhs)
{
    MathGmpFloat tmp;
    mpf_sub(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline MathGmpFloat operator/(const MathGmpFloat& lhs, const MathGmpFloat& rhs)
{
    MathGmpFloat tmp;
    mpf_div(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline MathGmpFloat operator*(const MathGmpFloat& lhs, const MathGmpFloat& rhs)
{
    MathGmpFloat tmp;
    mpf_mul(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline double to_double(const MathGmpFloat& rhs)
{
    return mpf_get_d(rhs.v);
}
