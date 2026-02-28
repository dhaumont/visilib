#pragma once

#include <mpfr.h>
#include <cmath>
struct MathMpfr;

inline std::ostream& operator<< (std::ostream& stream, const MathMpfr& val);
inline MathMpfr operator+(const MathMpfr& lhs, const MathMpfr& rhs);
inline MathMpfr operator-(const MathMpfr& lhs, const MathMpfr& rhs);
inline MathMpfr operator/(const MathMpfr& lhs, const MathMpfr& rhs);
inline MathMpfr operator*(const MathMpfr& lhs, const MathMpfr& rhs);

struct MathMpfr {
    mpfr_t v;
    void init()
    {
        mpfr_init2(v, get_default_precision());
    }
    inline void copy(const MathMpfr& rhs)
    {
        mpfr_set(v, rhs.v, MPFR_RNDN);
    }
    MathMpfr() {
        init();
        mpfr_set_d(v, 0, MPFR_RNDN);
    }
    MathMpfr(double d) {
        init();
        mpfr_set_d(v, d, MPFR_RNDN);
    }
    MathMpfr(const MathMpfr& rhs)
    {
        init();
        copy(rhs);
    }
    ~MathMpfr()
    {
        mpfr_clear(v);
    }
    const MathMpfr& operator=(const MathMpfr& rhs)
    {
        if (this != &rhs)
        {
            copy(rhs);
        }
        return *this;
    }
    inline MathMpfr operator-() const
    {
        MathMpfr tmp;
        mpfr_neg(tmp.v, v, MPFR_RNDN);
        return tmp;
    }
    inline const MathMpfr& operator*=(const MathMpfr& rhs)
    {
        *this = *this * rhs;
        return *this;
    }
    inline const MathMpfr& operator/=(const MathMpfr& rhs)
    {
        *this = *this / rhs;
        return *this;
    }
    inline const MathMpfr& operator-=(const MathMpfr& rhs)
    {
        *this = *this - rhs;
        return *this;
    }
    inline const MathMpfr& operator+=(const MathMpfr& rhs)
    {
        *this = *this + rhs;
        return *this;
    }
    inline bool operator<(const MathMpfr& rhs) const
    {
        return mpfr_cmp(v, rhs.v) < 0;
    }
    inline bool operator>(const MathMpfr& rhs) const
    {
        return mpfr_cmp(v, rhs.v) > 0;
    }
    inline bool operator<=(const MathMpfr& rhs) const
    {
        return mpfr_cmp(v, rhs.v) <= 0;
    }
    inline bool operator>=(const MathMpfr& rhs) const
    {
        return mpfr_cmp(v, rhs.v) >= 0;
    }
    inline MathMpfr abs() const
    {
        MathMpfr tmp;
        mpfr_abs(tmp.v, v, MPFR_RNDN);
        return tmp;
    }
    inline MathMpfr sqrt() const
    {
        MathMpfr tmp;
        mpfr_sqrt(tmp.v, v, MPFR_RNDN);
        return tmp;
    }
    inline bool isInfinite() const {
        return mpfr_inf_p(v);
    }
    inline bool isfinite() const {
        return !isInfinite();
    }
    inline static mp_rnd_t get_default_rounding_mode()
    {
       return (mp_rnd_t)(mpfr_get_default_rounding_mode());
    }
    inline static mp_prec_t get_default_precision()
    {
      return mpfr_get_default_prec();
    }
    inline static void set_default_precision(mp_prec_t prec)
    {
        mpfr_set_default_prec(prec);
    }
    inline static void set_default_rounding_mode(mp_rnd_t rnd_mode)
    {
        mpfr_set_default_rounding_mode(rnd_mode);
    }

    inline static mp_prec_t digits2bits(int d)
    {
        const double LOG2_10 = 3.3219280948873624;
        return mp_prec_t(std::ceil( d * LOG2_10 ));
    }

};

inline std::ostream& operator<< (std::ostream& stream, const MathMpfr& val)
{
    return stream << mpfr_get_d(val.v, MPFR_RNDN);
}
inline MathMpfr operator+(const MathMpfr& lhs, const MathMpfr& rhs)
{
    MathMpfr tmp;
    mpfr_add(tmp.v, lhs.v, rhs.v, MPFR_RNDN);
    return tmp;
}
inline MathMpfr operator-(const MathMpfr& lhs, const MathMpfr& rhs)
{
    MathMpfr tmp;
    mpfr_sub(tmp.v, lhs.v, rhs.v, MPFR_RNDN);
    return tmp;
}
inline MathMpfr operator/(const MathMpfr& lhs, const MathMpfr& rhs)
{
    MathMpfr tmp;
    mpfr_div(tmp.v, lhs.v, rhs.v, MPFR_RNDN);
    return tmp;
}
inline MathMpfr operator*(const MathMpfr& lhs, const MathMpfr& rhs)
{
    MathMpfr tmp;
    mpfr_mul(tmp.v, lhs.v, rhs.v, MPFR_RNDN);
    return tmp;
}
inline double to_double(const MathMpfr& rhs)
{
    return mpfr_get_d(rhs.v, MPFR_RNDN);
}

