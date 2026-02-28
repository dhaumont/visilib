#pragma once

#include <gmp.h>
struct MathGmpRational;

inline std::ostream& operator<< (std::ostream& stream, const MathGmpRational& val);
inline MathGmpRational operator+(const MathGmpRational& lhs, const MathGmpRational& rhs);
inline MathGmpRational operator-(const MathGmpRational& lhs, const MathGmpRational& rhs);
inline MathGmpRational operator/(const MathGmpRational& lhs, const MathGmpRational& rhs);
inline MathGmpRational operator*(const MathGmpRational& lhs, const MathGmpRational& rhs);

struct MathGmpRational {
    mpq_t v;
    void init()
    {
        mpq_init(v);
    }
    inline void copy(const MathGmpRational& rhs)
    {
        mpq_set(v, rhs.v);
    }
    MathGmpRational() {
        init();
        mpq_set_d(v, 0);
    }
    MathGmpRational(double d) {
        init();
        mpq_set_d(v, d);
    }
    MathGmpRational(const MathGmpRational& rhs)
    {
        init();
        copy(rhs);
    }
    ~MathGmpRational()
    {
        mpq_clear(v);
    }
    const MathGmpRational& operator=(const MathGmpRational& rhs)
    {
        if (this != &rhs)
        {
            copy(rhs);
        }
        return *this;
    }
    inline MathGmpRational operator-() const
    {
        MathGmpRational tmp;
        mpq_neg(tmp.v, v);
        return tmp;
    }
    inline const MathGmpRational& operator*=(const MathGmpRational& rhs)
    {
        *this = *this * rhs;
        return *this;
    }
    inline const MathGmpRational& operator/=(const MathGmpRational& rhs)
    {
        *this = *this / rhs;
        return *this;
    }
    inline const MathGmpRational& operator-=(const MathGmpRational& rhs)
    {
        *this = *this - rhs;
        return *this;
    }
    inline const MathGmpRational& operator+=(const MathGmpRational& rhs)
    {
        *this = *this + rhs;
        return *this;
    }
    inline bool operator<(const MathGmpRational& rhs) const
    {
        return mpq_cmp(v, rhs.v) < 0;
    }
    inline bool operator>(const MathGmpRational& rhs) const
    {
        return mpq_cmp(v, rhs.v) > 0;
    }
    inline bool operator<=(const MathGmpRational& rhs) const
    {
        return mpq_cmp(v, rhs.v) <= 0;
    }
    inline bool operator>=(const MathGmpRational& rhs) const
    {
        return mpq_cmp(v, rhs.v) >= 0;
    }
    inline MathGmpRational abs() const 
    {
        MathGmpRational tmp;
        mpq_abs(tmp.v, v);
        return tmp;
    }
    inline MathGmpRational sqrt() const 
    {
        MathGmpRational tmp;
        mpq_set_d(tmp.v, std::sqrt(mpq_get_d(v)));
        return tmp;
    }
    inline bool isInfinite() const {
        return false;
    }
    inline bool isfinite() const {
        return !isInfinite();
    }
};

inline std::ostream& operator<< (std::ostream& stream, const MathGmpRational& val)
{
    return stream << mpq_get_d(val.v);
}
inline MathGmpRational operator+(const MathGmpRational& lhs, const MathGmpRational& rhs)
{
    MathGmpRational tmp;
    mpq_add(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline MathGmpRational operator-(const MathGmpRational& lhs, const MathGmpRational& rhs)
{
    MathGmpRational tmp;
    mpq_sub(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline MathGmpRational operator/(const MathGmpRational& lhs, const MathGmpRational& rhs)
{
    MathGmpRational tmp;
    mpq_div(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline MathGmpRational operator*(const MathGmpRational& lhs, const MathGmpRational& rhs)
{
    MathGmpRational tmp;
    mpq_mul(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline double to_double(const MathGmpRational& rhs)
{
    return mpq_get_d(rhs.v);
}

