#pragma once

#include <gmp.h>
struct MathGmpRational;

inline std::ostream &operator<<(std::ostream &stream, const MathGmpRational &val);
inline MathGmpRational operator+(const MathGmpRational &lhs, const MathGmpRational &rhs);
inline MathGmpRational operator-(const MathGmpRational &lhs, const MathGmpRational &rhs);
inline MathGmpRational operator/(const MathGmpRational &lhs, const MathGmpRational &rhs);
inline MathGmpRational operator*(const MathGmpRational &lhs, const MathGmpRational &rhs);

struct MathGmpRational
{
    mpq_t v;
    void init()
    {
        mpq_init(v);
    }
    inline void copy(const MathGmpRational &rhs)
    {
        mpq_set(v, rhs.v);
    }
    MathGmpRational()
    {
        init();
        mpq_set_d(v, 0);
    }
    MathGmpRational(double d)
    {
        init();
        mpq_set_d(v, d);
    }
    MathGmpRational(const MathGmpRational &rhs)
    {
        init();
        copy(rhs);
    }
    ~MathGmpRational()
    {
        mpq_clear(v);
    }
    const MathGmpRational &operator=(const MathGmpRational &rhs)
    {
        if (this != &rhs)
        {
            copy(rhs);
        }
        return *this;
    }
    inline MathGmpRational operator-() const
    {
        MathGmpRational result;
        mpq_neg(result.v, v);
        return result;
    }
    inline const MathGmpRational &operator*=(const MathGmpRational &rhs)
    {
        *this = *this * rhs;
        return *this;
    }
    inline const MathGmpRational &operator/=(const MathGmpRational &rhs)
    {
        *this = *this / rhs;
        return *this;
    }
    inline const MathGmpRational &operator-=(const MathGmpRational &rhs)
    {
        *this = *this - rhs;
        return *this;
    }
    inline const MathGmpRational &operator+=(const MathGmpRational &rhs)
    {
        *this = *this + rhs;
        return *this;
    }
    inline bool operator<(const MathGmpRational &rhs) const
    {
        return mpq_cmp(v, rhs.v) < 0;
    }
    inline bool operator>(const MathGmpRational &rhs) const
    {
        return mpq_cmp(v, rhs.v) > 0;
    }
    inline bool operator<=(const MathGmpRational &rhs) const
    {
        return mpq_cmp(v, rhs.v) <= 0;
    }
    inline bool operator>=(const MathGmpRational &rhs) const
    {
        return mpq_cmp(v, rhs.v) >= 0;
    }
    inline MathGmpRational abs() const
    {
        MathGmpRational result;
        mpq_abs(result.v, v);
        return result;
    }
    inline MathGmpRational sqrt()
    {
        MathGmpRational result;
        mpz_t num, den;
        mpz_inits(num, den, NULL);

        mpq_canonicalize(v);
        mpq_get_num(num, v);
        mpq_get_den(den, v);

        if (mpz_perfect_square_p(num) && mpz_perfect_square_p(den))
        {
            mpz_sqrt(num, num);
            mpz_sqrt(den, den);

            mpq_set_num(result.v, num);
            mpq_set_den(result.v, den);
        }
        else
        {
            mpf_t res_real;
            mpf_init2(res_real, 2048);
            mpf_set_q(res_real, v); // Conversion q -> f
            mpf_sqrt(res_real, res_real);

            mpq_set_f(result.v, res_real);
            mpf_clear(res_real);
        }

        mpz_clears(num, den, NULL);
        return result;
    }
    inline bool isInfinite() const
    {
        return false;
    }
    inline bool isfinite() const
    {
        return !isInfinite();
    }
    inline double to_double() const
    {
        return mpq_get_d(v);
    }
};

inline std::ostream &operator<<(std::ostream &stream, const MathGmpRational &val)
{
    return stream << mpq_get_d(val.v);
}
inline MathGmpRational operator+(const MathGmpRational &lhs, const MathGmpRational &rhs)
{
    MathGmpRational result;
    mpq_add(result.v, lhs.v, rhs.v);
    return result;
}
inline MathGmpRational operator-(const MathGmpRational &lhs, const MathGmpRational &rhs)
{
    MathGmpRational result;
    mpq_sub(result.v, lhs.v, rhs.v);
    return result;
}
inline MathGmpRational operator/(const MathGmpRational &lhs, const MathGmpRational &rhs)
{
    MathGmpRational result;
    mpq_div(result.v, lhs.v, rhs.v);
    return result;
}
inline MathGmpRational operator*(const MathGmpRational &lhs, const MathGmpRational &rhs)
{
    MathGmpRational result;
    mpq_mul(result.v, lhs.v, rhs.v);
    return result;
}
