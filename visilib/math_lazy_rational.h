/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2021 by Denis Haumont

This file is part of Visilib.

Visilib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Visilib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Visilib. If not, see <http://www.gnu.org/licenses/>
*/

#pragma once

#if defined(ENABLE_GMP) && defined(ENABLE_MPFI)

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>

#include <gmp.h>
#include <mpfi.h>

#include "math_gmp_rational.h"

/* A lazy exact rational: the value of an expression of rationals, known by intervals that enclose it, and computed
exactly only when the intervals cannot decide a comparison.

In GMP rational, the exact coordinates of the deep vertices of the CSG grow to thousands of bits, and every exact
operation on them costs. Most decisions are signs far from zero, which an interval tells without the exact value:
- every value keeps the expression it comes from (a node of a DAG: an operation and its operands) and a double
  interval, computed with it: the bounds in round to nearest, widened by one ulp each (round to nearest errs by at most
  half an ulp), so that the interval always contains the exact value;
- a comparison the double intervals cannot decide asks the MPFI intervals (cIntervalPrecision bits), computed on
  demand from those of the operands and kept;
- a comparison those cannot decide either (an exact zero, or a difference beyond their precision) asks the exact GMP
  rational values, computed on demand and kept. A node whose exact value is known drops its operands: the DAG is
  pruned, and its intervals are tightened around the value.
The exact value of a node is then always the one GMP rational would have computed: the answers are those of
MathGmpRational, only the cost differs.

The square root keeps the meaning it has for MathGmpRational: a rational approximation (see MathGmpRational::sqrt),
taken from the exact value of its argument. */

/* The decisions taken at each stage, for the statistics of a query (see MathLazyRational::getStatistics).*/
struct MathLazyRationalStatistics
{
    size_t doubleDecisions = 0;     /**< @brief Comparisons decided by the double intervals*/
    size_t intervalDecisions = 0;   /**< @brief Comparisons decided by the MPFI intervals*/
    size_t exactDecisions = 0;      /**< @brief Comparisons decided by the exact values*/
    size_t exactEvaluations = 0;    /**< @brief Nodes whose exact value has been computed*/
};

struct MathLazyRational
{
    /** @brief The precision, in bits, of the MPFI intervals.*/
    static constexpr mpfr_prec_t cIntervalPrecision = 256;

    MathLazyRational() : mNode(zero())
    {
    }

    MathLazyRational(double d) : mNode(fromDouble(d))
    {
    }

    MathLazyRational(int i) : mNode(fromDouble(double(i)))
    {
    }

    MathLazyRational(const MathGmpRational& q) : mNode(fromRational(q))
    {
    }

    MathLazyRational operator-() const
    {
        return MathLazyRational(makeNode(NEG, mNode, nullptr));
    }

    MathLazyRational& operator+=(const MathLazyRational& rhs) { mNode = makeNode(ADD, mNode, rhs.mNode); return *this; }
    MathLazyRational& operator-=(const MathLazyRational& rhs) { mNode = makeNode(SUB, mNode, rhs.mNode); return *this; }
    MathLazyRational& operator*=(const MathLazyRational& rhs) { mNode = makeNode(MUL, mNode, rhs.mNode); return *this; }
    MathLazyRational& operator/=(const MathLazyRational& rhs) { mNode = makeNode(DIV, mNode, rhs.mNode); return *this; }

    friend MathLazyRational operator+(const MathLazyRational& a, const MathLazyRational& b) { return MathLazyRational(makeNode(ADD, a.mNode, b.mNode)); }
    friend MathLazyRational operator-(const MathLazyRational& a, const MathLazyRational& b) { return MathLazyRational(makeNode(SUB, a.mNode, b.mNode)); }
    friend MathLazyRational operator*(const MathLazyRational& a, const MathLazyRational& b) { return MathLazyRational(makeNode(MUL, a.mNode, b.mNode)); }
    friend MathLazyRational operator/(const MathLazyRational& a, const MathLazyRational& b) { return MathLazyRational(makeNode(DIV, a.mNode, b.mNode)); }

    friend bool operator<(const MathLazyRational& a, const MathLazyRational& b) { return compare(a, b) < 0; }
    friend bool operator>(const MathLazyRational& a, const MathLazyRational& b) { return compare(a, b) > 0; }
    friend bool operator<=(const MathLazyRational& a, const MathLazyRational& b) { return compare(a, b) <= 0; }
    friend bool operator>=(const MathLazyRational& a, const MathLazyRational& b) { return compare(a, b) >= 0; }
    friend bool operator==(const MathLazyRational& a, const MathLazyRational& b) { return compare(a, b) == 0; }
    friend bool operator!=(const MathLazyRational& a, const MathLazyRational& b) { return compare(a, b) != 0; }

    /** @brief -1, 0 or +1 as a is below, equal to or above b: by the double intervals, then the MPFI intervals, then
    the exact values.*/
    static int compare(const MathLazyRational& a, const MathLazyRational& b)
    {
        const Node& x = *a.mNode;
        const Node& y = *b.mNode;
        if (x.hi < y.lo)
        {
            getStatistics().doubleDecisions++;
            return -1;
        }
        if (x.lo > y.hi)
        {
            getStatistics().doubleDecisions++;
            return 1;
        }
        if (x.lo == x.hi && y.lo == y.hi && x.lo == y.lo)
        {
            // two intervals of width zero: each value is that double exactly
            getStatistics().doubleDecisions++;
            return 0;
        }
        const int interval = mpfi_cmp(intervalOf(a.mNode), intervalOf(b.mNode));
        if (interval != 0 && !mpfi_nan_p(intervalOf(a.mNode)) && !mpfi_nan_p(intervalOf(b.mNode)))
        {
            getStatistics().intervalDecisions++;
            return interval < 0 ? -1 : 1;
        }
        getStatistics().exactDecisions++;
        const int exact = mpq_cmp(exactOf(a.mNode).v, exactOf(b.mNode).v);
        return exact < 0 ? -1 : (exact > 0 ? 1 : 0);
    }

    /** @brief The sign of the value: -1, 0 or +1.*/
    int sign() const
    {
        return compare(*this, MathLazyRational());
    }

    MathLazyRational abs() const
    {
        return sign() < 0 ? -*this : *this;
    }

    /** @brief The approximate square root of MathGmpRational, from the exact value.*/
    MathLazyRational sqrt() const
    {
        MathGmpRational value = exactOf(mNode);
        return MathLazyRational(value.sqrt());
    }

    /** @brief The exact value.*/
    MathGmpRational exact() const
    {
        return exactOf(mNode);
    }

    /** @brief A double close to the value: the middle of the double interval, or of the MPFI interval when the double
    one is unbounded, or the exact value rounded when neither is bounded.*/
    double to_double() const
    {
        const Node& n = *mNode;
        if (n.exact)
        {
            return mpq_get_d(n.exact->v);
        }
        if (std::isfinite(n.lo) && std::isfinite(n.hi))
        {
            return n.lo == n.hi ? n.lo : 0.5 * n.lo + 0.5 * n.hi;
        }
        mpfr_t middle;
        mpfr_init2(middle, cIntervalPrecision);
        mpfi_mid(middle, intervalOf(mNode));
        const double result = mpfr_get_d(middle, MPFR_RNDN);
        mpfr_clear(middle);
        if (std::isfinite(result))
        {
            return result;
        }
        return mpq_get_d(exactOf(mNode).v);
    }

    bool isfinite() const
    {
        return true;
    }

    /** @brief The double interval of the value, for the tests.*/
    double lower() const { return mNode->lo; }
    double upper() const { return mNode->hi; }

    /** @brief Whether the exact value is known (computed, or the value of a leaf), for the tests.*/
    bool hasExact() const { return mNode->exact != nullptr; }

    /** @brief Whether the expression of the value still holds operands (it is pruned once the exact value is known).*/
    bool hasOperands() const { return mNode->a != nullptr; }

    /** @brief The MPFI interval of the value, for the tests.*/
    void getInterval(mpfi_t anInterval) const
    {
        mpfi_set(anInterval, intervalOf(mNode));
    }

    static MathLazyRationalStatistics& getStatistics()
    {
        static MathLazyRationalStatistics statistics;
        return statistics;
    }

private:
    enum Operation : unsigned char
    {
        LEAF,
        ADD,
        SUB,
        MUL,
        DIV,
        NEG
    };

    struct Node
    {
        Operation operation = LEAF;
        double lo = 0.0;
        double hi = 0.0;
        std::shared_ptr<Node> a;
        std::shared_ptr<Node> b;
        bool hasInterval = false;
        mpfi_t interval;
        std::unique_ptr<MathGmpRational> exact;

        ~Node()
        {
            if (hasInterval)
            {
                mpfi_clear(interval);
            }
        }
    };

    explicit MathLazyRational(const std::shared_ptr<Node>& aNode) : mNode(aNode)
    {
    }

    static double down(double x)
    {
        return std::nextafter(x, -std::numeric_limits<double>::infinity());
    }

    static double up(double x)
    {
        return std::nextafter(x, std::numeric_limits<double>::infinity());
    }

    static void setUnbounded(Node& n)
    {
        n.lo = -std::numeric_limits<double>::infinity();
        n.hi = std::numeric_limits<double>::infinity();
    }

    /** @brief The interval around a rational: mpq_get_d truncates, with an error below one ulp, and returns an infinity
    beyond the double range.*/
    static void setAround(Node& n, const MathGmpRational& q)
    {
        const double d = mpq_get_d(q.v);
        if (!std::isfinite(d))
        {
            setUnbounded(n);
            return;
        }
        if (isExactlyDouble(q, d))
        {
            n.lo = n.hi = d;
        }
        else
        {
            n.lo = down(d);
            n.hi = up(d);
        }
    }

    /** @brief Whether the rational is exactly the double d.*/
    static bool isExactlyDouble(const MathGmpRational& q, double d)
    {
        if (mpz_cmp_ui(mpq_denref(q.v), 1) != 0 && (mpz_scan1(mpq_denref(q.v), 0) + 1 != mpz_sizeinbase(mpq_denref(q.v), 2)))
        {
            return false; // a denominator that is not a power of two
        }
        mpq_t t;
        mpq_init(t);
        mpq_set_d(t, d);
        const bool equal = mpq_equal(t, q.v) != 0;
        mpq_clear(t);
        return equal;
    }

    static const std::shared_ptr<Node>& zero()
    {
        static const std::shared_ptr<Node> node = fromDouble(0.0);
        return node;
    }

    static std::shared_ptr<Node> fromDouble(double d)
    {
        std::shared_ptr<Node> n = std::make_shared<Node>();
        if (std::isfinite(d))
        {
            n->lo = n->hi = d;
        }
        else
        {
            setUnbounded(*n);
        }
        n->exact.reset(new MathGmpRational(d));
        return n;
    }

    static std::shared_ptr<Node> fromRational(const MathGmpRational& q)
    {
        std::shared_ptr<Node> n = std::make_shared<Node>();
        n->exact.reset(new MathGmpRational(q));
        setAround(*n, q);
        return n;
    }

    /** @brief A node of the expression, with its double interval.*/
    static std::shared_ptr<Node> makeNode(Operation anOperation, const std::shared_ptr<Node>& a, const std::shared_ptr<Node>& b)
    {
        std::shared_ptr<Node> n = std::make_shared<Node>();
        n->operation = anOperation;
        n->a = a;
        n->b = b;

        const double alo = a->lo, ahi = a->hi;
        double lo = 0.0, hi = 0.0;
        switch (anOperation)
        {
        case NEG:
            lo = -ahi;
            hi = -alo;
            break;
        case ADD:
            lo = down(alo + b->lo);
            hi = up(ahi + b->hi);
            break;
        case SUB:
            lo = down(alo - b->hi);
            hi = up(ahi - b->lo);
            break;
        case MUL:
        {
            const double p1 = alo * b->lo, p2 = alo * b->hi, p3 = ahi * b->lo, p4 = ahi * b->hi;
            lo = down(std::min(std::min(p1, p2), std::min(p3, p4)));
            hi = up(std::max(std::max(p1, p2), std::max(p3, p4)));
            break;
        }
        case DIV:
        {
            if (b->lo <= 0.0 && b->hi >= 0.0)
            {
                setUnbounded(*n);
                return n;
            }
            const double q1 = alo / b->lo, q2 = alo / b->hi, q3 = ahi / b->lo, q4 = ahi / b->hi;
            lo = down(std::min(std::min(q1, q2), std::min(q3, q4)));
            hi = up(std::max(std::max(q1, q2), std::max(q3, q4)));
            break;
        }
        default:
            break;
        }
        // an infinity times zero, or an overflow: the interval says nothing, the next stage decides
        if (std::isnan(lo) || std::isnan(hi))
        {
            setUnbounded(*n);
        }
        else
        {
            n->lo = lo;
            n->hi = hi;
        }
        return n;
    }

    /** @brief The MPFI interval of the node, computed with those of its operands and kept.*/
    static mpfi_srcptr intervalOf(const std::shared_ptr<Node>& aNode)
    {
        Node& n = *aNode;
        if (n.hasInterval)
        {
            return n.interval;
        }
        mpfi_init2(n.interval, cIntervalPrecision);
        if (n.exact)
        {
            mpfi_set_q(n.interval, n.exact->v);
        }
        else
        {
            switch (n.operation)
            {
            case NEG: mpfi_neg(n.interval, intervalOf(n.a)); break;
            case ADD: mpfi_add(n.interval, intervalOf(n.a), intervalOf(n.b)); break;
            case SUB: mpfi_sub(n.interval, intervalOf(n.a), intervalOf(n.b)); break;
            case MUL: mpfi_mul(n.interval, intervalOf(n.a), intervalOf(n.b)); break;
            case DIV: mpfi_div(n.interval, intervalOf(n.a), intervalOf(n.b)); break;
            default: break;
            }
        }
        n.hasInterval = true;
        return n.interval;
    }

    /** @brief The exact value of the node, computed with those of its operands and kept: the operands are then dropped
    and the intervals tightened.*/
    static const MathGmpRational& exactOf(const std::shared_ptr<Node>& aNode)
    {
        Node& n = *aNode;
        if (n.exact)
        {
            return *n.exact;
        }
        std::unique_ptr<MathGmpRational> value(new MathGmpRational());
        switch (n.operation)
        {
        case NEG: mpq_neg(value->v, exactOf(n.a).v); break;
        case ADD: mpq_add(value->v, exactOf(n.a).v, exactOf(n.b).v); break;
        case SUB: mpq_sub(value->v, exactOf(n.a).v, exactOf(n.b).v); break;
        case MUL: mpq_mul(value->v, exactOf(n.a).v, exactOf(n.b).v); break;
        case DIV: mpq_div(value->v, exactOf(n.a).v, exactOf(n.b).v); break;
        default: break;
        }
        n.exact = std::move(value);
        getStatistics().exactEvaluations++;

        n.a.reset();
        n.b.reset();
        setAround(n, *n.exact);
        if (n.hasInterval)
        {
            mpfi_set_q(n.interval, n.exact->v);
        }
        return *n.exact;
    }

    std::shared_ptr<Node> mNode;
};

inline std::ostream& operator<<(std::ostream& stream, const MathLazyRational& val)
{
    return stream << val.to_double();
}

#endif
