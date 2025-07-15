#include <gtest/gtest.h>

#include <ginac/ginac.h>

#include "Sym/Expression.h"

using namespace sym;

using T = double;
//using T = GiNaC::numeric;

TEST(TestSym, Symbol)
{
    Expression<T> const e{Symbol{"x"}};
}

TEST(TestSym, Number)
{
    Expression<T> const e{1.2};
}

TEST(TestSym, AddSymbolAndNumber)
{
    Symbol const x{"x"};

    Expression<T> const e1{x + T{1.2}};

    std::cout << e1.simplify().str() << std::endl;
}

TEST(TestSym, AddSymbols)
{
    Symbol const x{"x"};
    Symbol const y{"y"};

    Expression<T> const e{operator+<T>(x, y)};

    std::cout << e.str() << std::endl;
}

TEST(TestSym, 2xAdd3x)
{
    Symbol const x{"x"};

    Expression<T> const e{T{2} * x + T{3} * x};

    std::cout << e.simplify().str() << std::endl;
}

TEST(TestSym, 2xSub2x)
{
    Symbol const x{"x"};

    Expression<T> const e{T{2} * x - T{2} * x};

    std::cout << e.simplify().str() << std::endl;
}

TEST(TestSym, xMulx)
{
    Symbol const x{"x"};

    Expression<T> const e1{operator*<T>(x, x)};

    //std::cout << e1.str() << std::endl;
    std::cout << e1.simplify().str() << std::endl;

    Expression<T> const e2{e1 * x};
    Expression<T> const e3{e1.simplify() * x};

    //std::cout << e2.str() << std::endl;
    std::cout << e2.simplify().str() << std::endl;
    //std::cout << e3.str() << std::endl;
    std::cout << e3.simplify().str() << std::endl;
}

TEST(TestSym, xPowaMulxPowa)
{
    Symbol const x{"x"}, a{"a"};

    Mul<T> const x2{x, a};
    Expression<T> const e{x2 * x2};

    std::cout << e.simplify().str() << std::endl;
}

TEST(TestSym, xPowaMulxPowb)
{
    Symbol const x{"x"}, a{"a"}, b{"b"};

    Mul<T> const xa{x, a};
    Mul<T> const xb{x, b};
    Expression<T> const e{xa * xb};

    std::cout << e.simplify().str() << std::endl;
}

TEST(TestSym, AddEquality)
{
    Symbol const x{"x"}, y{"y"};

    Expression<T> const e1{operator+<T>(x, y)};
    Expression<T> const e2{operator+<T>(y, x)};

    EXPECT_TRUE(e1 == e2);
}

TEST(TestSym, MulEquality)
{
    Symbol const x{"x"}, y{"y"};

    Expression<T> const e1{operator*<T>(x, y)};
    Expression<T> const e2{operator*<T>(y, x)};

    EXPECT_TRUE(e1 == e2);
}

TEST(TestSym, 2xPow2Addx)
{
    Symbol const x{"x"};

    Expression<T> const e{T{2} * x + operator*<T>(x, x)};

    std::cout << e.simplify().str() << std::endl;
}

TEST(TestSym, 1AddxFactor2Addy)
{
    Symbol const x{"x"}, y{"y"};

    Expression<T> const e{(T{1} + x) * (T{2} + y)};

    std::cout << e.simplify().str() << std::endl;
    std::cout << e.expand().simplify().str() << std::endl;
}

TEST(TestSym, 1AddxPow2)
{
    Symbol const x{"x"};

    Expression<T> const e{(T{1} + x) * (T{1} + x)};

    std::cout << e.simplify().str() << std::endl;
    std::cout << e.expand().str() << std::endl;
    std::cout << e.expand().simplify().str() << std::endl;
}
