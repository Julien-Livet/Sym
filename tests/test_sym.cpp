#include <gtest/gtest.h>

#include "Sym/Expression.h"

using namespace sym;

using T = double;

TEST(TestSym, Symbol)
{
    Expression<T> const e{Symbol{"x"}};

    EXPECT_TRUE(e.str() == "x");
}

TEST(TestSym, Number)
{
    Expression<T> const e{1.2};

    EXPECT_TRUE(e.str() == "1.2");
}

TEST(TestSym, AddSymbolAndNumber)
{
    Symbol const x{"x"};

    Expression<T> const e1{x + T{1.2}};
    Expression<T> const e2{x - T{0.8}};

    EXPECT_TRUE(e1.str() == "x+1.2");
    EXPECT_TRUE(e2.simplify().str() == "x-0.8");
}

TEST(TestSym, AddSymbols)
{
    Symbol const x{"x"};
    Symbol const y{"y"};

    Expression<T> const e{operator+<T>(x, y)};

    EXPECT_TRUE(e.str() == "x+y");
}

TEST(TestSym, 2xAdd3x)
{
    Symbol const x{"x"};

    Expression<T> const e{T{2} * x + T{3} * x};

    EXPECT_TRUE(e.str() == "5*x");
    EXPECT_TRUE(e.simplify().str() == "5*x");
}

TEST(TestSym, 2xSub2x)
{
    Symbol const x{"x"};

    Expression<T> const e{T{2} * x - T{2} * x};

    EXPECT_TRUE(e.str() == "0");
    EXPECT_TRUE(e.simplify().str() == "0");
}

TEST(TestSym, xMulx)
{
    Symbol const x{"x"};

    Expression<T> const e1{operator*<T>(x, x)};

    EXPECT_TRUE(e1.str() == "x**2");

    Expression<T> const e2{e1 * x};
    Expression<T> const e3{e1.simplify() * x};

    EXPECT_TRUE(e2.str() == "x**3");
    EXPECT_TRUE(e3.simplify().str() == "x**3");
}

TEST(TestSym, xPowaMulxPowa)
{
    Symbol const x{"x"}, a{"a"};

    Mul<T> const x2{x, a};
    Expression<T> const e{x2 * x2};

    EXPECT_TRUE(e.str() == "x**(2*a)");
    EXPECT_TRUE(e.simplify().str() == "x**(2*a)");
}

TEST(TestSym, xPowaMulxPowb)
{
    Symbol const x{"x"}, a{"a"}, b{"b"};

    Mul<T> const xa{x, a};
    Mul<T> const xb{x, b};
    Expression<T> const e{xa * xb};

    EXPECT_TRUE(e.str() == "x**(a+b)");
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

TEST(TestSym, xPow2Add2x)
{
    Symbol const x{"x"};

    Expression<T> const e{T{2} * x + operator*<T>(x, x)};

    EXPECT_TRUE(e.str() == "2*x+x**2" || e.str() == "x**2+2*x");
}

TEST(TestSym, 1AddxFactor2Addy)
{
    Symbol const x{"x"}, y{"y"};

    Expression<T> const e{(T{1} + x) * (T{2} + y)};

    EXPECT_TRUE(e.simplify().str() == "(1+x)*(2+y)");
    EXPECT_TRUE(e.expand().simplify().str() == "2+2*x+y+x*y");
}

TEST(TestSym, 1AddxPow2)
{
    Symbol const x{"x"};

    Expression<T> const e{(T{1} + x) * (T{1} + x)};

    EXPECT_TRUE(e.str() == "(1+x)**2");
    EXPECT_TRUE(e.simplify().str() == "(1+x)**2");
    EXPECT_TRUE(e.expand().simplify().str() == "x**2+2*x+1");
}

TEST(TestSym, exp)
{
    EXPECT_TRUE(sym::exp(Expression<T>{Symbol{"x"}}).str() == "exp(x)");
}

TEST(TestSym, expLog)
{
    EXPECT_TRUE(sym::exp(sym::log(Expression<T>{Symbol{"x"}})).str() == "x");
}

TEST(TestSym, 2xMul3x)
{
    Symbol const x{"x"};

    Expression<T> const e{T{2} * x * T{3} * x};

    EXPECT_TRUE(e.simplify().str() == "6*x**2");
}

TEST(TestSym, ComplexNullExpression)
{
    Symbol const y{"y"};

    Expression<T> const e = (0.996642*(((1.38482e-10*(operator*<T>(y, y)))+(0.00191499*(y)))*((1.38472e-10*(operator*<T>(y, y)))+(0.00191499*(y)))))+((-3.65487e-06*(operator*<T>(y, y)))+(-9.6662e-13*(y)))+4.66616e-12;

    EXPECT_TRUE(e.isNull());
}

TEST(TestSym, sin0AddLog1AddCos0)
{
    Expression<T> const e{sym::sin(Expression(0.0)) + sym::log(Expression(1.0)) + sym::cos(Expression(0.0))};

    EXPECT_TRUE(e.simplify().str() == "1");
}

TEST(TestSym, EvalMinLogExp)
{
    Symbol const x{"x"};

    Expression<T> const e{sym::min(sym::log(Expression<T>(x)), sym::exp(Expression<T>(x)))};

    std::map<std::string, T> map;
    map["x"] = 2.0;

    EXPECT_TRUE(e.eval(map) == std::min(std::log(map["x"]), std::exp(map["x"])));
}
