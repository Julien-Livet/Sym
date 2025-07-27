#ifndef SYM_EXPRESSION_H
#define SYM_EXPRESSION_H

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <type_traits>
#include <vector>

#include <boost/algorithm/string/replace.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace sym
{
    template <typename T>
    class Expression;

    template <typename T>
    class Add;

    template <typename T>
    class Mul;

    template <typename T>
    class Composition;

    struct Symbol;
}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
sym::Expression<T> operator+(sym::Symbol const& lhs, T const& rhs);
template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
sym::Expression<T> operator+(T const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator+(sym::Symbol const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator+(sym::Expression<T> const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator+(sym::Symbol const& lhs, sym::Expression<T> const& rhs);
template <typename T>
sym::Expression<T> operator+(sym::Expression<T> const& lhs, T const& rhs);
template <typename T>
sym::Expression<T> operator+(T const& lhs, sym::Expression<T> const& rhs);
template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
sym::Expression<T> operator-(sym::Symbol const& lhs, T const& rhs);
template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
sym::Expression<T> operator-(T const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator-(sym::Symbol const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator-(sym::Expression<T> const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator-(sym::Symbol const& lhs, sym::Expression<T> const& rhs);
template <typename T>
sym::Expression<T> operator-(sym::Expression<T> const& lhs, T const& rhs);
template <typename T>
sym::Expression<T> operator-(T const& lhs, sym::Expression<T> const& rhs);
template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
sym::Expression<T> operator*(sym::Symbol const& lhs, T const& rhs);
template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
sym::Expression<T> operator*(T const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator*(sym::Symbol const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator*(sym::Expression<T> const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator*(sym::Symbol const& lhs, sym::Expression<T> const& rhs);
template <typename T>
sym::Expression<T> operator*(sym::Expression<T> const& lhs, T const& rhs);
template <typename T>
sym::Expression<T> operator*(T const& lhs, sym::Expression<T> const& rhs);
template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
sym::Expression<T> operator/(sym::Symbol const& lhs, T const& rhs);
template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
sym::Expression<T> operator/(T const& lhs, sym::Symbol const& rsh);
template <typename T>
sym::Expression<T> operator/(sym::Symbol const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator/(sym::Expression<T> const& lhs, sym::Symbol const& rhs);
template <typename T>
sym::Expression<T> operator/(sym::Symbol const& lhs, sym::Expression<T> const& rhs);
template <typename T>
sym::Expression<T> operator/(sym::Expression<T> const& lhs, T const& rhs);
template <typename T>
sym::Expression<T> operator/(T const& lhs, sym::Expression<T> const& rhs);
template <typename T>
sym::Mul<T> operator-(sym::Expression<T> const& expression);
template <typename T>
sym::Mul<T> operator-(sym::Add<T> const& add);
template <typename T>
sym::Mul<T> operator-(sym::Mul<T> const& mul);
template <typename T>
sym::Mul<T> operator-(sym::Composition<T> const& composition);
template <typename T>
sym::Add<T> operator+(sym::Mul<T> const& mul1, sym::Mul<T> const& mul2);
template <typename T>
sym::Mul<T> operator*(sym::Add<T> const& add1, sym::Add<T> const& add2);

namespace sym
{
    class Function
    {
        public:
            std::string name;

            Function(std::string const& name, size_t args = 1) : name{name}, args_{args}
            {
            }

            size_t args() const
            {
                return args_;
            }

            bool operator==(Function const& other) const = default;
            bool operator!=(Function const& other) const = default;

        private:
            size_t args_;
    };

    template <typename T>
    class Mul
    {
        public:
            Mul()
            {
            }

            Mul(Expression<T> const& expression, Expression<T> power = T{1})
            {
                if (expression == T{1} || (expression.isMul() && expression.mul().expressions().empty()))
                    power = T{1};

                if (expression.isNull())
                {
                    expressions_.emplace_back(T{0});
                    powers_.emplace_back(T{1});
                }
                else if (power != T{0})
                {
                    if ((power == T{1} || (power.isMul() && power.mul().expressions().empty()))
                        && expression.isMul())
                    {
                        expressions_ = expression.mul().expressions_;
                        powers_ = expression.mul().powers_;
                    }
                    else
                    {
                        expressions_.emplace_back(expression);
                        powers_.emplace_back(power);
                    }
                }
            }

            Mul& operator=(Expression<T> const& expression)
            {
                expressions_.clear();
                powers_.clear();

                if (expression.isNull())
                {
                    expressions_.emplace_back(T{0});
                    powers_.emplace_back(T{1});
                }
                else if (expression.isMul())
                {
                    expressions_ = expression.mul().expressions_;
                    powers_ = expression.mul().powers_;
                }
                else
                {
                    expressions_.emplace_back(expression);
                    powers_.emplace_back(T{1});
                }

                return *this;
            }

            Mul<T>& inverse()
            {
                for (auto& p: powers_)
                    p = -p;
            }

            Mul<T> inverse() const
            {
                auto m{*this};

                return m.inverse();
            }

            bool operator==(Mul const& other) const
            {
                size_t count{0};

                for (size_t i{0}; i < other.expressions_.size(); ++i)
                {
                    auto const& e{other.expressions_[i]};

                    auto const it{std::find(expressions_.begin(),
                                            expressions_.end(),
                                            e)};

                    if (it == expressions_.end())
                        return false;

                    auto const j{std::distance(expressions_.begin(), it)};

                    if (powers_[j] == other.powers_[i])
                        ++count;
                }

                return count == expressions_.size();
            }

            bool operator!=(Mul const& other) const
            {
                return !operator==(other);
            }

            bool isNumber() const
            {
                if (expressions_.empty() || isNull())
                    return true;

                for (size_t i{0}; i < expressions_.size(); ++i)
                {
                    if (!expressions_[i].isNumber() || !powers_[i].isNumber())
                        return false;
                }

                return true;
            }

            T number() const
            {
                assert(isNumber());
                
                if (isNull())
                    return T{0};

                T x{1};

                for (size_t i{0}; i < expressions_.size(); ++i)
                    x *= std::pow(expressions_[i].number(), powers_[i].number());

                return x;
            }

            T numberFactor() const
            {
                T x{1};

                for (size_t i{0}; i < expressions_.size(); ++i)
                {
                    if (powers_[i].isNumber())
                    {
                        if (expressions_[i].isNumber())
                            x *= std::pow(expressions_[i].number(), powers_[i].number());
                        else if (expressions_[i].isMul())
                            x *= std::pow(expressions_[i].mul().numberFactor(), powers_[i].number());
                    }
                }

                return x;
            }

            Mul<T> simplify() const
            {
                Mul<T> m;

                if (numberFactor() != T{1})
                {
                    m.expressions_.emplace_back(numberFactor());
                    m.powers_.emplace_back(T{1});
                }

                for (size_t i{0}; i < expressions_.size(); ++i)
                {
                    if (!expressions_[i].isNumber() || !powers_[i].isNumber())
                        m *= Mul<T>{expressions_[i].simplify(), powers_[i].simplify()};
                }

                return m;
            }

            bool isNull(T eps = 1e-4) const
            {
                if (std::abs(numberFactor()) <= eps)
                    return true;

                for (auto const& e: expressions_)
                {
                    if (e.isNull(eps))
                        return true;
                }

                return false;
            }

            Expression<T> operator+(Expression<T> const& other) const
            {
                auto tmp{other};

                return tmp += other;
            }

            Mul<T> operator*(Mul<T> const& other) const
            {
                auto tmp{*this};

                return tmp *= other;
            }

            Mul<T> operator/(Mul<T> const& other) const
            {
                auto tmp{*this};

                return tmp /= other;
            }

            Mul<T>& operator*=(Expression<T> const& other)
            {
                if (isNull())
                    return *this;

                if (other.isNull())
                {
                    expressions_.clear();
                    powers_.clear();
                    expressions_.emplace_back(T{0});
                    powers_.emplace_back(T{1});

                    return *this;
                }

                if (other.isMul())
                    return *this *= other.mul();

                auto const it{std::find(expressions_.begin(),
                                        expressions_.end(),
                                        other)};

                if (it != expressions_.end())
                {
                    auto const i{std::distance(expressions_.begin(), it)};

                    powers_[i] += T{1};
                }
                else
                {
                    expressions_.emplace_back(other);
                    powers_.emplace_back(T{1});
                }

                return *this;
            }

            Mul<T>& operator*=(Mul<T> const& other)
            {
                for (size_t i{0}; i < other.expressions_.size(); ++i)
                {
                    auto const& e{other.expressions_[i]};
                    auto const it{std::find(expressions_.begin(),
                                            expressions_.end(),
                                            e)};

                    if (it != expressions_.end())
                    {
                        auto const j{std::distance(expressions_.begin(), it)};

                        powers_[j] += other.powers_[i];
                    }
                    else
                    {
                        expressions_.emplace_back(e);
                        powers_.emplace_back(other.powers_[i]);
                    }
                }

                return *this;
            }

            Mul<T>& operator/=(Expression<T> const& other)
            {
                auto const it{std::find(expressions_.begin(),
                                        expressions_.end(),
                                        other)};

                if (it != expressions_.end())
                {
                    auto const i{std::distance(expressions_.begin(), it)};

                    powers_[i] -= 1;
                }
                else
                {
                    expressions_.emplace_back(other);
                    powers_.emplace_back(-1);
                }

                return *this;
            }

            Mul<T>& operator/=(Mul<T> const& other)
            {
                for (size_t i{0}; i < other.expressions_.size(); ++i)
                {
                    auto const& e{other.expressions_[i]};
                    auto const it{expressions_.find(e)};

                    if (it != expressions_.end())
                    {
                        auto const j{std::distance(expressions_.begin(), it)};

                        powers_[j] -= other.powers_[i];
                    }
                    else
                    {
                        expressions_.emplace_back(e);
                        powers_.emplace_back(-other.powers_[i]);
                    }
                }

                return *this;
            }

            auto const& expressions() const
            {
                return expressions_;
            }

            auto const& powers() const
            {
                return powers_;
            }

            std::string str() const
            {
                if (expressions_.empty())
                    return "1";

                if (isNull())
                    return "0";

                std::string s;

                if (expressions_.front().isNumber())
                {
                    std::ostringstream oss;

                    oss << expressions_.front().number();

                    s += oss.str();
                }
                else if (powers_.front().isComplex() && (expressions_.front().isMul()
                                                         || expressions_.front().isAdd())
                         || (expressions_.front().isAdd() && (expressions_.size() > 1
                                                              || (powers_.front().isNumber()
                                                                  && powers_.front().number() != T{1}))))
                    s += "(" + expressions_.front().str() + ")";
                else
                    s += expressions_.front().str();

                if (!(powers_.front().isNumber() && powers_.front().number() == 1))
                {
                    s += "**";

                    if (powers_.front().isNumber())
                    {
                        std::ostringstream oss;

                        oss << powers_.front().number();

                        s += oss.str();
                    }
                    else if (powers_.front().isAdd() || powers_.front().isMul())
                        s += "(" + powers_.front().str() + ")";
                    else
                        s += powers_.front().str();
                }

                for (size_t i{1}; i < expressions_.size(); ++i)
                {
                    s += "*";

                    if (expressions_[i].isNumber())
                    {
                        std::ostringstream oss;

                        oss << expressions_[i].number();

                        s += oss.str();
                    }
                    else if (powers_[i].isComplex() && (expressions_[i].isMul()
                                                        || expressions_[i].isAdd())
                            || expressions_[i].isAdd())
                        s += "(" + expressions_[i].str() + ")";
                    else
                        s += expressions_[i].str();

                    if (!(powers_[i].isNumber() && powers_[i].number() == 1))
                    {
                        s += "**";

                        if (powers_[i].isNumber())
                        {
                            std::ostringstream oss;

                            oss << powers_[i].number();

                            s += oss.str();
                        }
                        else if (powers_[i].isAdd() || powers_[i].isMul())
                            s += "(" + powers_[i].str() + ")";
                        else
                            s += powers_[i].str();
                    }
                }

                return s;
            }

            Mul<T> factor() const
            {
                Mul<T> e{*this};

                for (auto& expression: e.expressions_)
                    expression = expression.factor();

                return e;
            }

            Expression<T> expand() const
            {
                if (expressions_.empty())
                    return Add<T>{};

                {
                    Expression<T> const* term{nullptr};
                    Expression<T> const* power{nullptr};

                    for (size_t i{0}; i < expressions_.size(); ++i)
                    {
                        if (powers_[i].isNumber()
                            && powers_[i].number() == int(powers_[i].number())
                            && std::abs(powers_[i].number()) > 1
                            && expressions_[i].isAdd()
                            && expressions_[i].add().expressions().size() > 1)
                        {
                            term = &expressions_[i];
                            power = &powers_[i];
                            break;
                        }
                    }

                    if (term && power)
                    {
                        Expression<T> const term1{term->add().expressions().front()};
                        Expression<T> term2;

                        for (size_t i{1}; i < term->add().expressions().size(); ++i)
                            term2 += term->add().expressions()[i];

                        auto const n{static_cast<size_t>(std::abs(power->number()))};

                        Add<T> sum;

                        for (size_t i{0}; i <= n; ++i)
                            sum.append(Mul<T>{static_cast<T>(boost::math::binomial_coefficient<double>(n, i))}
                                       * Mul<T>{term1, static_cast<T>(i)}
                                       * Mul<T>{term2, static_cast<T>(n - i)});

                        if (power->number() < 0)
                            return Mul<T>{sum, -1}.expand();
                        else
                            return sum.expand();
                    }
                }

                auto const& first{expressions_.front()};
                Expression<T> const* second{nullptr};
                size_t index{0};

                for (size_t i{1}; i < expressions_.size(); ++i)
                {
                    if (expressions_[i].isAdd()
                        && expressions_[i].add().expressions().size() > 1
                        && powers_[i].isNumber()
                        && std::abs(powers_[i].number()) == 1)
                    {
                        second = &expressions_[i];
                        index = i;
                        break;
                    }
                }

                if (!second)
                    return *this;

                Mul<T> term;

                for (size_t i{1}; i < expressions_.size(); ++i)
                {
                    if (index != i)
                        term *= Mul<T>(expressions_[i], powers_[i]);
                }

                Add<T> sum;

                if (first.isAdd())
                {
                    for (size_t i{0}; i < second->add().expressions().size(); ++i)
                    {
                        for (size_t j{0}; j < first.add().expressions().size(); ++j)
                            sum.append(first.add().expressions()[j] * second->add().expressions()[i]);
                    }
                }
                else
                {
                    for (size_t i{0}; i < second->add().expressions().size(); ++i)
                        sum.append(first * second->add().expressions()[i]);
                }

                if (term.expressions_.empty())
                    return sum;
                else
                    return sum * term;
            }
            
            T eval(std::map<std::string, T> const& map) const
            {
                T value{1};

                for (size_t i{0}; i < expressions_.size(); ++i)
                    value *= std::pow(expressions_[i].eval(map), powers_[i].eval(map));

                return value;
            }

        private:
            std::vector<Expression<T> > expressions_;
            std::vector<Expression<T> > powers_;
    };

    template <typename T>
    class Add
    {
        public:
            Add()
            {
            }

            Add(Expression<T> const& expression)
            {
                if (expression.isAdd())
                    expressions_ = expression.add().expressions_;
                else if (!expression.isNull())
                    expressions_.emplace_back(expression);
            }

            void append(Expression<T> const& expression)
            {
                if (!expression.isNull())
                    expressions_.emplace_back(expression);
            }

            Add& operator=(Expression<T> const& expression)
            {
                expressions_.clear();

                if (expression.isAdd())
                    expressions_ = expression.add().expressions_;
                else if (!expression.isNull())
                    expressions_.emplace_back(expression);

                return *this;
            }

            bool operator==(Add const& other) const
            {
                size_t count{0};

                for (size_t i{0}; i < other.expressions_.size(); ++i)
                {
                    auto const& e{other.expressions_[i]};

                    auto const it{std::find(expressions_.begin(),
                                            expressions_.end(),
                                            e)};

                    if (it == expressions_.end())
                        return false;

                    ++count;
                }

                return count == expressions_.size();
            }

            bool operator!=(Add const& other) const
            {
                return !operator==(other);
            }

            auto const& expressions() const
            {
                return expressions_;
            }

            bool isNumber() const
            {
                if (expressions_.empty() || isNull())
                    return true;

                for (auto const& e: expressions_)
                {
                    if (!e.isNumber())
                        return false;
                }

                return true;
            }

            T number() const
            {
                assert(isNumber());

                if (isNull())
                    return T{0};

                T x{0};

                for (auto const& e: expressions_)
                    x += e.number();

                return x;
            }
            
            T additiveNumber() const
            {
                T x{0};

                for (auto const& e: expressions_)
                    if (e.isNumber())
                        x += e.number();

                return x;
            }

            Add<T> factor() const
            {
                Add<T> e;

                for (auto& expression: expressions_)
                    e += expression.factor();

                return e;
            }

            Add<T> simplify() const
            {
                if (expressions_.size() == 1)
                {
                    auto const& mul{expressions_.front()};

                    if (mul.expressions().size() == 1
                        && mul.powers().front().isNumber()
                        && mul.powers().front().number() == T{1}
                        && mul.powers().front().isAdd())
                        return mul.expressions().front().add();
                }

                Add<T> a;

                if (isNumber())
                    a.expressions_.emplace_back(number());
                else
                {
                    for (auto const& e: expressions_)
                        a.expressions_.emplace_back(e.simplify());
                }

                return a;
            }

            bool isNull(T eps = 1e-4) const
            {
                if (expressions_.empty())
                    return true;

                for (auto const& e: expressions_)
                {
                    if (!e.isNull(eps))
                        return false;
                }

                return true;
            }

            Expression<T> operator*(Expression<T> const& other) const
            {
                auto tmp{other};

                return tmp *= other;
            }

            Add<T> operator+(Add<T> const& other) const
            {
                auto d{*this};

                return d += other;
            }

            Add<T> operator-(Add<T> const& other) const
            {
                auto d{*this};

                return d -= other;
            }

            Add<T>& operator+=(Expression<T> const& other)
            {
                if (other.isAdd())
                    return *this += other.add();

                if (other.isNull())
                    return *this;

                if (other.isMul())
                {
                    auto const& mul{other.mul()};

                    std::vector<size_t> occurences(expressions_.size());

                    for (size_t i{0}; i < expressions_.size(); ++i)
                    {
                        for (auto const& e: mul.expressions())
                        {
                            //TODO: add test on same powers too?
                            if (std::find(expressions_[i].expressions().begin(),
                                          expressions_[i].expressions().end(),
                                          e) != expressions_[i].expressions().end())
                                ++occurences[i];
                        }
                    }

                    auto const it{std::max_element(occurences.begin(), occurences.end())};

                    if (it != occurences.end() && *it)
                    {
                        auto const i{std::distance(occurences.begin(), it)};

                        std::vector<Expression<T> > factors1;
                        std::vector<Expression<T> > terms;

                        auto const mulExpressions{mul.expressions()};
                        auto const expressions{expressions_[i].expressions()};

                        for (size_t j{0}; j < mulExpressions.size(); ++j)
                        {
                            auto const& e{mulExpressions[j]};
                            bool isTerm{false};
                            auto const itBis{std::find(expressions.begin(),
                                                       expressions.end(),
                                                       e)};

                            if (itBis != expressions.end())
                            {
                                auto const k{std::distance(expressions.begin(), itBis)};

                                isTerm = (expressions_[i].powers()[k] == mul.powers()[j]);
                            }

                            Mul<T> const expr{mul.expressions()[j], mul.powers()[j]};

                            if (isTerm)
                                terms.emplace_back(expr);
                            else
                                factors1.emplace_back(expr);
                        }

                        std::vector<Expression<T> > factors2;

                        for (size_t j{0}; j < expressions.size(); ++j)
                        {
                            Mul<T> const expr{expressions[j], expressions_[i].powers()[j]};

                            if (std::find(terms.begin(),
                                          terms.end(),
                                          expr) == terms.end())
                                factors2.emplace_back(expr);
                        }

                        Mul<T> factor1, factor2;

                        for (auto const& f: factors1)
                            factor1 *= f;

                        for (auto const& f: factors2)
                            factor2 *= f;

                        Add<T> sum{factor1};
                        sum.append(factor2);

                        Mul<T> factor{sum};

                        for (auto const& f: terms)
                            factor *= f;

                        expressions_[i] = factor;
                    }
                    else
                        expressions_.emplace_back(other);

                    return *this;
                }

                std::vector<size_t> occurences(expressions_.size());

                for (size_t i{0}; i < expressions_.size(); ++i)
                {
                    //TODO: add test on same powers too?
                    if (std::find(expressions_[i].expressions().begin(),
                                  expressions_[i].expressions().end(),
                                  other) != expressions_[i].expressions().end())
                        ++occurences[i];
                }

                auto const it{std::max_element(occurences.begin(), occurences.end())};

                if (it != occurences.end() && *it)
                {
                    auto const i{std::distance(occurences.begin(), it)};

                    std::vector<Expression<T> > factors;
                    Expression<T> const term{other};

                    for (auto const& e: expressions_[i].expressions())
                    {
                        if (e != term)
                            factors.emplace_back(e);
                    }

                    Mul<T> factor;

                    for (auto const& f: factors)
                        factor *= f;

                    if (factor.expressions().empty())
                    {
                        factor = T{2};

                        if (term != T{1})
                            expressions_[i] = factor * term;
                        else
                            expressions_[i] = factor;
                    }
                    else
                        expressions_[i] = (factor + T{1}) * term;
                }
                else
                    expressions_.emplace_back(other);

                return *this;
            }

            Add<T>& operator+=(Add<T> const& other)
            {
                for (auto const& e: other.expressions_)
                    *this += e;

                return *this;
            }

            Add<T>& operator-=(Expression<T> const& other)
            {
                return *this += -other;
            }

            Add<T>& operator-=(Add<T> const& other)
            {
                for (auto const& e: other.expressions_)
                    *this -= e;

                return *this;
            }

            std::string str() const
            {
                if (expressions_.empty())
                    return "0";

                std::string s;

                s += expressions_.front().str();

                for (size_t i{1}; i < expressions_.size(); ++i)
                {
                    s += "+";
                    s += expressions_[i].str();
                }

                return s;
            }

            Add<T> expand() const
            {
                Add<T> e;

                for (auto const& expression: expressions_)
                    e.append(expression.expand());

                return e;
            }
            
            T eval(std::map<std::string, T> const& map) const
            {
                T value{0};

                for (auto const& e: expressions_)
                    value += e.eval(map);

                return value;
            }

        private:
            std::vector<Mul<T> > expressions_;
    };

    template <typename T>
    class Composition
    {
        public:
            Composition(Function const& function,
                        std::vector<Expression<T> > const& expressions,
                        std::function<T(std::vector<T>)> const& eval,
                        std::function<bool()> const& isNumber = [] () {return false;},
                        std::function<T()> const& number = [] () {return T{0};})
                : function_{function}, expressions_{expressions}, eval_{eval}, isNumber_{isNumber}, number_{number}
            {
                assert(function.args() == expressions.size());
            }

            auto const& function() const
            {
                return function_;
            }

            auto const& expressions() const
            {
                return expressions_;
            }

            std::string str() const
            {
                std::string s;

                s += function_.name + "(";

                if (expressions_.size())
                {
                    s += expressions_.front().str();

                    for (size_t i{1}; i < expressions_.size(); ++i)
                    {
                        s += ", ";
                        s += expressions_[i].str();
                    }
                }

                s += ")";

                return s;
            }

            bool operator==(Composition<T> const& other) const
            {
                return function_ == other.function_ && expressions_ == other.expressions_;
            }

            bool operator!=(Composition<T> const& other) const
            {
                return !operator==(other);
            }

            Composition<T> expand() const
            {
                auto e{*this};

                for (auto& expression: e.expressions_)
                    expression = expression.expand();

                return e;
            }

            Composition<T> simplify() const
            {
                auto e{*this};

                for (auto& expression: e.expressions_)
                    expression = expression.simplify();

                return e;
            }


            Composition<T> factor() const
            {
                auto e{*this};

                for (auto& expression: e.expressions_)
                    expression = expression.factor();

                return e;
            }

            bool isNumber() const
            {
                return isNumber_();
            }
            
            T number() const
            {
                assert(isNumber());
                
                return number_();
            }
            
            T eval(std::map<std::string, T> const& map) const
            {
                std::vector<T> params;
                params.reserve(expressions_.size());

                for (auto const& e: expressions_)
                    params.emplace_back(e.eval(map));

                return eval_(params);
            }

        private:
            Function function_;
            std::vector<Expression<T> > expressions_;
            std::function<T(std::vector<T>)> eval_;
            std::function<bool()> isNumber_;
            std::function<T()> number_;
    };

    struct Symbol
    {
        std::string name;

        bool operator==(Symbol const& other) const = default;
        bool operator!=(Symbol const& other) const = default;
    };

    template <typename T>
    class Expression
    {
        public:
            enum Type
            {
                NumberType,
                SymbolType,
                AddType,
                MulType,
                CompositionType
            };

            Expression() : type_{AddType}, add_{std::make_unique<Add<T> >()}
            {
            }

            Expression(T const& number) : type_{NumberType}, number_{std::make_unique<T>(number)}
            {
            }

            Expression(Symbol const& symbol) : type_{SymbolType}, symbol_{std::make_unique<Symbol>(symbol)}
            {
            }

            Expression(Add<T> const& add) : type_{AddType}, add_{std::make_unique<Add<T> >(add)}
            {
            }

            Expression(Mul<T> const& mul) : type_{MulType}, mul_{std::make_unique<Mul<T> >(mul)}
            {
            }

            Expression(Composition<T> const& composition) : type_{CompositionType}, composition_{std::make_unique<Composition<T> >(composition)}
            {
            }

            Type type() const
            {
                return type_;
            }

            Expression(Expression const& other) : type_{other.type_}
            {
                if (type_ == NumberType)
                    number_ = std::make_unique<T>(*other.number_);
                else if (type_ == SymbolType)
                    symbol_ = std::make_unique<Symbol>(*other.symbol_);
                else if (type_ == AddType)
                    add_ = std::make_unique<Add<T> >(*other.add_);
                else if (other.isMul())
                {
                    mul_ = std::make_unique<Mul<T> >(other.mul());
                    type_ = MulType;
                }
                else// if (type_ == CompositionType)
                    composition_ = std::make_unique<Composition<T> >(*other.composition_);
            }

            Expression& operator=(Expression const& other)
            {
                type_ = other.type_;
                number_.reset();
                symbol_.reset();
                add_.reset();
                composition_.reset();

                if (type_ == NumberType)
                    number_ = std::make_unique<T>(*other.number_);
                else if (type_ == SymbolType)
                    symbol_ = std::make_unique<Symbol>(*other.symbol_);
                else if (type_ == AddType)
                    add_ = std::make_unique<Add<T> >(*other.add_);
                else if (other.isMul())
                {
                    mul_ = std::make_unique<Mul<T> >(other.mul());
                    type_ = MulType;
                }
                else// if (type_ == CompositionType)
                    composition_ = std::make_unique<Composition<T> >(*other.composition_);

                return *this;
            }

            bool operator==(Expression const& other) const
            {
                if (type_ != other.type_)
                    return false;

                if (type_ == NumberType)
                    return *number_ == *other.number_;
                else if (type_ == SymbolType)
                    return *symbol_ == *other.symbol_;
                else if (type_ == AddType)
                    return *add_ == *other.add_;
                else if (type_ == MulType)
                    return *mul_ == *other.mul_;
                else// if (type_ == CompositionType)
                    return *composition_ == *other.composition_;
            }

            bool operator!=(Expression const& other) const
            {
                return !operator==(other);
            }

            bool isNull(T eps = 1e-4) const
            {
                if (type_ == NumberType)
                    return std::abs(*number_) <= eps;
                else if (type_ == SymbolType)
                    return false;
                else if (type_ == AddType)
                    return add_->isNull(eps);
                else if (type_ == MulType)
                    return mul_->isNull(eps);
                else// if (type_ == CompositionType)
                {
                    if (composition_->isNumber())
                        return std::abs(composition_->number()) <= eps;

                    return false;
                }
            }

            bool isMul() const
            {
                if (isAdd()
                    && add().expressions().size() == 1)
                    return true;

                return type_ == MulType;
            }

            bool isAdd() const
            {
                return type_ == AddType;
            }

            bool isComposition() const
            {
                return type_ == CompositionType;
            }

            bool isSymbol() const
            {
                return type_ == SymbolType;
            }

            bool isNumber() const
            {
                if (type_ == NumberType)
                    return true;
                else if (type_ == SymbolType)
                    return false;
                else if (type_ == AddType)
                    return add_->isNumber();
                else if (type_ == MulType)
                    return mul_->isNumber();
                else// if (type_ == CompositionType)
                    return composition_->isNumber();
            }

            T number() const
            {
                assert(isNumber());

                T x{0};

                if (type_ == NumberType)
                    x = *number_;
                else if (type_ == AddType)
                    x = add_->number();
                else if (type_ == MulType)
                    x = mul_->number();
                else// if (type_ == CompositionType)
                    x = composition_->number();

                return x;
            }

            auto const& symbol() const
            {
                assert(type_ == SymbolType);

                return *symbol_;
            }

            auto const& add() const
            {
                assert(type_ == AddType);

                return *add_;
            }

            auto const& mul() const
            {
                if (isAdd()
                    && add().expressions().size() == 1)
                    return add().expressions().front();

                assert(type_ == MulType);

                return *mul_;
            }

            auto const& composition() const
            {
                assert(type_ == CompositionType);

                return *composition_;
            }

            Expression<T> operator+(Expression<T> const& other) const
            {
                auto e{*this};

                return e += other;
            }

            Expression<T> operator-(Expression<T> const& other) const
            {
                auto e{*this};

                return e -= other;
            }

            Expression<T> operator*(Expression<T> const& other) const
            {
                auto e{*this};

                return e *= other;
            }

            Expression<T> operator/(Expression<T> const& other) const
            {
                auto e{*this};

                return e /= other;
            }

            Expression<T>& operator+=(Expression<T> const& other)
            {
                if (isAdd())
                {
                    *add_ += other;

                    return *this;
                }
                else if (other.isAdd())
                {
                    auto e{other.add()};

                    e += *this;

                    *this = e;

                    return *this;
                }

                Add<T> add{*this};
                add += other;

                *this = add;

                return *this;
            }

            Expression<T>& operator-=(Expression<T> const& other)
            {
                return *this += -other;
            }

            Expression<T>& operator*=(Expression<T> const& other)
            {
                if (isMul())
                {
                    *this = mul() * other;

                    return *this;
                }
                else if (other.isMul())
                {
                    auto e{other.mul()};

                    e *= *this;

                    *this = e;

                    return *this;
                }

                Mul<T> mul{*this};
                mul *= other;

                *this = mul;

                return *this;
            }

            Expression<T>& operator/=(Expression<T> const& other)
            {
                if (isMul())
                {
                    *mul_ *= Mul<T>{other, T{-1}};

                    return *this;
                }
                else if (other.isMul())
                {
                    auto e{other.mul()};
                    e.inverse();
                    e *= *this;

                    *this = e;

                    return *this;
                }

                Mul<T> mul{*this};
                mul *= Mul<T>{other, T{-1}};

                *this = mul;

                return *this;
            }

            Expression<T> expand() const
            {
                auto e{*this};

                if (e.isMul())
                    e = e.mul().expand();
                else if (e.isAdd())
                    e = e.add().expand();
                else if (e.isComposition())
                    e = e.composition().expand();

                return e;
            }

            Expression<T> factor() const
            {
                auto e{*this};

                if (e.isMul())
                    e.mul_ = std::make_unique<Mul<T> >(e.mul().factor());
                else if (e.isAdd())
                    e.add_ = std::make_unique<Add<T> >(e.add_->factor());
                else if (e.isComposition())
                    e.composition_ = std::make_unique<Composition<T> >(e.composition_->factor());

                return e;
            }

            Expression<T> simplify() const
            {
                auto e{*this};

                if (e.isMul())
                {
                    e.mul_ = std::make_unique<Mul<T> >(mul().simplify());
                    e.type_ = MulType;
                    e.add_.reset();
                }
                else if (e.isAdd())
                    e.add_ = std::make_unique<Add<T> >(e.add_->simplify());
                else if (e.isComposition())
                    e.composition_ = std::make_unique<Composition<T> >(e.composition_->simplify());

                return e;
            }

            std::string str() const
            {
                std::string s;

                if (type_ == NumberType)
                {
                    std::ostringstream oss;

                    oss << *number_;

                    s += oss.str();
                }
                else if (type_ == SymbolType)
                    s += symbol_->name;
                else if (type_ == AddType)
                    s += add_->str();
                else if (isMul())
                    s += mul().str();
                else// if (type_ == CompositionType)
                    s += composition_->str();

                boost::replace_all(s, "+-", "-");

                return s;
            }

            bool isComplex() const
            {
                if (isNumber() || isComposition() || isSymbol())
                    return false;
                else
                    return true;
            }

            T eval(std::map<std::string, T> const& map) const
            {
                if (isSymbol())
                    return map.at(symbol_->name);
                else if (isNumber())
                    return number();
                else if (isAdd())
                    return add_->eval(map);
                else if (isMul())
                    return mul().eval(map);
                else// if (isComposition())
                    return composition_->eval(map);
            }

        private:
            Type type_;
            std::unique_ptr<T> number_;
            std::unique_ptr<Symbol> symbol_;
            std::unique_ptr<Add<T> > add_;
            std::unique_ptr<Mul<T> > mul_;
            std::unique_ptr<Composition<T> > composition_;
    };

    template <typename T>
    Expression<T> abs(Expression<T> const& e)
    {
        return Composition<T>(Function("abs"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::abs(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::abs(e.number());});
    }

    template <typename T>
    Expression<T> ceil(Expression<T> const& e)
    {
        return Composition<T>(Function("ceil"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::ceil(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::ceil(e.number());});
    }

    template <typename T>
    Expression<T> floor(Expression<T> const& e)
    {
        return Composition<T>(Function("floor"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::floor(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::floor(e.number());});
    }

    template <typename T>
    Expression<T> min(Expression<T> const& x, Expression<T> const& y)
    {
        return Composition<T>(Function("min", 2), std::vector<Expression<T> >{x, y},
                              [] (std::vector<T> const& params) {return std::min(params[0], params[1]);},
                              [x, y] () {return x.isNumber() && y.isNumber();},
                              [x, y] () {return std::min(x.number(), y.number());});
    }

    template <typename T>
    Expression<T> max(Expression<T> const& x, Expression<T> const& y)
    {
        return Composition<T>(Function("max", 2), std::vector<Expression<T> >{x, y},
                              [] (std::vector<T> const& params) {return std::max(params[0], params[1]);},
                              [x, y] () {return x.isNumber() && y.isNumber();},
                              [x, y] () {return std::max(x.number(), y.number());});
    }

    template <typename T>
    Expression<T> exp(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "log")
            return e.composition().expressions().front();

        return Composition<T>(Function("exp"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::exp(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::exp(e.number());});
    }

    template <typename T>
    Expression<T> log(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "exp")
            return e.composition().expressions().front();

        return Composition<T>(Function("log"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::log(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::log(e.number());});
    }

    template <typename T>
    Expression<T> sin(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "asin")
            return e.composition().expressions().front();

        return Composition<T>(Function("sin"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::sin(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::sin(e.number());});
    }

    template <typename T>
    Expression<T> asin(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "sin")
            return e.composition().expressions().front();

        return Composition<T>(Function("asin"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::asin(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::asin(e.number());});
    }

    template <typename T>
    Expression<T> cos(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "acos")
            return e.composition().expressions().front();

        return Composition<T>(Function("cos"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::cos(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::cos(e.number());});
    }

    template <typename T>
    Expression<T> acos(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "cos")
            return e.composition().expressions().front();

        return Composition<T>(Function("acos"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::acos(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::acos(e.number());});
    }

    template <typename T>
    Expression<T> tan(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "atan")
            return e.composition().expressions().front();

        return Composition<T>(Function("tan"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::tan(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::tan(e.number());});
    }

    template <typename T>
    Expression<T> atan(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "tan")
            return e.composition().expressions().front();

        return Composition<T>(Function("atan"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::atan(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::atan(e.number());});
    }

    template <typename T>
    Expression<T> sinh(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "asinh")
            return e.composition().expressions().front();

        return Composition<T>(Function("sinh"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::sinh(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::sinh(e.number());});
    }

    template <typename T>
    Expression<T> asinh(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "sinh")
            return e.composition().expressions().front();

        return Composition<T>(Function("asinh"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::asinh(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::asinh(e.number());});
    }

    template <typename T>
    Expression<T> cosh(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "acosh")
            return e.composition().expressions().front();

        return Composition<T>(Function("cosh"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::cosh(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::cosh(e.number());});
    }

    template <typename T>
    Expression<T> acosh(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "cosh")
            return e.composition().expressions().front();

        return Composition<T>(Function("acosh"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::acosh(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::acosh(e.number());});
    }

    template <typename T>
    Expression<T> tanh(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "atanh")
            return e.composition().expressions().front();

        return Composition<T>(Function("tanh"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::tanh(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::tanh(e.number());});
    }

    template <typename T>
    Expression<T> atanh(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "tanh")
            return e.composition().expressions().front();

        return Composition<T>(Function("atanh"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::atanh(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::atanh(e.number());});
    }

    template <typename T>
    Expression<T> cot(Expression<T> const& e)
    {
        return Composition<T>(Function("cot"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return T{1} / std::tan(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return T{1} / std::tan(e.number());});
    }

    template <typename T>
    Expression<T> sqrt(Expression<T> const& e)
    {
        return Composition<T>(Function("sqrt"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return std::sqrt(params.front());},
                              [e] () {return e.isNumber();},
                              [e] () {return std::sqrt(e.number());});
    }

    template <typename T>
    Expression<T> inverse(Expression<T> const& e)
    {
        if (e.isComposition() && e.composition().function().name == "inverse")
            return e.composition().expressions().front();

        return Composition<T>(Function("inverse"), std::vector<Expression<T> >{e},
                              [] (std::vector<T> const& params) {return T{1} / params.front();},
                              [e] () {return e.isNumber();},
                              [e] () {return T{1} / e.number();});
    }
}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >*>
sym::Expression<T> operator+(sym::Symbol const& lhs, T const& rhs)
{
    sym::Expression<T> e{lhs};

    return e += rhs;
}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >*>
sym::Expression<T> operator+(T const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e += rhs;
}

template <typename T>
sym::Expression<T> operator+(sym::Symbol const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e += rhs;
}

template <typename T>
sym::Expression<T> operator+(sym::Expression<T> const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e += rhs;
}

template <typename T>
sym::Expression<T> operator+(sym::Symbol const& lhs, sym::Expression<T> const& rhs)
{
    sym::Expression<T> e{lhs};

    return e += rhs;
}

template <typename T>
sym::Expression<T> operator+(sym::Expression<T> const& lhs, T const& rhs)
{
    sym::Expression<T> e{lhs};

    return e += rhs;
}

template <typename T>
sym::Expression<T> operator+(T const& lhs, sym::Expression<T> const& rhs)
{
    sym::Expression<T> e{lhs};

    return e += rhs;
}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >*>
sym::Expression<T> operator-(sym::Symbol const& lhs, T const& rhs)
{
    sym::Expression<T> e{lhs};

    return e -= rhs;
}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >*>
sym::Expression<T> operator-(T const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e -= rhs;
}

template <typename T>
sym::Expression<T> operator-(sym::Symbol const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e -= rhs;
}

template <typename T>
sym::Expression<T> operator-(sym::Expression<T> const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e -= rhs;
}

template <typename T>
sym::Expression<T> operator-(sym::Symbol const& lhs, sym::Expression<T> const& rhs)
{
    sym::Expression<T> e{lhs};

    return e -= rhs;
}

template <typename T>
sym::Expression<T> operator-(sym::Expression<T> const& lhs, T const& rhs)
{
    sym::Expression<T> e{lhs};

    return e -= rhs;
}

template <typename T>
sym::Expression<T> operator-(T const& lhs, sym::Expression<T> const& rhs)
{
    sym::Expression<T> e{lhs};

    return e -= rhs;
}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >*>
sym::Expression<T> operator*(sym::Symbol const& lhs, T const& rhs)
{
    sym::Expression<T> e{lhs};

    return e *= rhs;
}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >*>
sym::Expression<T> operator*(T const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e *= rhs;
}

template <typename T>
sym::Expression<T> operator*(sym::Symbol const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e *= rhs;
}

template <typename T>
sym::Expression<T> operator*(sym::Expression<T> const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e *= rhs;
}

template <typename T>
sym::Expression<T> operator*(sym::Symbol const& lhs, sym::Expression<T> const& rhs)
{
    sym::Expression<T> e{lhs};

    return e *= rhs;
}

template <typename T>
sym::Expression<T> operator*(sym::Expression<T> const& lhs, T const& rhs)
{
    sym::Expression<T> e{lhs};

    return e *= rhs;
}

template <typename T>
sym::Expression<T> operator*(T const& lhs, sym::Expression<T> const& rhs)
{
    sym::Expression<T> e{lhs};

    return e *= rhs;
}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >*>
sym::Expression<T> operator/(sym::Symbol const& lhs, T const& rhs)
{
    sym::Expression<T> e{lhs};

    return e /= rhs;
}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T> >*>
sym::Expression<T> operator/(T const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e /= rhs;
}

template <typename T>
sym::Expression<T> operator/(sym::Symbol const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e /= rhs;
}

template <typename T>
sym::Expression<T> operator/(sym::Expression<T> const& lhs, sym::Symbol const& rhs)
{
    sym::Expression<T> e{lhs};

    return e /= rhs;
}

template <typename T>
sym::Expression<T> operator/(sym::Symbol const& lhs, sym::Expression<T> const& rhs)
{
    sym::Expression<T> e{lhs};

    return e /= rhs;
}

template <typename T>
sym::Expression<T> operator/(sym::Expression<T> const& lhs, T const& rhs)
{
    sym::Expression<T> e{lhs};

    return e /= rhs;
}

template <typename T>
sym::Expression<T> operator/(T const& lhs, sym::Expression<T> const& rhs)
{
    sym::Expression<T> e{lhs};

    return e /= rhs;
}

template <typename T>
sym::Mul<T> operator-(sym::Expression<T> const& expression)
{
    sym::Mul<T> e{expression};

    e *= T{-1};

    return e;
}

template <typename T>
sym::Mul<T> operator-(sym::Add<T> const& add)
{
    sym::Mul<T> e{add};

    e *= T{-1};

    return e;
}

template <typename T>
sym::Mul<T> operator-(sym::Mul<T> const& mul)
{
    sym::Mul<T> e{mul};

    e *= T{-1};

    return e;
}

template <typename T>
sym::Mul<T> operator-(sym::Composition<T> const& composition)
{
    sym::Mul<T> e{composition};

    e *= T{-1};

    return e;
}

template <typename T>
sym::Add<T> operator+(sym::Mul<T> const& mul1, sym::Mul<T> const& mul2)
{
    sym::Add<T> e{mul1};

    return e += mul2;
}

template <typename T>
sym::Mul<T> operator*(sym::Add<T> const& add1, sym::Add<T> const& add2)
{
    sym::Mul<T> e{add1};

    return e *= add2;
}

#endif // SYM_EXPRESSION_H
