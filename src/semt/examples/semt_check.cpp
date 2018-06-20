/*
 * semt_example.cpp
 *
 * Created on: Nov 19, 2009
 * Author: stefan
 */

#include <numeric>
#include <iostream>
#include <algorithm>

#include "semt/Semt.h"
#include "semt/Common.h"

using namespace SEMT;
using namespace std;

vector<string> & split(const string &s, char delim, vector<string> & elems)
{
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim))
    {
        if (item.length() > 0)
            elems.push_back(item);
    }
    return elems;
}

void test(const string& msg, const string& got, const string& expect)
{
    vector<string> got_lines, expect_lines;
    split(got, '\n', got_lines);
    split(expect, '\n', expect_lines);

    if (got_lines.size() != expect_lines.size())
    {
        cout << "\t\tError: " << msg << ": unequal length\n";
    }
    else
    {
        for (int i = 0; i < got_lines.size(); ++i)
        {
            std::string& gotl = got_lines[i], expectl = expect_lines[i];
            gotl.erase(remove_if(gotl.begin(), gotl.end(), ::isspace), gotl.end());
            expectl.erase(remove_if(expectl.begin(), expectl.end(), ::isspace), expectl.end());

            if (gotl != expectl)
            {
                cout << "\t\tError: " << msg << " in line " << i <<
                        "\ngot:\n" << gotl << ",\n expected:\n" << expectl << "\n\n";
            }
        }
    }
}

void test(const string& msg, double got, double expect)
{
    if (fabs(got - expect) > SEMT_EPS)
    {
        cout << "\t\tError: " << msg << ", got: "
                << got << ", expected: " << expect << ", diff = " << fabs(got - expect) << "\n";
    }
}

void test_exception(const VectorExpr& f, CAR x, const string& msg)
{
    try
    {
        f(x);
        cout << "\t\tError, expected exception: " << msg << endl;
    }
    catch (const runtime_error&)
    {
    }
}

template<class T>
string stream(const T& t)
{
    ostringstream os;
    os << t;
    return os.str();
}

template<typename e>
void semt_test_f(SEMT::Expr<e> ex,
                 const string& f_str,
                 const string& df_str,
                 SEMT::CAR x_0,
                 SEMT_PRECISION f_x0,
                 int line)
{
    const size_t vars = SEMT::Expr<e> ::LastVar + 1;
    SEMT::DifferentiableVectorExpr<vars, 1> F;
    F.push_back(ex);
    ostringstream msg;
    msg << "(line " << line << ") " << ex.toString();
    test(msg.str(), stream(F.get_function()), "[\t" + f_str + "\t]");
    test(msg.str(), stream(F.get_derivative(1)), "[\t" + df_str + "\t]");
    test(msg.str(), F(x_0)[0], f_x0);
}

#define TEST_F(name, expr, f_str, df_str, f_x0) \
    auto name = expr;                           \
    semt_test_f(expr, f_str, df_str, x_0, f_x0, (int)__LINE__);

struct alpha
{
    constexpr static SEMT_PRECISION value = 5.234;
};
Expr<Literal<alpha> > a;
/*
 * You can hack your way around SEMT and use this to return some global variable (you have been warned).
 */
inline SEMT_PRECISION alphaf()
{
    return 1.13423;
}
Expr<Func<alphaf> > af;

Expr<Fold_t<Plus, IntIterator<Variable, 0, 1> >::Result> SumVars;
Expr<Fold_t<Times, IntIterator<Integer, 1, 10> >::Result> Fac10;

#define DISABLE_ELABORATE_TEST 0

#if not DISABLE_ELABORATE_TEST
void test_diffable_expr()
{
    cout << "\tTesting high level methods...\n";
    DifferentiableVectorExpr<2, 1> F;
    F.push_back(x0 + x1);
    F.push_back(x0 - x1);
    F.push_back(x0 * x1);
    F.push_back(x0 / (x1 + One));
    F.push_back(pow(x0, x1));

    F.push_back(sin(x0 + x1));
    F.push_back(cos(x0 + x1));
    F.push_back(cos(sin(x0 * x1)) + a * sin(cos(x0 + x1)));
    F.push_back(pow(cos(x0 * x1), a * sin(x0 + x1)));
    F.push_back(tan(x0 * x1 + ThreeQuarters));

    F.push_back(sinh(x0 * x1 + a));
    F.push_back(cosh(x0 / (x1 + af)));
    F.push_back(exp(x0 * x1) + ln(ThreeQuarters));
    F.push_back(ln(x0 / (x1 + One)));
    F.push_back(tanh(pow(x0, x1)));

    F.push_back(Fac10 * SumVars);

    ostringstream astr, afstr;
    astr << alpha::value;
    afstr << alphaf();

    string Fstr =
        "[(x0 + x1),\n"
        "(x0 - x1),\n"
        "(x0 * x1),\n"
        "(x0 / (x1 + 1)),\n"
        "(x0)^(x1),\n"
        "sin((x0 + x1)),\n"
        "cos((x0 + x1)),\n"
        "(cos(sin((x0 * x1))) + (" + astr.str() + " * sin(cos((x0 + x1))))),\n"
        "(cos((x0 * x1)))^((" + astr.str() + " * sin((x0 + x1)))),\n"
        "tan(((x0 * x1) + 3/4)),\n"
        "sinh(((x0 * x1) + " + astr.str() + ")),\n"
        "cosh((x0 / (x1 + " + afstr.str() + "))),\n"
        "(exp((x0 * x1)) + ln(3/4)),\n"
        "ln((x0 / (x1 + 1))),\n"
        "tanh((x0)^(x1)),\n"
        "(3628800 * (x0 + x1))]";

    string DFstr =
        "[1,\n"
        "1,\n"
        "1,\n"
        "-1,\n"
        "x1,\n"
        "x0,\n"
        "((x1 + 1) / ((x1 + 1))^(2)),\n"
        "((0 - x0) / ((x1 + 1))^(2)),\n"
        "((x0)^(x1) * (x1 / x0)),\n"
        "{ ((x0)^(x1) * ln(x0)), if x0 != 0 },\n"
        "cos((x0 + x1)),\n"
        "cos((x0 + x1)),\n"
        "(-1 * sin((x0 + x1))),\n"
        "(-1 * sin((x0 + x1))),\n"
        "((-1 * (sin(sin((x0 * x1))) * (cos((x0 * x1)) * x1))) + (" + astr.str() + " * (cos(cos((x0 + x1))) * (-1 * sin((x0 + x1)))))),\n"
        "((-1 * (sin(sin((x0 * x1))) * (cos((x0 * x1)) * x0))) + (" + astr.str() + " * (cos(cos((x0 + x1))) * (-1 * sin((x0 + x1)))))),\n"
        "((cos((x0 * x1)))^((" + astr.str() + " * sin((x0 + x1)))) * (((" + astr.str() + " * cos((x0 + x1))) * ln(cos((x0 * x1)))) + (((" + astr.str() + " * sin((x0 + x1))) * (-1 * (sin((x0 * x1)) * x1))) / cos((x0 * x1))))),\n"
        "((cos((x0 * x1)))^((" + astr.str() + " * sin((x0 + x1)))) * (((" + astr.str() + " * cos((x0 + x1))) * ln(cos((x0 * x1)))) + (((" + astr.str() + " * sin((x0 + x1))) * (-1 * (sin((x0 * x1)) * x0))) / cos((x0 * x1))))),\n"
        "(x1 * (1 + (tan(((x0 * x1) + 3/4)))^(2))),\n"
        "(x0 * (1 + (tan(((x0 * x1) + 3/4)))^(2))),\n"
        "(cosh(((x0 * x1) + " + astr.str() + ")) * x1),\n"
        "(cosh(((x0 * x1) + " + astr.str() + ")) * x0),\n"
        "(sinh((x0 / (x1 + " + afstr.str() + "))) * ((x1 + " + afstr.str() + ") / ((x1 + " + afstr.str() + "))^(2))),\n"
        "(sinh((x0 / (x1 + " + afstr.str() + "))) * ((0 - x0) / ((x1 + " + afstr.str() + "))^(2))),\n"
        "(exp((x0 * x1)) * x1),\n"
        "(exp((x0 * x1)) * x0),\n"
        "(((x1 + 1) / ((x1 + 1))^(2)) / (x0 / (x1 + 1))),\n"
        "(((0 - x0) / ((x1 + 1))^(2)) / (x0 / (x1 + 1))),\n"
        "(((x0)^(x1) * (x1 / x0)) * (1 - (tanh((x0)^(x1)))^(2))),\n"
        "{ (((x0)^(x1) * ln(x0)) * (1 - (tanh((x0)^(x1)))^(2))), if x0 != 0 },\n"
        "3628800,\n"
        "3628800]";

    test("DifferentiableVectorExpr function's string", stream(F.get_function()), Fstr);
    test("DifferentiableVectorExpr derivative's string", stream(F.get_derivative(1)), DFstr);
}

void test_simplifications()
{
    cout << "\tTesting simplifications...\n";
    ostringstream os1, os2;
    os1 << Zero + Zero - Zero - Zero + Zero + One + Zero + Zero * Zero * One * Zero
            + Zero + x0 + Zero + Zero + Zero + Zero + Zero + One + Zero + Zero
            + Zero + Zero + Zero;
    os2 << (One * One * Ten * pow(One, Ten) * One / (One + Zero)) + (Ten - Ten) + (x0 * One / x0)
            - One;
    test("simplification", os1.str(), string("((1 + x0) + 1)"));
    test("simplification", os2.str(), string("{10, if x0 != 0}"));
}

void test_partial()
{
    Array x_0(3);
    x_0[0] = 0.386434;
    x_0[1] = 3.2342342;
    x_0[2] = -4.0;

    ostringstream astr, afstr;
    astr << setw(SEMT_STRLEN) << alpha::value;
    afstr << setw(SEMT_STRLEN) << alphaf();

    string f1str =
            "((((x0)^(2) + (x1)^(2)) - (x2)^(2)) / (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + "
                    + astr.str() + ") + " + afstr.str() + "))";
    string df1str =
            "((((2 * x0) * (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + " + astr.str() + ") + "
                    + afstr.str()
                    + ")) - ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) * ((((x0 + x1) + x0) * x1) * x2))) / ((((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + "
                    + astr.str() + ") + " + afstr.str() + "))^(2)),\n"
                            "\t((((2 * x1) * (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + "
                    + astr.str() + ") + " + afstr.str()
                    + ")) - ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) * (((x0 * x1) + (x0 * (x0 + x1))) * x2))) / ((((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + "
                    + astr.str() + ") + " + afstr.str() + "))^(2)),\n"
                            "\t((((0 - (2 * x2)) * (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + "
                    + astr.str() + ") + " + afstr.str()
                    + ")) - ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) * (x2 + (((x0 * (x0 + x1)) * x1) + x2)))) / ((((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + "
                    + astr.str() + ") + " + afstr.str() + "))^(2))";
    cout << "\tTesting initialization, simplification and constant expressions...\n";

    TEST_F(f1,
            (pow(x0, Two) + x1 * x1 - x2 * x2) / ((x0 * (x0 + x1) * x1 + x2) * x2 + pow(One, Ten) + a + af),
            f1str,
            df1str,
            ((((x_0[0]*x_0[0]) + (x_0[1] * x_0[1])) - (x_0[2] * x_0[2])) / (((((((x_0[0] * (x_0[0] + x_0[1])) * x_0[1]) + x_0[2]) * x_0[2]) + 1) + alpha::value) + alphaf())));

    TEST_F(p1,
            One * One *One * (Zero + Zero + Zero ) * One * (Zero + Zero + Zero +One + Zero + Zero) + Ten + ThreeQuarters,
            "(10 + 3/4)",
            "0",
            10.75);
    TEST_F(p2,
            pow(One, Zero) * Zero * Zero * One * One + x0 + x1,
            "(x0 + x1)",
            "1,\n\t1",
            (x_0[0]+x_0[1]));
    TEST_F(p3, x0 + a, "(x0 + " + astr.str() + ")", "1", (x_0[0]+alpha::value));
    TEST_F(p4, x0 + af, "(x0 + " + afstr.str() + ")", "1", (x_0[0]+alphaf()));
    TEST_F(p5, x0 + ThreeQuarters, "(x0 + 3/4)", "1", (x_0[0]+3.0/4));
    TEST_F(p6, a + x0, "(" + astr.str() + " + x0)", "1", (x_0[0]+alpha::value));
    TEST_F(p7, af + x0, "(" + afstr.str() + " + x0)", "1", (x_0[0]+alphaf()));
    TEST_F(p8, ThreeQuarters + x0, "(3/4 + x0)", "1", (x_0[0]+3.0/4));
    TEST_F(p9, Zero + Zero + Zero + Zero + Zero + One + Zero + Zero + Zero + Zero
            + Zero + x0 + Zero + Zero + Zero + Zero + Zero + One + Zero + Zero
            + Zero + Zero + Zero, "((1 + x0) + 1)", "1", x_0[0]+2);

    cout << "\tTesting power...\n";
    TEST_F(epow1, pow(x0, x0), "(x0)^(x0)", "{((x0)^(x0) * (ln(x0) + 1)), if x0 != 0}", pow(x_0[0], x_0[0]));
    TEST_F(epow2, pow(x0, x1), "(x0)^(x1)",
            "((x0)^(x1) * (x1 / x0)),\n"
            "{((x0)^(x1) * ln(x0)), if x0 != 0}",
            pow(x_0[0], x_0[1]));

    TEST_F(pow1, pow(ThreeQuarters, Ten), "(3/4)^(10)", "0", pow(0.75,10));
    TEST_F(pow2, pow(a, Ten), "(" + astr.str() + ")^(10)", "0", pow(alpha::value,10));
    TEST_F(pow3, pow(af, Ten), "(" + afstr.str() + ")^(10)", "0", pow(alphaf(),10));
    //todo: somehow, the next one produces a slight error
    TEST_F(pow4, pow(x0, Ten), "(x0)^(10)", "(10 * (x0)^(9))", pow(x_0[0],10.0));
    TEST_F(pow5, pow(p1, Ten), "("+p1.toString()+")^(10)", "0", pow(p1.apply(x_0),10));
    TEST_F(pow6, pow(p2, Ten), "("+p2.toString()+")^(10)", "(10 * ((x0 + x1))^(9)),\n"
    "\t(10 * ((x0 + x1))^(9))", pow(p2.apply(x_0),10));
    TEST_F(pow7,
            pow(p3, Ten),
            "("+p3.toString()+")^(10)",
            "(10 * ((x0 + " + astr.str() + "))^(9))",
            pow(p3.apply(x_0),10));
    TEST_F(pow8,
            pow(p4, Ten),
            "("+p4.toString()+")^(10)",
            "(10 * ((x0 + " + afstr.str() + "))^(9))",
            pow(p4.apply(x_0),10));
    TEST_F(pow9,
            pow(p5, Ten),
            "("+p5.toString()+")^(10)",
            "(10 * ((x0 + 3/4))^(9))",
            pow(p5.apply(x_0),10));
    TEST_F(newpower, pow(x0 + x1, Tousand), "((x0 + x1))^(1000)", "(1000 * ((x0 + x1))^(999)),\n"
    "\t(1000 * ((x0 + x1))^(999))", pow(SEMT_PRECISION(x_0[0]+x_0[1]),1000));

    TEST_F(compound1, (x2 + (af + ThreeQuarters) * (pow(x0 + x1 + x2, Ten))),
            "(x2 + ((" + afstr.str() + " + 3/4) * (((x0 + x1) + x2))^(10)))",
            "((" + afstr.str() + " + 3/4) * (10 * (((x0 + x1) + x2))^(9))),\n"
            "\t((" + afstr.str() + " + 3/4) * (10 * (((x0 + x1) + x2))^(9))),\n"
            "\t(1 + ((" + afstr.str() + " + 3/4) * (10 * (((x0 + x1) + x2))^(9))))",
            x_0[2] + (alphaf() + 0.75) * pow(x_0[0] + x_0[1] + x_0[2],10));
    TEST_F(compound2,
            (x2 + (af + ThreeQuarters) * pow(x0 + x1 + x2, Ten))/(x2 + (af + ThreeQuarters) * pow(x0 + x1 + x2, Ten)),
            "{1, if (x2 + ((1.13423 + 3/4) * (((x0 + x1) + x2))^(10))) != 0}",
            "{0, if (x2 + ((1.13423 + 3/4) * (((x0 + x1) + x2))^(10))) != 0}",
            1.0);

    cout << "\tTesting trigoniometric functions...\n";
    TEST_F(trig1, sin(x0), "sin(x0)", "cos(x0)", sin(x_0[0]));
    TEST_F(trig2, cos(x0), "cos(x0)", "(-1 * sin(x0))", cos(x_0[0]));
    TEST_F(trig3, tan(x0), "tan(x0)", "(1 + (tan(x0))^(2))", tan(x_0[0]));

    TEST_F(trig4, sinh(x0), "sinh(x0)", "cosh(x0)", sinh(x_0[0]));
    TEST_F(trig5, cosh(x0), "cosh(x0)", "sinh(x0)", cosh(x_0[0]));
    TEST_F(trig6, tanh(x0), "tanh(x0)", "(1 - (tanh(x0))^(2))", tanh(x_0[0]));

    cout << "\tTesting inverse trigoniometric functions...\n";
    TEST_F(itrig1, asin(x0), "arcsin(x0)", "(1 / ((1 - (x0)^(2)))^(1/2))", asin(x_0[0]));
    TEST_F(itrig2, acos(x0), "arccos(x0)", "(-1 / ((1 - (x0)^(2)))^(1/2))", acos(x_0[0]));
    TEST_F(itrig3, atan(x0), "arctan(x0)", "(1 / (1 + (x0)^(2)))", atan(x_0[0]));

    TEST_F(itrig4, asinh(x0), "arsinh(x0)", "(1 / ((1 + (x0)^(2)))^(1/2))", asinh(x_0[0]));
    TEST_F(itrig5, acosh(x1), "arcosh(x1)",
            "0,\n"
            "\t(1 / (((x1)^(2) - 1))^(1/2))", acosh(x_0[1]));
    TEST_F(itrig6, atanh(x0), "artanh(x0)", "(1 / (1 - (x0)^(2)))", atanh(x_0[0]));

    cout << "\tTesting exponential & logarithm...\n";
    TEST_F(exp1, exp(x0), "exp(x0)", "exp(x0)", exp(x_0[0]));
    TEST_F(ln1, ln(x0), "ln(x0)", "(1 / x0)", log(x_0[0]));

    string ifdabsf1 =
        "{ ((((2 * x0) * (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423)) - ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) * ((((x0 + x1) + x0) * x1) * x2))) / ((((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423))^(2)), if ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) / (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423)) > 0 },\n"
        "{ ((((2 * x1) * (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423)) - ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) * (((x0 * x1) + (x0 * (x0 + x1))) * x2))) / ((((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423))^(2)), if ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) / (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423)) > 0 },\n"
        "{ ((((0 - (2 * x2)) * (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423)) - ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) * (x2 + (((x0 * (x0 + x1)) * x1) + x2)))) / ((((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423))^(2)), if ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) / (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423)) > 0 }";

    try
    {
        TEST_F(expln, exp(ln(f1)), "{ "+f1.toString()+", if > 0 }", ifdabsf1, f1.apply(x_0));
        cout << "missing exception for exp(ln) at \n";
    } catch (const runtime_error& e)
    {
    }

    TEST_F(lnexp, ln(exp(f1)), f1str, df1str, f1.apply(x_0));

    TEST_F(intln,
            x0*ln(x0)-x0,
            "((x0 * ln(x0)) - x0)",
            "((ln(x0) + 1) - 1)",
            x_0[0]*log(x_0[0]) - x_0[0]);

    cout << "\tTesting absolute value...\n";

    string dabsf1 = "({ sgn(((((x0)^(2) + (x1)^(2)) - (x2)^(2)) / (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423))), if != 0 } * ((((2 * x0) * (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423)) - ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) * ((((x0 + x1) + x0) * x1) * x2))) / ((((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423))^(2))),\n"
    "({ sgn(((((x0)^(2) + (x1)^(2)) - (x2)^(2)) / (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423))), if != 0 } * ((((2 * x1) * (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423)) - ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) * (((x0 * x1) + (x0 * (x0 + x1))) * x2))) / ((((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423))^(2))),\n"
    "({ sgn(((((x0)^(2) + (x1)^(2)) - (x2)^(2)) / (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423))), if != 0 } * ((((0 - (2 * x2)) * (((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423)) - ((((x0)^(2) + (x1)^(2)) - (x2)^(2)) * (x2 + (((x0 * (x0 + x1)) * x1) + x2)))) / ((((((((x0 * (x0 + x1)) * x1) + x2) * x2) + 1) + 5.234) + 1.13423))^(2)))";
    TEST_F(abs1, abs(f1), string("abs("+f1str+")"), dabsf1, fabs(f1.apply(x_0)));

    TEST_F(abs2, abs(Zero), "abs(0)", "0", 0.0);
    TEST_F(abs3, abs(x0 - x1), "abs((x0 - x1))",
            "{ sgn((x0 - x1)), if != 0 },\n"
            "({ sgn((x0 - x1)), if != 0 } * -1)", fabs(x_0[0] - x_0[1]));
    TEST_F(abs4, abs(x0), "abs(x0)", "{ sgn(x0), if != 0 }", fabs(x_0[0]));

    DifferentiableVectorExpr<1, 1> abs_deriv;
    abs_deriv.push_back(abs4);
    VectorExpr dabs4 = abs_deriv.get_derivative(1);
    x_0[0] = 2.213;
    test("Error with absolute value:", dabs4(x_0).at(0), 1.0);
    x_0[0] = -2.213;
    test("Error with absolute value:", dabs4(x_0).at(0), -1.0);

    x_0[0] = -1.0;

    try
    {
        diff_at(abs(x0 + One), x0, x_0);
        cout << "\t\tDifferentiating " << abs(x0 + One)
                << " for x0 at -1 should have thrown an exception.\n";
    } catch (const runtime_error& msg)
    {
    }
}

void test_conditionals()
{
    cout << "\tTesting conditionals...\n";

    Array x(3, 1.0);
    DifferentiableVectorExpr<3, 2> DVE;

    auto f1 = abs(x0 + x1) + Ten;
    auto f2 = sgn(x0 + x1 + One) + x2;
    auto f3 = exp(ln((x0 + x1 + Two))) + x2;
    DVE.push_back(f1);
    DVE.push_back(f2);
    DVE.push_back(f3);

    VectorExpr F = DVE.get_function();
    VectorExpr dF = DVE.get_derivative(1);
    VectorExpr ddF = DVE.get_derivative(2);

    dF(x);
    x[0] = -1.0;
    test_exception(dF, x, "differentiate abs at zero.");
    x[0] = -2.0;
    test_exception(dF, x, "differentiate sgn at zero.");
    x[0] = -3.0;
    test_exception(F, x, "evaluate exp(ln(0)).");

    string f_str =
        "[(abs((x0 + x1)) + 10),\n"
        "(sgn(((x0 + x1) + 1)) + x2),\n"
        "({ ((x0 + x1) + 2), if > 0 } + x2)]";

    string df_str =
        "[{ sgn((x0 + x1)), if != 0 },\n"
        "{ sgn((x0 + x1)), if != 0 },\n"
        "0,\n"
        "{ 0, if ((x0 + x1) + 1) != 0 },\n"
        "{ 0, if ((x0 + x1) + 1) != 0 },\n"
        "1,\n"
        "{ 1, if ((x0 + x1) + 2) > 0 },\n"
        "{ 1, if ((x0 + x1) + 2) > 0 },\n"
        "1]";

    string ddf_str =
        "[{ { 0, if (x0 + x1) != 0 }, if sgn((x0 + x1)) != 0 },\n"
        "{ { 0, if (x0 + x1) != 0 }, if sgn((x0 + x1)) != 0 },\n"
        "0,\n"
        "{ { 0, if (x0 + x1) != 0 }, if sgn((x0 + x1)) != 0 },\n"
        "{ { 0, if (x0 + x1) != 0 }, if sgn((x0 + x1)) != 0 },\n"
        "0,\n"
        "0,\n"
        "0,\n"
        "0,\n"
        "{ 0, if ((x0 + x1) + 1) != 0 },\n"
        "{ 0, if ((x0 + x1) + 1) != 0 },\n"
        "0,\n"
        "{ 0, if ((x0 + x1) + 1) != 0 },\n"
        "{ 0, if ((x0 + x1) + 1) != 0 },\n"
        "0,\n"
        "0,\n"
        "0,\n"
        "0,\n"
        "{ 0, if ((x0 + x1) + 2) > 0 },\n"
        "{ 0, if ((x0 + x1) + 2) > 0 },\n"
        "0,\n"
        "{ 0, if ((x0 + x1) + 2) > 0 },\n"
        "{ 0, if ((x0 + x1) + 2) > 0 },\n"
        "0,\n"
        "0,\n"
        "0,\n"
        "0]";

    test("Defined_if function:", stream(F), f_str);
    test("Defined_if 1st derivative:", stream(dF), df_str);
    test("Defined_if 2nd derivative:", stream(ddF), ddf_str);

}

template<int i> struct ipp_times_xi
{
    typedef Expr<Times<Integer<i + 1> , Variable<i> > > simple_type;
};

void test_folding()
{
    cout << "\tTesting fold...\n";
    Array x_0(10, 1.0);
    auto sum_ixi = foldr<Plus, ipp_times_xi, 0, 9>();
    semt_test_f(
            sum_ixi,
            "(x0 + ((2 * x1) + ((3 * x2) + ((4 * x3) + ((5 * x4) + ((6 * x5) + ((7 * x6) + ((8 * x7) + ((9 * x8) + (10 * x9))))))))))",
            "1,\n\t2,\n\t3,\n\t4,\n\t5,\n\t6,\n\t7,\n\t8,\n\t9,\n\t10",
            x_0, 55.0, __LINE__);
    semt_test_f(
            foldl<Minus, ipp_times_xi, 0, 9>(),
            "(((((((((x0 - (2 * x1)) - (3 * x2)) - (4 * x3)) - (5 * x4)) - (6 * x5)) - (7 * x6)) - (8 * x7)) - (9 * x8)) - (10 * x9))",
            "1,\n\t-2,\n\t-3,\n\t-4,\n\t-5,\n\t-6,\n\t-7,\n\t-8,\n\t-9,\n\t-10",
            x_0, -53.0, __LINE__);
}
#else
#define test_diffable_expr();
#define test_partial();
#define test_simplifications();
#define test_folding();
#endif // DISABLE_ELABORATE_TEST
int main()
{
    //cout.setf(ios::scientific, ios::floatfield);
    //cout.setf(ios::right);
    //cout.precision(SEMT_STRPREC);

    cout << "Self testing... \n";
    test_diffable_expr();
    test_partial();
    test_simplifications();
    test_conditionals();
    test_folding();
    cout << "done\n";

    vector<SEMT_PRECISION> x_0(3);
    x_0[0] = -1.0;
    x_0[1] = 2.0021;
    x_0[2] = 2.3;

    cout << "sizeof( Var, Int, Expr<Var> , Expr<Int> ) = " << sizeof(Variable<1000>) << ", "
            << sizeof(Integer<1000>) << ", " << sizeof(x0) << ", " << sizeof(Ten) << "\n ";

    // These: Zero / Zero , Zero / (One * Zero + One * Zero)
    // produces the following error while compiling:
    // "error: incomplete type ‘SEMT::Divison_by_Zero’ used in nested name specifier"
    cout << Ten / (x0 + One + x0 * Zero - One - x0)
            << " should produce a compile time error, but assocativity is rather hard.\n";

    return 0;
}

