from sympy import diff, pprint, Interval, Abs, solve
from sympy.plotting import plot
from sympy.calculus.util import maximum
from sympy.abc import x, y
from math import factorial


def basic_checks(X0: float, Xk: float, funct, h: float) -> tuple[list, list]:
    if funct == None or X0 == None or Xk == None or h == None:
        raise (ValueError, "Недостатньо аргументів")

    n = int((Xk-X0)/h)

    X = []
    Y = []

    for i in range(n+1):
        X.append(X0+i*h)
        Y.append(round(funct.subs(x, X[i]), 4))

    return X, Y, n


def left_rectangles_method(X0: float, Xk: float, funct, h: float) -> float:
    X, Y, n = basic_checks(X0, Xk, funct, h)
    approximate_integral = 0

    for i in range(1, n+1):
        approximate_integral += Y[i-1]

    approximate_integral = approximate_integral * h

    return round(approximate_integral, 3)


def trapezoid_method(X0: float, Xk: float, funct, h: float) -> float:
    X, Y, n = basic_checks(X0, Xk, funct, h)
    approximate_integral = Y[0] + Y[-1]

    for i in range(1, n):
        approximate_integral += 2 * Y[i]

    approximate_integral = approximate_integral * (h/2)

    return round(approximate_integral, 3)


def simpson_method(X0: float, Xk: float, funct, h: float) -> float:
    X, Y, n = basic_checks(X0, Xk, funct, h)
    m = int(n/2)
    approximate_integral = Y[0] + Y[-1]

    for i in range(1, m+1):
        approximate_integral += 4 * Y[2*i-1]

    for j in range(1, m):
        approximate_integral += 2 * Y[2*j]

    approximate_integral = approximate_integral * (h/3)

    return round(approximate_integral, 3)


def RungeRomberg_specification(Ih1: float, Ih2: float, h:float, p: int) -> str:
    I = Ih2+((Ih2-Ih1)/(2**p-1))
    I = str(I) + " + O(" + str(h**(p+1)) + ")"
    return I


def numerical_integration(X0: float, Xk: float, funct, h1: float, h2: float) -> None:
    resultsRect1 = left_rectangles_method(X0, Xk, funct, h1)
    resultsTrap1 = trapezoid_method(X0, Xk, funct, h1)
    resultsSimps1 = simpson_method(X0, Xk, funct, h1)

    resultsRect2 = left_rectangles_method(X0, Xk, funct, h2)
    resultsTrap2 = trapezoid_method(X0, Xk, funct, h2)
    resultsSimps2 = simpson_method(X0, Xk, funct, h2)

    RungRombRect = RungeRomberg_specification(resultsRect1, resultsRect2, h1, 1)
    RungRombTrap = RungeRomberg_specification(resultsTrap1, resultsTrap2, h1, 2)
    RungRombSimps = RungeRomberg_specification(resultsSimps1, resultsSimps2, h1, 4)


    print(f"Наближені значення інтегралу \n∫({X0}->{Xk}){funct}\n\n\n")
    print(
        f"Обчислені методом прямокутників при\nh1={h1}: {resultsRect1}\nh2={h2}: {resultsRect2}\n\nУточнення за процедурою Рунге-Ронберга: {RungRombRect}\n\n\n")
    print(
        f"Обчислені методом трапецій при\nh1={h1}: {resultsTrap1}\nh2={h2}: {resultsTrap2}\n\nУточнення за процедурою Рунге-Ронберга: {RungRombTrap}\n\n\n")
    print(
        f"Обчислені методом Сімпсона при\nh1={h1}: {resultsSimps1}\nh2={h2}: {resultsSimps2}\n\nУточнення за процедурою Рунге-Ронберга: {RungRombSimps}\n")


if __name__ == "__main__":
    funct = (25+x**2)**(1/2)
    X0 = -2.0
    Xk = 2.0
    h1 = 1.0
    h2 = 0.5

    numerical_integration(X0, Xk, funct, h1, h2)
