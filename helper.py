from sympy import Wild, Symbol, symbols, Eq, Add, Mul, Integer, Pow


def sortingset(ele):
    return len(ele.args)


def constructEquation(number, lengths):
    coefficients = []

    def addNumber(length, letter):
        neweq = Add(*(Mul(Integer(2 ** n), Symbol(f"{letter}{n}")) for n in range(length)))
        for i in range(length):
            globals()[f"{letter}{i}"] = symbols(f"{letter}{i}")
            coefficients.append(symbols(f"{letter}{i}"))
        return neweq

    eq = Pow(number - Mul(addNumber(lengths[0], "p"), addNumber(lengths[1], "q")), 2)

    return eq, coefficients


def createSubstitutionList(eq, coefficients, p, q, debug=False):
    copyCoefficients = coefficients
    substitutions = []
    x = Wild("x")
    runningNumber = 1
    print(eq)
    for i in p:
        for n in q:
            if isinstance(i, Symbol) and isinstance(n, Symbol):
                print(f"{i}, {n}")
                query = x * i * n
                result = eq.find(query)
                remove = []
                for elem in result:
                    if isinstance(elem, Integer):
                        remove.append(elem)
                result = result.difference(remove)
                if result:
                    print(result)

                    globals()[f"t{runningNumber}"] = symbols(f"t{runningNumber}")
                    coefficients.append(symbols(f"t{runningNumber}"))
                    x3 = symbols(f"t{runningNumber}")
                    runningNumber += 1
                    substitutions.append((i, n, x3))

    if debug:
        print("")
        print(substitutions)
        print("")
        print("")

    return substitutions, coefficients


def twoLocalFieldReduction(eq, substitutions, debug=False):
    a = Wild("a", properties=[lambda x: isinstance(x, Symbol)])
    x = Wild("x")
    for elem in substitutions:
        # print(elem)
        # print("")
        query = x * a * elem[0] * elem[1]
        result = sorted(eq.find(query), key=sortingset, reverse=True)
        print(result)
        x1 = elem[0]
        x2 = elem[1]

        x4 = elem[2]
        for val in result:
            x3 = 1
            for var in val.args:
                if isinstance(var, Symbol) and var != elem[0] and var != elem[1]:
                    x3 = Mul(x3, var)

            print(f"{x1}, {x2}, {x3}, {x4}")
            print("")

            if isinstance(val.args[0], Integer) and val.args[0] < 0:
                eq = eq.subs(-x1 * x2 * x3, (-x4 * x3 + 2 * (x1 * x2 - 2 * x1 * x4 - 2 * x2 * x4 + 3 * x4)))
            else:
                eq = eq.subs(x1 * x2 * x3, (x4 * x3 + 2 * (x1 * x2 - 2 * x1 * x4 - 2 * x2 * x4 + 3 * x4)))

            # print(eq)
            # print("")
            # print("")
        # result = sorted(eq.find(query), key=sortingset, reverse=True)
        # print(result)

    return eq

def substituteExponent(eq, coefficients, debug=False):
    # Substitute qn**2 with qn, n being any number
    substitutions = []
    for elem in coefficients:
        substitutions.append((elem ** 2, elem))
    eq = eq.subs(substitutions)

    eq.expand()
    return eq

def dict2ltxtab(d: dict, coefficients, bare=False, headrow = None):
    if not bare:
        print(r"\begin{center}")
        print(r"\begin{tabular}{|c|c|}")
        print(r"\hline")

    for n in range(len(coefficients)-1):
        out = ""
        for i in range(0, n):
            out = out + "&"
        for p in range(n+1, len(coefficients)):
            if d.get((coefficients[n], coefficients[p])) is not None:
                out = out + " " + str(d.get((coefficients[n], coefficients[p]))) + "&"
            elif d.get((coefficients[p], coefficients[n])) is not None:
                out = out + " " + str(d.get((coefficients[p], coefficients[n]))) + "&"
            else:
                out = out + " 0 &"
        print(out + "\\\\")
    if not bare:
        print(r"\hline")
        print(r"\end{tabular}")
        print(r"\end{center}")