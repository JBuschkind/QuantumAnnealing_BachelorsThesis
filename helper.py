from sympy import Wild, Symbol, symbols, Eq, Add, Mul, Integer, Pow


def sortingset(ele):
    return len(ele.args)


def constructEquation(number, lengths):
    coefficients = []

    def addNumber(length):
        neweq = Add(*(Mul(Integer(2 ** (n)), Symbol(f"x{len(coefficients) + n}")) for n in range(length)))
        for i in range(len(coefficients), len(coefficients) + length):
            globals()[f"x{i}"] = symbols(f"x{i}")
            coefficients.append(symbols(f"x{i}"))
        return neweq

    eq1 = Pow(number - Mul( *(addNumber(i) for i in lengths)), 2)

    return eq1, coefficients

def twoLocalFieldReduction(eq1, assumptions, coefficients, debug=False):
    a = Wild("a", properties=[lambda x: isinstance(x, Symbol)])
    b = Wild("b", properties=[lambda x: isinstance(x, Symbol)])
    c = Wild("c", properties=[lambda x: isinstance(x, Symbol)])
    x = Wild("x")

    if debug:
        print("")
        print("Before Reduction:")
        print(eq1)
        print("")

    query = x * a * b * c
    runningNumber = 0
    result = sorted(eq1.find(query), key=sortingset, reverse=True)
    print(result)
    actions = 0
    while result:

        val = result.pop(0)
        if debug:
            print(f"Action {actions}: {val}")

        x1 = 1
        x2 = 1
        x3 = 1
        n = 1

        for elem in val.args:
            if isinstance(elem, Symbol):
                if n == 1:
                    x1 = elem
                    n += 1
                elif n == 2:
                    x2 = elem
                    n += 1
                else:
                    x3 = x3 * elem
        x4 = symbols(f"c{runningNumber}")
        coefficients.append(symbols(f"c{runningNumber}"))
        runningNumber += 1

        if debug:
            print(f"No:{actions} first: {x1}, second: {x2}, third: {x3}, fourth: {x4}")
            actions += 1
        eq1 = eq1.subs(x1 * x2 * x3, (x4 * x3 + 2 * (x1 * x2 - 2 * x1 * x4 - 2 * x2 * x4 + 3 * x4)))
        assumptions.append(Eq(x4, x1 * x2))
        eq1 = eq1.expand()

        result = sorted(eq1.find(query), key=sortingset, reverse=True)

    if debug:
        print("")
        print("After reduction:")
        print(eq1)
        print("")

    return eq1, assumptions, coefficients


def substituteExponent(eq1, coefficients, debug=False):
    # Substitute qn**2 with qn, n being any number
    substitutions = []
    for elem in coefficients:
        substitutions.append((elem ** 2, elem))
    eq1 = eq1.subs(substitutions)

    eq1.expand()
    return eq1
