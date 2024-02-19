from sympy import Wild, Symbol, symbols, Eq, Add, Mul, Integer, Pow, sign

exp = Mul(-8, symbols("x1"))

print(exp.is_negative)

