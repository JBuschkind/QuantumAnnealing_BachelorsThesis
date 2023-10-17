from sympy import symbols, Eq, Add, Integer, Symbol, Mul, Pow, Wild, simplify
from helper import *

debug = True

num = 143
l_num = num.bit_length()
l_p = 4
l_q = 4

l_block = 2

if debug:
    print(f"Number: {num}, Length: {l_num}")

num_exp = []
num_bits = []
for i in range(l_num):
    num_bits.append((num >> i) & 1)
p = []
q = []
coefficients = []

for i in range(l_p):
    globals()[f"p{i}"] = symbols(f"p{i}")
    coefficients.append(symbols(f"p{i}"))
    p.append(symbols(f"p{i}"))

for i in range(l_q):
    globals()[f"q{i}"] = symbols(f"q{i}")
    coefficients.append(symbols(f"q{i}"))
    q.append(symbols(f"q{i}"))

for i in range(l_num):
    num_exp.append(Integer(0))

if num_bits[0]:
    p[0] = 1
    q[0] = 1

if num_bits[l_num-1]:
    p[l_p-1] = 1
    q[l_q-1] = 1

for i in range(l_p):
    for n in range(l_q):
        num_exp[i+n] = Add(num_exp[i+n], Mul(p[i], q[n]))

start = 0
while num_exp[start] == 1 or num_exp[start] == 0:
    start += 1

counter = 1
first = True
for i in range(start, l_num, l_block):
    if not first:
        if i + l_block <= l_num:
            for n in range(l_block):
                globals()[f"c{counter}"] = symbols(f"c{counter}")
                coefficients.append(symbols(f"c{counter}"))
                num_exp[i+n] = Add(num_exp[i+n], symbols(f"c{counter}"))
                counter += 1
    else:
        first = False


start = 0
while num_exp[start] == 1 or num_exp[start] == 0:
    start += 1

new_counter = 1
eqs = []
last = False
for i in range(start, l_num, l_block):
    if i + l_block <= l_num:
        lhs = Integer(0)
        rhs = Integer(0)
        for n in range(l_block):
            lhs = Add(lhs, Mul(num_exp[i+n], Pow(2, n)))
            if new_counter < counter:
                rhs = Add(rhs, Add(Mul(symbols(f"c{new_counter}"), Pow(2, 2+n)), Mul(num_bits[i+n], Pow(2, n))))
                new_counter += 1
            else:
                rhs = Add(rhs, Mul(num_bits[i + n], Pow(2, n)))
                last = True
        if last:
            rest = i + l_block
            while rest < l_num:
                rhs = Add(rhs, Mul(num_bits[rest], Pow(2, rest-i)))
                rest += 1
        eqs.append(lhs - rhs)

eq = Integer(0)

for elem in eqs:
    eq = Add(eq, Pow(elem, 2))

substitutions = []

a = Wild("a", properties=[lambda x: isinstance(x, Symbol)])
b = Wild("b", properties=[lambda x: isinstance(x, Symbol)])
x = Wild("x")

query = x * a * b
result = sorted(eq.find(query), key=sortingset, reverse=True)
runningNumber = 1

for val in result:
    #print(val.args)

    x1 = 1
    x2 = 1
    n = 1

    for elem in val.args:
        if isinstance(elem, Symbol):
            if n == 1:
                x1 = elem
                n += 1
            elif n == 2:
                x2 = elem
                n += 1
    x3 = symbols(f"t{runningNumber}")
    runningNumber += 1
    substitutions.append((x1, x2, x3))

substitutions = [(Symbol("p1"), Symbol("q1"), Symbol("t1")), (Symbol("p1"), Symbol("q2"), Symbol("t2")), (Symbol("p2"), Symbol("q2"), Symbol("t3")), (Symbol("p2"), Symbol("q1"), Symbol("t4"))]
print(substitutions)
print(eq)
eq = eq.expand()
#print(eq)

eq = substituteExponent(eq, coefficients, debug)

print(eq)
a = Wild("a", properties=[lambda x: isinstance(x, Symbol)])
x = Wild("x")
for elem in substitutions:
    print(elem)
    query = x * a * elem[0] * elem[1]
    result = sorted(eq.find(query), key=sortingset, reverse=True)
    print(result)
    x1 = elem[0]
    x2 = elem[1]

    x4 = elem[2]
    while result:

        for val in result:
            x3 = 1
            for var in val.args:
                if isinstance(var, Symbol) and var != elem[0] and var != elem[1]:
                    x3 = Mul(x3, var)

            print(f"{x1}, {x2}, {x3}, {x4}")
            if isinstance(val.args[0], Integer) and val.args[0] < 0:
                eq = eq.subs(x1 * x2 * x3, (x4 * x3 + 2 * (x1 * x2 - 2 * x1 * x4 - 2 * x2 * x4 + 3 * x4)))
            else:
                eq = eq.subs(x1 * x2 * x3, (x4 * x3 + 2 * (x1 * x2 - 2 * x1 * x4 - 2 * x2 * x4 + 3 * x4)))

            #print(eq)
            result.remove(val)
        #result = sorted(eq.find(query), key=sortingset, reverse=True)
        #print(result)
eq = eq.expand()
print(eq)

for val in coefficients:
    eq = eq.replace(val, (1 - val) / 2)
eq = eq.expand()

print(eq)