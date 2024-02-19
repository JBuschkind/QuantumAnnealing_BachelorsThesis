from dwave.system import DWaveSampler, EmbeddingComposite
from helper import *


debug = True

num = 4489
l_num = num.bit_length()
numFactors = 2
num_bits = []
#l_p = math.floor((math.log2(num) + 1) / 2)
#l_q = math.floor((math.log2(num) + 1) / 2)
l_p = 7
l_q = 7

p = []
q = []

lengths = [l_p, l_q]

for i in range(l_num):
    num_bits.append((num >> i) & 1)

for i in range(l_p):
    p.append(symbols(f"p{i}"))

for i in range(l_q):
    q.append(symbols(f"q{i}"))


if debug:
    counterLength = 0
    for i in lengths:
        print(f"Length Number {counterLength}: {i}")
        counterLength += 1

# l1 = 2
# l2 = 2

numbers = [None] * sum(lengths)

"""
coefficients = []
for i in range(1, l1 + 1):
    globals()[f"p{i}"] = symbols(f"p{i}")
    coefficients.append(symbols(f"p{i}"))

for i in range(1, l2 + 1):
    globals()[f"q{i}"] = symbols(f"q{i}")
    coefficients.append(symbols(f"q{i}"))

eq1 = Eq(Pow(num - Add(*(Mul(Integer(2 ** (n - 1)), Symbol(f"p{n}")) for n in range(l1 + 1))) * Add(*(Mul(Integer(2 ** (n - 1)), Symbol(f"q{n}")) for n in range(l2 + 1))), 2), y)
print(eq1)
print(coefficients)
"""

eq, coefficients = constructEquation(num, lengths)
print(eq)
print(coefficients)


if num_bits[0]:
    eq = eq.subs({Symbol('p0'): 1, Symbol('q0'): 1})
    p[0] = 1
    q[0] = 1
    coefficients.remove(Symbol("p0"))
    coefficients.remove(Symbol("q0"))
    if debug:
        print("")
        print("Removing the First Bits")
        print(eq)
        print(coefficients)
        print("")



if num_bits[l_num - 1]:
    eq = eq.subs({Symbol(f"p{l_p - 1}"): 1, Symbol(f"q{l_q - 1}"): 1})
    p[l_p - 1] = 1
    q[l_q - 1] = 1
    coefficients.remove(Symbol(f"p{l_p - 1}"))
    coefficients.remove(Symbol(f"q{l_q - 1}"))
    if debug:
        print("")
        print("Removing the Last Bits")
        print(eq)
        print(coefficients)
        print("")


substitutions, coefficients = createSubstitutionList(eq, coefficients, p, q, debug)
print(substitutions)
eq = eq.expand()
print(eq)
# Substitute xn**2 with xn where xn is every possible Symbol
eq = substituteExponent(eq, coefficients, debug)

print(eq)

# Reduce the local Fields to have a maximum of 2 coefficients and save the constraints
assumptions = []
eq = twoLocalFieldReduction(eq, substitutions, debug)


for val in coefficients:
    eq = eq.replace(val, (1 - val) / 2)
eq = eq.expand()

# local fields
h = {}
# couplings
J = {}

for val in eq.args:
    if debug:
        print(f"{val}: ")
    if isinstance(val, Mul):
        if len(val.args) == 2:
            if isinstance(val.args[0], Symbol):
                if debug:
                    print(f"1,  {val.args[0]},  {val.args[1]}")
                J[val.args[0], val.args[1]] = 1.0
            else:
                if debug:
                    print(f"{val.args[0]},  {val.args[1]}")
                h[val.args[1]] = val.args[0] / 2
        elif len(val.args) == 3:
            if debug:
                print(f"{val.args[0]},  {val.args[1]},  {val.args[2]}")
            J[val.args[1], val.args[2]] = val.args[0] / 2

if debug:
    print("")
    print(f"local Fields: {h}")
    print(f"couplings: {J}")
    print("")
# print(dict2ltxtab(J, coefficients))
# send and retrieve Data from DWave solvers


sampler = EmbeddingComposite(DWaveSampler())
response = sampler.sample_ising(h, J, num_reads=100, annealing_time=2000)


print(response)
print(response.record)
print(response.record[0][0])

for n in range(len(response.record)):
    copy_p = p.copy()
    copy_q = q.copy()
    for i in range(len(copy_p)):
        if isinstance(copy_p[i], Symbol):
            if response.samples()[n][copy_p[i]] == 1:
                copy_p[i] = 0
            else:
                copy_p[i] = 1

    for i in range(len(copy_q)):
        if isinstance(copy_q[i], Symbol):
            if response.samples()[n][copy_q[i]] == 1:
                copy_q[i] = 0
            else:
                copy_q[i] = 1

    num_p = 0
    num_q = 0
    counter = 0
    for elem in copy_p:
        num_p += elem * 2**counter
        counter += 1

    counter = 0
    for elem in copy_q:
        num_q += elem * 2**counter
        counter += 1

    print(f"num: {n}, p: {num_p}, q: {num_q}, occ: {response.record[n][2]}")
