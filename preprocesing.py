import math
from sympy import symbols, Eq, Add, Integer, Symbol, Mul, Pow, Wild
from dwave.system import DWaveSampler, EmbeddingComposite
from helper import *
from numpy import prod

debug = True

num = 15
numFactors = 2
l1 = math.floor((math.log2(num) + 1) / 2)
l2 = math.ceil((math.log2(num) + 1) / 2)

lengths = [4, 5]
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

eq1, coefficients = constructEquation(num, lengths)
print(eq1)
print(coefficients)

"""
if num & (1 << 0):
    eq1 = eq1.subs({Symbol('p1'): 1, Symbol('q1'): 1})
    p[0] = 1
    q[0] = 1
    coefficients.remove(Symbol('p1'))
    coefficients.remove(Symbol('q1'))
"""

eq1 = eq1.expand()

# Substitute xn**2 with xn where xn is every possible Symbol
eq1 = substituteExponent(eq1, coefficients, debug)

# Reduce the local Fields to have a maximum of 2 coefficients and save the constraints
assumptions = []
eq1, assumptions, coefficients = twoLocalFieldReduction(eq1, assumptions, coefficients, debug)


for val in coefficients:
    eq1 = eq1.replace(val, (1 - val) / 2)
eq1 = eq1.expand()

# local fields
h = {}
# couplings
J = {}

for val in eq1.args:
    if debug:
        print(f"{val}: ")
    # if isinstance(val, Integer):
        # eventually bias handling
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
# send and retrieve Data from DWave solvers
sampler = EmbeddingComposite(DWaveSampler())
response = sampler.sample_ising(h, J, num_reads=50, annealing_time=2000)

sampleNum = 0
for sample in response.samples():
    copyAssumptions = assumptions
    copyEq1 = eq1
    for val in sample:
        # print(val, sample[val])
        for elem in range(len(copyAssumptions)):
            copyAssumptions[elem] = copyAssumptions[elem].subs(val, math.floor((sample[val] - 1) / -2))
        copyEq1 = copyEq1.subs(val, math.floor((sample[val] - 1) / -2))
        if val.name[0] == "x":
            name = val.name
            name = name.replace("x", "")
            numbers[int(name)] = math.floor((sample[val] - 1) / -2)

    test = True
    for elem in copyAssumptions:
        if not (elem):
            test = False


    print(numbers)
    addedNumbers = [0] * len(lengths)
    numberBits = [""] * numFactors
    counter = 0
    counter2 = 0
    for i in lengths:
        for n in range(i):
            addedNumbers[counter2] += numbers[counter]*2**n
            numberBits[counter2] += numbers[counter]
            counter += 1
        counter2 += 1

# (test or sampleNum == 0) and
    if prod(addedNumbers) == num:
        print(f"Sample Number {sampleNum}")

        for elem in range(numFactors):
            print(f"{numberBits[elem]}, {addedNumbers[elem]}")

    sampleNum += 1
