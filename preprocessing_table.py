import math

from dwave.system import EmbeddingComposite, DWaveSampler

from helper import *

debug = True

num = 143
l_num = num.bit_length()
l_p = math.ceil(math.log2(num) / 2)
l_q = math.ceil(math.log2(num) / 2)

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
    coefficients.remove(Symbol("p0"))
    coefficients.remove(Symbol("q0"))

if num_bits[l_num - 1]:
    p[l_p - 1] = 1
    q[l_q - 1] = 1
    coefficients.remove(Symbol(f"p{l_p - 1}"))
    coefficients.remove(Symbol(f"q{l_q - 1}"))

print(p)
print(q)
for i in range(l_p):
    for n in range(l_q):
        num_exp[i + n] = Add(num_exp[i + n], Mul(p[i], q[n]))

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
                num_exp[i + n] = Add(num_exp[i + n], symbols(f"c{counter}"))
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
            lhs = Add(lhs, Mul(num_exp[i + n], Pow(2, n)))
            if new_counter < counter:
                rhs = Add(rhs, Add(Mul(symbols(f"c{new_counter}"), Pow(2, 2 + n)), Mul(num_bits[i + n], Pow(2, n))))
                new_counter += 1
            else:
                rhs = Add(rhs, Mul(num_bits[i + n], Pow(2, n)))
                last = True
        if last:
            rest = i + l_block
            while rest < l_num:
                rhs = Add(rhs, Mul(num_bits[rest], Pow(2, rest - i)))
                rest += 1
        eqs.append(lhs - rhs)

eq = Integer(0)

for elem in eqs:
    eq = Add(eq, Pow(elem, 2))

substitutions = []
substitutions, coefficients = createSubstitutionList(eq, coefficients, p, q, debug)
eq = eq.expand()

eq = substituteExponent(eq, coefficients, debug)

print(eq)
print("")
print("")

eq = twoLocalFieldReduction(eq, substitutions, debug)

print(eq)
print(coefficients)
for val in coefficients:
    eq = eq.replace(val, (1 - val) / 2)
eq = eq.expand()

print(eq)
eq = Mul(eq, 2)
print(eq)

# local fields
h = {}
# couplings
J = {}

for val in eq.args:
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
                h[val.args[1]] = val.args[0]
        elif len(val.args) == 3:
            if debug:
                print(f"{val.args[0]},  {val.args[1]},  {val.args[2]}")
            J[val.args[1], val.args[2]] = val.args[0]

if debug:
    print("")
    print(f"local Fields: {h}")
    print(f"couplings: {J}")
    print("")

exit()
sampler = EmbeddingComposite(DWaveSampler())
response = sampler.sample_ising(h, J, num_reads=25, annealing_time=200)

print(response)

for n in range(len(response.samples())):
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

    print(f"num: {n}, p: {num_p}, q: {num_q}")

exit()
print(response)
print(response.samples()[0])
for i in range(len(p)):
    if isinstance(p[i], Symbol):
        if response.samples()[0][p[i]] == 1:
            p[i] = 0
        else:
            p[i] = 1

for i in range(len(q)):
    if isinstance(q[i], Symbol):
        if response.samples()[0][q[i]] == 1:
            q[i] = 0
        else:
            q[i] = 1

num_p = 0
num_q = 0
counter = 0
for elem in p:
    num_p += elem * 2**counter
    counter += 1

counter = 0
for elem in q:
    num_q += elem * 2**counter
    counter += 1

print(f"p: {num_p},q: {num_q}")