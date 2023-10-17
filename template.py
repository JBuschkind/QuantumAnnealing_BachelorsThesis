import dwave.cloud.solver
from dimod import ConstrainedQuadraticModel, Integer
from dwave.system import LeapHybridSampler, EmbeddingComposite, DWaveSampler

M = 100
N = 10


coefficients = {
    'p1': 2,
    'p2': 2,
    'q1': 2,
    'q2': 2,
    'c1': -4,
    'c2': -8,
    'c3': -4,
    'c4': -8
}

qubo = {}

for var in coefficients:
    qubo[(var,var)] = M * coefficients[var]

for var1 in coefficients:
    for var2 in coefficients:
        if var1 != var2:
            qubo[(var1,var2)] = 2 * M * coefficients[var1] * coefficients[var2]


for avr in coefficients:
    qubo[(var,var)] += N

qubo[('p1', 'q1')] = -4 * N
qubo[('p2', 'q2')] = -4 * N
qubo[('p1', 'q2')] = 4 * N
qubo[('p2', 'q1')] = 4 * N

#qubo['constant'] = 9

for var in coefficients:
    qubo[(var,var)] /= N

#scaling_factor = 1.0 / (N * N)
#for key in qubo:
#    qubo[key] *= scaling_factor

print(qubo)

sampler = EmbeddingComposite(DWaveSampler())

response = sampler.sample_qubo(qubo, num_reads=10)

for sample in response.samples():
    print(sample)
#cqm = ConstrainedQuadraticModel()
#cqm.from_quadratic_model((2*p2 + 2*p1*q1 + 2*q2 - 8*c2 - 4*c1 + p1 + q1 - 3) + (2*q1 + 2*p2*q2 + 2*p1 + 2*c2 - 8*c4 -4*c3 + p2*q1 + p1*q2 + c1 + 1) + (q2 + p2 + c3 + 2*c4 - 2))
#sampler = EmbeddingComposite(DWaveSampler())
#sampler.sample(cqm, num_reads=1000)
#answer = dwave.cloud.solver.BQMSolver.sample_qubo((2*p2 + 2*p1*q1 + 2*q2 - 8*c2 - 4*c1 + p1 + q1 - 3) + (2*q1 + 2*p2*q2 + 2*p1 + 2*c2 - 8*c4 -4*c3 + p2*q1 + p1*q2 + c1 + 1) + (q2 + p2 + c3 + 2*c4 - 2))
#print(answer)

