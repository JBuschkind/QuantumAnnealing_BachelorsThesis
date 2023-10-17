from dwave.system import DWaveSampler, EmbeddingComposite

h = {'a': 58, 'b': 50, 'c': 12, 'd': -80}
J = {('a', 'b'): 25, ('a', 'c'): -6, ('a', 'd'): -64, ('b', 'c'): 2, ('b', 'd'): -64, ('c', 'd'): 16}

sampler = EmbeddingComposite(DWaveSampler())

response = sampler.sample_ising(h, J, num_reads=200)

for sample in response.samples():
    print(sample)

print(response)