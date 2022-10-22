import sys
import time
import numpy as np


from ase.calculators.emt import EMT
from ase.io import read
from ase import Atoms

from gofee import GOFEE
from gofee.surrogate.gpr import GPR
from gofee.candidates import CandidateGenerator, StartGenerator
from gofee.candidates import RattleMutation, RattleMutation2, PermutationMutation

# --- DEFINE SYSTEM ---
# (a) supercell
template = Atoms('', cell = np.eye(3)*12)

# (b) confinement box
confinement_cell = np.eye(3)*6
confinement_corner = np.array([3, 3, 3])
box = [confinement_corner, confinement_cell]

# (c) GOFEE atoms
stoichiometry = 6*[1] + 6*[6]

# --- GOFEE STRUCTURE GENERATION & MUTATION
# (a) start generator
sg = StartGenerator(template, stoichiometry, box)

# (b) rattle mutation
n_to_optimize = len(stoichiometry)
permutation = PermutationMutation(n_to_optimize, Npermute=2)
rattle = RattleMutation(n_to_optimize, Nrattle=3, rattle_range=2)
candidate_generator = CandidateGenerator([0.2, 0.2, 0.6],
                                         [sg, permutation, rattle])

# --- CALCULATOR ---
calc= EMT()

# --- GOFEE SEARCH ---
start_time = time.time()
search = GOFEE(calc=calc,
               startgenerator=sg,
               candidate_generator=candidate_generator,
               max_steps=20,
               population_size=5)

search.run()

end_time = time.time()
elapsed_time = end_time - start_time
print("Execution time: ", elapsed_time, 'seconds')
