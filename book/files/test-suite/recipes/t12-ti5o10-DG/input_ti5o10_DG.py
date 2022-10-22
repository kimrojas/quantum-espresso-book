import sys
import time
import numpy as np
import os

from ase.io import read
from ase.calculators.dftb import Dftb
from ase import Atoms

from gofee import GOFEE
from gofee.surrogate.gpr import GPR
from gofee.candidates import CandidateGenerator, StartGenerator
from gofee.candidates import RattleMutation, RattleMutation2, PermutationMutation
from gofee.utils import OperationConstraint

# --- SET THREADS ---
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-t', '--threads', type=int, default=1)
args = parser.parse_args()


# --- CALCULATOR ---
os.environ['DFTB_PREFIX'] = '/home/krojas/share/lib/slakos/matsci'
os.environ['ASE_DFTB_COMMAND'] = f'OMP_NUM_THREADS={args.threads} dftb+ > PREFIX.out'
# os.environ['ASE_DFTB_COMMAND'] = f'dftb+ > PREFIX.out'

calc = Dftb(label='ti5o10',
            Hamiltonian_SCC='No',
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_Ti='"d"',
            Hamiltonian_MaxAngularMomentum_O='"p"',
            Hamiltonian_Charge='0.000000',
            Hamiltonian_Filling='Fermi {',
            Hamiltonian_Filling_empty='Temperature [Kelvin] = 0.000000',
            kpts=[1, 1, 1])



# --- DEFINE SYSTEM ---
# (a) supercell
template = Atoms('',
            cell=[20,20,20],
            pbc=[0, 0, 0])

# (b) confinement box
confinement_cell = 4*np.eye(3)
confinement_corner = np.array((8.0, 8.0, 8.0))
box = [confinement_corner, confinement_cell]

# (c) GOFEE atoms
stoichiometry = 5*[22]+10*[8]


# --- GOFEE STRUCTURE GENERATION & MUTATION
# (a) start generator
sg = StartGenerator(template, stoichiometry, box)

# (b) rattle mutation
n_to_optimize = len(stoichiometry)
rattle = RattleMutation(n_to_optimize, Nrattle=3, rattle_range=4)
candidate_generator = CandidateGenerator(probabilities=[0.2, 0.8],
                                         operations=[sg, rattle])


# --- GOFEE SEARCH ---
start_time = time.time()
search = GOFEE(calc=calc,
                startgenerator=sg,
                candidate_generator=candidate_generator,
                max_steps=50,
                population_size=5)

search.run()

end_time = time.time()
elapsed_time = end_time - start_time
print("Execution time: ", elapsed_time, 'seconds')
