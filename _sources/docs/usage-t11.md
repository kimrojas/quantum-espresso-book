# HB sheet system

A GOFEE calculation of `HB sheet` with DFTB+ as a fast and relatively accurate evaluator. Here we use a parallelized GOFEE implementation and a **Non-parallelized** DFTB+ calculation (Parallelization mixing will be discussed in the future).

```{warning}
Parameters used are for testing and tutorial purposes only. It is not for production.
```

```python
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

calc = Dftb(label='hb',
            Hamiltonian_SCC='No',
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_B='"p"',
            Hamiltonian_MaxAngularMomentum_H='"s"',
            Hamiltonian_Charge='0.000000',
            Hamiltonian_Filling='Fermi {',
            Hamiltonian_Filling_empty='Temperature [Kelvin] = 0.000000',
            kpts=[8, 8, 1])



# --- DEFINE SYSTEM ---
# (a) supercell
template = read('../../scaffold.traj')
# template = Atoms('', cell = np.eye(3)*12)

# (b) confinement box
confinement_cell = template.get_cell()
confinement_corner = np.array((0.000, 0.000, 0.000))
box = [confinement_corner, confinement_cell]
box_constraint = OperationConstraint(box=box)

# (c) GOFEE atoms
stoichiometry = 4*[1]+3*[5]


# --- GOFEE STRUCTURE GENERATION & MUTATION
# (a) start generator
sg = StartGenerator(template, stoichiometry, box)

# (b) rattle mutation
n_to_optimize = len(stoichiometry)
permutation = PermutationMutation(n_to_optimize, Npermute=2)
rattle = RattleMutation(n_to_optimize, Nrattle=3, rattle_range=4)
candidate_generator = CandidateGenerator([0.3, 0.3, 0.4],
                                         [sg, permutation, rattle])


# --- GOFEE SEARCH ---
start_time = time.time()
search = GOFEE(calc=calc,
               startgenerator=sg,
               candidate_generator=candidate_generator,
               max_steps=20,
               population_size=5,
               position_constraint=box_constraint)

search.run()

end_time = time.time()
elapsed_time = end_time - start_time
print("Execution time: ", elapsed_time, 'seconds')

```