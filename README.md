# Graph state generation in the presence of noise
This repository contains sample programs to carry out various graph state operations in the presence of noise, including CZ gates, local complementations and measurement in different bases (X,Y,Z). This includes functions to calculate the fidelity of the obtained graph state after some set of operations. Specifically, there are 3 files to carry out simulations to compare the performance of two graph state distribution algorithms: the subgraph complementation protocol and the factory node + teleportation protocol.

1. _ghz_simulate.py_ generates graphs of error probability/number of qubits vs fidelity for the two protocols when creating a GHZ state.
2. _bip_simulate.py_ generates graphs of error probability/number of qubits vs fidelity for the two protocols when creating a complete bipartite graph state.
3. _ghz_simulate_parallel.py_ generates graphs of error probability/number of qubits vs fidelity for the two protocols when creating a GHZ state incorporating parallelization in the protocol.
