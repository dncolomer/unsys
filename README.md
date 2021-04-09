# ehsim
Entanglement Hypergraph Simulator

# Installation

## Full Installation

You can install this package like this:

```shell
pip3 install git+dncolomer/ehsim_prototype.git
```

Or you can clone this repository and run `pip3 install -e .` from inside the top level directory.

## Manual Installation

You will need to install prerequisite libraries before using `ehsim`. You can either install them like this,

```shell
pip3 install -r ./requirements.txt
```

Or you can install them manually:

```shell
pip3 install numpy hypernetx networkx
```

# Get Started

You can quickly get started with this library by creating your own quantum hypergraph.

Let's start with a simple example: [The Bell state](https://en.wikipedia.org/wiki/Bell_state).

![](https://upload.wikimedia.org/wikipedia/commons/f/fc/The_Hadamard-CNOT_transform_on_the_zero-state.png)

Let's create a new hypergraph:

```python
import ehsim.visualization as vis

from ehsim import Hypergraph
from ehsim.gates import H, X, SWAP, CX, CCX

# Let's start with 2 qubits
system = Hypergraph(2)
vis.print_raw(system)
```

Now, let's apply the Hadamard gate to qubit 0:

```python
system.rewrite(H,["q0"])
vis.print_raw(system)
```

And now we can apply a controlled-NOT gate:

```python
system.rewrite(CX,["q0","q1"])
vis.print_raw(system)
```

We can now factor the qubits and see what happens:

```python
circuit.factorQubits(["q0", "q1"])
vis.print_raw(system)
```

Congratulations! You just ran your first quantum hypergraph simulation!