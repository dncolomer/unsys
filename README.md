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
from ehsim import Hypergraph
from ehsim.gates import H, X

# start with 2 qubits:
circuit = Hypergraph(2)
```

Now, let's apply the Hadamard gate to qubit 0:

```python
circuit.apply("q0", H)
```

And now we can apply a controlled-NOT gate:

```python
circuit.apply("q1", X, controls=["q0"])
```

We can now factor the qubits and see what happens:

```python
circuit.factorQubits(["q0", "q1"])
```

Congratulations! You just ran your first quantum hypergraph simulation!