# Installation
```
conda create --name py27 python=2.7
conda activate py27
pip install decorator==3.3.1
pip install networkx==2.0
pip install numpy
pip install scipy
```

# To calculate proximity
general
```
from toolbox import wrappers
file_name = "data/toy.sif"
network = wrappers.get_network(file_name, only_lcc = True)
nodes_from = ["A", "C"]
nodes_to = ["B", "D", "E"]
d, z, (mean, sd) = wrappers.calculate_proximity(network, nodes_from, nodes_to, min_bin_size = 2, seed=452456)
```

running
```
from toolbox import wrappers
file_name = "proximity/data/cor0.6.sif""
network = wrappers.get_network(file_name, only_lcc = True)
nodes_from = ["A", "C"]
nodes_to = ["B", "D", "E"]
d, z, (mean, sd) = wrappers.calculate_proximity(network, nodes_from, nodes_to, min_bin_size = 2, seed=452456)
```
