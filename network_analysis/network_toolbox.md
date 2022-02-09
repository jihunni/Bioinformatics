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
file_name = "proximity/data/cor0.6.sif"
network = wrappers.get_network(file_name, only_lcc = True)
nodes_from = ["A", "C"]
nodes_to = ["B", "D", "E"]
d, z, (mean, sd) = wrappers.calculate_proximity(network, nodes_from, nodes_to, min_bin_size = 2, seed=452456)

mt = ["ENSG00000168827","ENSG00000254093","ENSG00000145982","ENSG00000074582","ENSG00000131368","ENSG00000072506","ENSG00000214026","ENSG00000179271","ENSG00000100890","ENSG00000158042","ENSG00000136522","ENSG00000168924","ENSG00000154719","ENSG00000185608","ENSG00000171861","ENSG00000137513","ENSG00000140521","ENSG00000086504","ENSG00000166902","ENSG00000135776","ENSG00000115286","ENSG00000110717","ENSG00000180992","ENSG00000133983","ENSG00000259494","ENSG00000163607","ENSG00000029639","ENSG00000172172","ENSG00000256525"]
purine = ["ENSG00000128059","ENSG00000138031","ENSG00000136877","ENSG00000163655","ENSG00000241186","ENSG00000130348","ENSG00000138363","ENSG00000035687","ENSG00000143774","ENSG00000117118","ENSG00000170190","ENSG00000159131","ENSG00000095059","ENSG00000239900","ENSG00000100714","ENSG00000178921","ENSG00000243678","ENSG00000183955","ENSG00000128050"]
glycolysis = ["ENSG00000102144","ENSG00000105220","ENSG00000163931","ENSG00000153574","ENSG00000171314","ENSG00000074800"]
```
