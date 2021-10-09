install.packages("igraph") # Download and install the package
library(igraph) # Load package
demo(package="igraph")

# Create networks
g1 <- graph( edges=c(1,2, 2,3, 3, 1), n=3, directed=F )
plot(g1) # A simple plot of the network - we'll talk more about plots later
class(g1) ## [1] "igraph"
g1
    # IGRAPH 84b8fdd U--- 3 3 -- 
    #   + edges from 84b8fdd:
    #   [1] 1--2 2--3 1--3

# Now with 10 vertices, and directed by default:
g2 <- graph( edges=c(1,2, 2,3, 3, 1), n=10 )
plot(g2)
g2
    # IGRAPH 9b4d87a D--- 10 3 -- 
    #   + edges from 9b4d87a:
    #   [1] 1->2 2->3 3->1


g3 <- graph( c("John", "Jim", "Jim", "Jill", "Jill", "John")) # named vertices
# When the edge list has vertex names, the number of nodes is not needed
plot(g3)
g3
    # IGRAPH d43af75 DN-- 3 3 -- 
    #   + attr: name (v/c)
    # + edges from d43af75 (vertex names):
    #   [1] John->Jim  Jim ->Jill Jill->John

g4 <- graph( c("John", "Jim", "Jim", "Jack", "Jim", "Jack", "John", "John"),
             isolates=c("Jesse", "Janis", "Jennifer", "Justin") )
  # In named graphs we can specify isolates by providing a list of their names.
plot(g4, edge.arrow.size=.5, vertex.color="gold", vertex.size=15,
     vertex.frame.color="gray", vertex.label.color="black",
     vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.2)

# Small graphs can also be generated with a description of this kind: - for undirected tie, +- or -+
#   for directed ties pointing left & right, ++ for a symmetric tie, and “:” for sets of vertices.
plot(graph_from_literal(a---b, b---c)) # the number of dashes doesn't matter
plot(graph_from_literal(a--+b, b+--c))
plot(graph_from_literal(a+-+b, b+-+c))
plot(graph_from_literal(a:b:c---c:d:e))

gl <- graph_from_literal(a-b-c-d-e-f, a-g-h-b, h-e:f:i, j)
plot(gl)


E(g4) # The edges of the object
V(g4) # The vertices of the object
g4[]
g4[1,]

V(g4)$name # automatically generated when we created the network.

V(g4)$gender <- c("male", "male", "male", "male", "female", "female", "male")
E(g4)$type <- "email" # Edge attribute, assign "email" to all edges
E(g4)$weight <- 10
# Edge weight, setting all existing edges to 10
edge_attr(g4)
vertex_attr(g4)
graph_attr(g4)
graph_attr_names(g4)

plot(g4, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
     vertex.color=c( "pink", "skyblue")[1+(V(g4)$gender=="male")] )

#simplify our graph to remove loops & multiple edges between the same nodes.
g4s <- simplify( g4, remove.multiple = T, remove.loops = F,
                 edge.attr.comb=c(weight="sum", type="ignore") )
plot(g4s, vertex.label.dist=1.5)

## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)
print(g, e=TRUE, v=TRUE)
plot(g)

## The opposite operation
as_data_frame(g, what="vertices")
as_data_frame(g, what="edges")

#sub-graph############################################
#1. Cliques (complete subgraphs of an undirected graph)
cliques(net) # list of cliques
sapply(cliques(net), length) # clique sizes
largest_cliques(net) # cliques with max number of nodes

vcol <- rep("grey80", vcount(net))
vcol[unlist(largest_cliques(net))] <- "gold"
plot(as.undirected(net), vertex.label=V(net)$name, vertex.color=vcol)


# 2. community detection
#Community detection based on edge betweenness (Newman-Girvan)
net = gl
ceb <- cluster_edge_betweenness(net) # list
dendPlot(ceb, mode="hclust")
plot(ceb, net)

class(ceb) ## [1] "communities"
length(ceb) # number of communities
membership(ceb) # community membership for each node
modularity(ceb) # how modular the graph partitioning is
crossing(ceb, net) # boolean vector: TRUE for edges across communities

# Community detection based on based on propagating labels
clp <- cluster_label_prop(net)
plot(clp, net)

# Community detection based on greedy optimization of modularity
cfg <- cluster_fast_greedy(as.undirected(net))
plot(cfg, as.undirected(net))

V(net)$community <- cfg$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(net, vertex.color=colrs[V(net)$community])


# Community strucure via short random walks
cwalk = walktrap.community(graph, weights = E(graph)$weight, steps = 4, merges =
                     TRUE, modularity = FALSE, labels = TRUE)
cwalk = walktrap.community(net, weights = E(net)$weight)
plot(cwalk, net)

# 3. K-core decomposition
  # The k-core is the maximal subgraph in which every node has degree of at least k. This also means
  # that the (k+1)-core will be a subgraph of the k-core.
  # The result here gives the coreness of each vertex in the network. A node has coreness D if it belongs
  # to a D-core but not to (D+1)-core.
kc <- coreness(net, mode="all")
plot(net, vertex.size=kc*6, vertex.label=kc, vertex.color=colrs[kc])

assortativity_nominal(net, V(net)$media.type, directed=F)
assortativity(net, V(net)$audience.size, directed=F)
assortativity_degree(net, directed=F)
