## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=8, fig.height=8)

## ----setup---------------------------------------------------------------
library(PAFway)

set.seed(123)
#Make 300 nodes
nodes=paste("node", c(1:300))
#First 3 node names:
print(nodes[1:3])

#Assign them random GO terms
randomGO=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N")[sample(c(1:14), 300, replace=T)]
names(randomGO)=nodes
#First 3 nodes and associated GO terms:
print(randomGO[1:3])

#Make 1000 edges
edgesRandom=t(sapply(c(1:1000), function(i){
    nodes[sample(300, 2)]
 }))
#First 3 edges:
print(edgesRandom[1:3,])

## ------------------------------------------------------------------------
#Assign each node a second GO term, separated by a '_' symbol
randomGO2=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N")[sample(c(1:14), 300, replace=T)]
randomGO_multiple=sapply(c(1:300), function(i){paste(randomGO[i], randomGO2[i], sep="_")})
names(randomGO_multiple)=nodes

#print first 5 elements, as a demo
print(randomGO_multiple[1:5])

## ------------------------------------------------------------------------
#Interesting GO terms:
GO_interesting=c("B", "C", "D", "F")


## ------------------------------------------------------------------------
random_edge_weights=rnorm(length(edgesRandom[,1]), 1, 0.001)
print(random_edge_weights[1:5])
print(length(random_edge_weights))

## ------------------------------------------------------------------------
#This will run pafway, with no edge weights, for all the GO terms 
a=pafway(randomGO, edgesRandom, unique(randomGO))
print(a[1:5, 1:5])

## ------------------------------------------------------------------------
draw_network(a)

draw_heatmap(a)

## ------------------------------------------------------------------------
draw_network(a, adjMethod = "bonferroni")
draw_heatmap(a, adjMethod = "bonferroni")

## ------------------------------------------------------------------------
#This will run pafway, with no edge weights, for all the GO terms 
b=pafway_edge_weights(randomGO, cbind(edgesRandom, random_edge_weights), unique(randomGO))
print(b[1:5, 1:5])

## ------------------------------------------------------------------------
draw_network(b)

draw_heatmap(b)

## ------------------------------------------------------------------------
draw_network(b, adjMethod = "bonferroni")
draw_heatmap(b, adjMethod = "bonferroni")

## ------------------------------------------------------------------------
#Multiple GO terms are in: randomGO_multiple
b=pafway(randomGO_multiple, edgesRandom, unique(randomGO), exact=F)
print(b[1:5], b[1:5])

## ------------------------------------------------------------------------
#GO terms of interest are in GO_interesting
b=pafway(randomGO, edgesRandom, GO_interesting)
print(b[1:5], b[1:5])

## ------------------------------------------------------------------------
#This will run pafway, with no edge weights, for all the GO terms 
b=pafway_edge_weights(randomGO, cbind(edgesRandom, random_edge_weights), unique(randomGO), step=0.001, thresholdZero=0.0001)
print(b[1:5, 1:5])

