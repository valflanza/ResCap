library("tidyverse")
library("MCL")
library("igraph")




##### AlleleNetwork Alternative ###############################################

Net.abr = read_delim("Net.abr.net", col_names = FALSE, delim = "\t")
Net.bac = read_delim("Net.bac.net", col_names = FALSE, delim = "\t")
Net.rel = read_delim("Net.rel.net", col_names = FALSE, delim = "\t")

colnames(Net.abr) = c("Source", "Target", "value")
colnames(Net.bac) = c("Source", "Target", "value")
colnames(Net.rel) = c("Source", "Target", "value")

Net.abr.matrix = Net.abr %>% filter(!is.na(Target), !is.na(Source)) %>% spread(Source, value, fill = 0)
Net.bac.matrix = Net.bac %>% filter(!is.na(Target), !is.na(Source)) %>% spread(Source, value, fill = 0)
Net.rel.matrix = Net.rel %>% filter(!is.na(Target), !is.na(Source)) %>% spread(Source, value, fill = 0)

tmp = Net.abr.matrix$Target
Net.abr.matrix = Net.abr.matrix[, -1]
rownames(Net.abr.matrix) = tmp

tmp = Net.bac.matrix$Target
Net.bac.matrix = Net.bac.matrix[, -1]
rownames(Net.bac.matrix) = tmp


tmp = Net.rel.matrix$Target
Net.rel.matrix = Net.rel.matrix[, -1]
rownames(Net.rel.matrix) = tmp


Net.abr.graph = graph_from_adjacency_matrix(
  as.matrix(Net.abr.matrix),
  mode = "upper",
  weighted = TRUE,
  diag = FALSE
)
Net.bac.graph = graph_from_adjacency_matrix(
  as.matrix(Net.abr.matrix),
  mode = "upper",
  weighted = TRUE,
  diag = FALSE
)
Net.rel.graph = graph_from_adjacency_matrix(
  as.matrix(Net.abr.matrix),
  mode = "upper",
  weighted = TRUE,
  diag = FALSE
)

################  Using MCL algorithm for clustering classification #########################
Net.abr.mcl = mcl(as_adjacency_matrix(Net.abr.graph),
                  addLoops = TRUE,
                  allow1 = TRUE)
Net.bac.mcl = mcl(as_adjacency_matrix(Net.bac.graph),
                  addLoops = TRUE,
                  allow1 = TRUE)
Net.rel.mcl = mcl(as_adjacency_matrix(Net.rel.graph),
                  addLoops = TRUE,
                  allow1 = TRUE)




Net.abr.cl = data_frame(
  Gene = rownames(Net.abr.matrix),
  Cluster = Net.abr.mcl$Cluster,
  DataSet = "abr"
)
Net.bac.cl = data_frame(
  Gene = rownames(Net.bac.matrix),
  Cluster = Net.bac.mcl$Cluster,
  DataSet = "bac"
)
Net.rel.cl = data_frame(
  Gene = rownames(Net.rel.matrix),
  Cluster = Net.rel.mcl$Cluster,
  DataSet = "rel"
)

Net.all.cl = bind_rows(Net.abr.cl,
                       Net.bac.cl,
                       Net.rel.cl)
############ END MCL algorithm for clustering classification   ################################

############ Clustering by Components (Calculate the maximal (weakly or strongly) connected components of a graph) #############
############ This is an alternative if MCL fails because the size of the AlleleNetwork ####################

Net.abr.comp = components(Net.abr.graph)
Net.bac.comp = components(Net.bac.graph)
Net.rel.comp = components(Net.rel.graph)

Net.abr.cl = data_frame(
  Gene = rownames(as.data.frame(Net.abr.comp$membership)),
  Cluster = Net.abr.comp$membership,
  DataSet = "abr"
)
Net.bac.cl = data_frame(
  Gene = rownames(as.data.frame(Net.bac.comp$membership)),
  Cluster = Net.bac.comp$membership,
  DataSet = "bac"
)
Net.rel.cl = data_frame(
  Gene = rownames(as.data.frame(Net.rel.comp$membership)),
  Cluster = Net.rel.comp$membership,
  DataSet = "rel"
)
Net.all.cl = bind_rows(Net.abr.cl,
                       Net.bac.cl,
                       Net.rel.cl)

########### END of Clustering by components ################################################



Connections.abr = bind_rows(Net.abr %>% select(X = Target, value),
                            Net.abr %>% select(X = Source, value)) %>% group_by(X) %>% summarise(N = sum(value))
Connections.bac = bind_rows(Net.bac %>% select(X = Target, value),
                            Net.bac %>% select(X = Source, value)) %>% group_by(X) %>% summarise(N = sum(value))
Connections.rel = bind_rows(Net.rel %>% select(X = Target, value),
                            Net.rel %>% select(X = Source, value)) %>% group_by(X) %>% summarise(N = sum(value))
Connections.all = bind_rows(Connections.abr, Connections.bac, Connections.rel)
colnames(Connections.all) = c("Gene", "Connections")


##### END AlleleNetwork ######



#### Load data #######################
lista = dir()[grep("csv", dir())]

j = 0
for (i in lista)
{
  tabla = read_delim(i, delim = "\t", col_names = FALSE)
  colnames(tabla) = c("Gene", "Reads", "RPK", "Uniq", "Coverage")
  tabla$Sample = i
  
  if (j > 0)
  {
    Full.table = bind_rows(Full.table, tabla)
  } else{
    Full.table = tabla
  }
  j = j + 1
}
#### END load data ######################




##### Join AlleleNetwork and Data abundance ############

Full.table = Full.table %>% separate(Sample, c("Sample", "DataSet", "kk"), sep =
                                       "\\.") %>% select(-kk) %>% full_join(., Net.all.cl)
Full.table = Full.table %>% left_join(., Connections.all)
Full.table$Connections[is.na(Full.table$Connections)] = 0
reads = read_delim("../Nreads.txt", col_names = TRUE, delim = ",")   #### File with the reads count per Sample

Full.table = left_join(Full.table, reads)
Full.table = Full.table %>% mutate(FinalClust = ifelse(is.na(Cluster), Gene, paste(DataSet, Cluster, sep =
                                                                                     "")))
RepresentativesOfCluster = Full.table %>% group_by(FinalClust, Gene) %>% summarise(conn = max(Connections)) %>% group_by(FinalClust) %>% mutate(top = max(conn)) %>% filter(conn == top) %>% group_by(FinalClust) %>% summarise(Representative = first(Gene))
Full.table = full_join(Full.table, RepresentativesOfCluster)

##### END Join AlleleNetwork and Data abundance ############



Full.table = Full.table %>% mutate(RPKM = RPK * 1e6 / TotalReads,
                                   UpM = Uniq * 1e6 / TotalReads)


