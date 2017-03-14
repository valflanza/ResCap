library("tidyverse")
library("MCL")
library("igraph")




##### AlleleNetwork Alternative ###############################################



Net.abr.graph = graph_from_edgelist(as.matrix(Net.abr), directed = FALSE)
Net.bac.graph = graph_from_edgelist(as.matrix(Net.bac), directed = FALSE)
Net.rel.graph = graph_from_edgelist(as.matrix(Net.rel), directed = FALSE)

Net.abr.matrix = as_adjacency_matrix(Net.abr.graph)
Net.abr.mcl = mcl(Net.abr.matrix, addLoops = TRUE, allow1 = TRUE)

Net.bac.matrix = as_adjacency_matrix(Net.bac.graph)
Net.bac.mcl = mcl(Net.bac.matrix, addLoops = TRUE, allow1 = TRUE)

Net.rel.matrix = as_adjacency_matrix(Net.rel.graph)
Net.rel.mcl = mcl(Net.rel.matrix, addLoops = TRUE, allow1 = TRUE)


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

Connections.abr = bind_rows(Net.abr %>% select(X = Gen, N), Net.abr %>% select(X = siguienteG, N)) %>% group_by(X) %>% summarise(N = sum(N))
Connections.bac = bind_rows(Net.bac %>% select(X = Gen, N), Net.bac %>% select(X = siguienteG, N)) %>% group_by(X) %>% summarise(N = sum(N))
Connections.rel = bind_rows(Net.rel %>% select(X = Gen, N), Net.rel %>% select(X = siguienteG, N)) %>% group_by(X) %>% summarise(N = sum(N))
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
reads = read_table("../count.reads",col_names = FALSE)
colnames(reads) = c("TotalReads","Sample")
reads$Sample = gsub(".R1.fastq","",reads$Sample)
reads$TotalReads = reads$TotalReads/4
Full.table = left_join(Full.table, Nreads)
Full.table = Full.table %>% mutate(FinalClust = ifelse(is.na(Cluster), Gene, paste(DataSet, Cluster, sep =
                                                                                     "")))
RepresentativesOfCluster = Full.table %>% group_by(FinalClust, Gene) %>% summarise(conn = max(Connections)) %>% group_by(FinalClust) %>% mutate(top = max(conn)) %>% filter(conn == top) %>% group_by(FinalClust) %>% summarise(Representative = first(Gene))
Full.table = full_join(Full.table, RepresentativesOfCluster)
##### END Join AlleleNetwork and Data abundance ############
