################
#	
#	R script with the pipeline to perform the Allele Network (alternatively to perl script)
#	and perform the clustering analysis to create the Mapping Gene Clusters. Moreover, the
#	script includes the loading process of count data and the fusion of both parts (allele network and count data)
#
#
################


library("tidyverse")
library("MCL")
library("igraph")



##### AlleleNetwork Alternative ###############################################
lista = dir()[grep("abr.drs",dir())]

Net = data.frame(Gen = "",siguiente = "",N =0)
for(i in lista)
{
  tabla = read_delim(i, delim = "\t", col_names = FALSE)
  colnames(tabla) = c("Read","Gen","X1","X2","X3")
  Net = bind_rows(Net,tabla %>% mutate(siguienteR = lead(Read), siguienteG = (lead(Gen))) %>% mutate(link = (Read == siguienteR)) %>% filter(link == TRUE) %>% group_by(Gen,siguienteG) %>% summarise(N = n()))
}

Net.abr = Net %>% group_by(Gen,siguienteG) %>% summarise(N = sum(N))

lista = dir()[grep("bac.drs",dir())]

Net = data.frame(Gen = "",siguiente = "",N =0)
for(i in lista)
{
  tabla = read_delim(i, delim = "\t", col_names = FALSE)
  colnames(tabla) = c("Read","Gen","X1","X2","X3")
  Net = bind_rows(Net,tabla %>% mutate(siguienteR = lead(Read), siguienteG = (lead(Gen))) %>% mutate(link = (Read == siguienteR)) %>% filter(link == TRUE) %>% group_by(Gen,siguienteG) %>% summarise(N = n()))
}

Net.bac = Net %>% group_by(Gen,siguienteG) %>% summarise(N = sum(N))

lista = dir()[grep("rel.drs",dir())]

Net = data.frame(Gen = "",siguiente = "",N =0)
for(i in lista)
{
  tabla = read_delim(i, delim = "\t", col_names = FALSE)
  colnames(tabla) = c("Read","Gen","X1","X2","X3")
  Net = bind_rows(Net,tabla %>% mutate(siguienteR = lead(Read), siguienteG = (lead(Gen))) %>% mutate(link = (Read == siguienteR)) %>% filter(link == TRUE) %>% group_by(Gen,siguienteG) %>% summarise(N = n()))
}

Net.rel = Net %>% group_by(Gen,siguienteG) %>% summarise(N = sum(N))


Net.abr.graph = graph_from_edgelist(as.matrix(Net.abr %>% filter(N>0) %>% select(Gen,siguienteG)), directed = FALSE)
Net.bac.graph = graph_from_edgelist(as.matrix(Net.bac %>% filter(N>0) %>% select(Gen,siguienteG)), directed = FALSE)
Net.rel.graph = graph_from_edgelist(as.matrix(Net.rel %>% filter(N>0) %>% select(Gen,siguienteG)), directed = FALSE)

Net.abr.matrix = as_adjacency_matrix(Net.abr.graph)
Net.abr.mcl = mcl(Net.abr.matrix, addLoops = TRUE, allow1 = TRUE)

Net.bac.matrix = as_adjacency_matrix(Net.bac.graph)
Net.bac.mcl = mcl(Net.bac.matrix, addLoops = TRUE, allow1 = TRUE)

Net.rel.matrix = as_adjacency_matrix(Net.rel.graph)
Net.rel.mcl = mcl(Net.rel.matrix, addLoops = TRUE, allow1 = TRUE)


Net.abr.cl = data_frame(Gene = rownames(Net.abr.matrix),Cluster = Net.abr.mcl$Cluster, DataSet = "abr")
Net.bac.cl = data_frame(Gene = rownames(Net.bac.matrix),Cluster = Net.bac.mcl$Cluster, DataSet = "bac")
Net.rel.cl = data_frame(Gene = rownames(Net.rel.matrix),Cluster = Net.rel.mcl$Cluster, DataSet = "rel")

tmp = bind_rows(Net.abr.cl,
                Net.bac.cl,
                Net.rel.cl)

Connections.abr = bind_rows(Net.abr %>% select(X = Gen,N), Net.abr %>% select(X = siguienteG,N)) %>% group_by(X) %>% summarise(N = sum(N))
Connections.bac = bind_rows(Net.bac %>% select(X = Gen,N), Net.bac %>% select(X = siguienteG,N)) %>% group_by(X) %>% summarise(N = sum(N))
Connections.rel = bind_rows(Net.rel %>% select(X = Gen,N), Net.rel %>% select<Q(X = siguienteG,N)) %>% group_by(X) %>% summarise(N = sum(N))
Connections.all = bind_rows(Connections.abr, Connections.bac, Connections.rel)


##### END AlleleNetwork ######



#### Load data #######################
lista = dir()[grep("csv",dir())]

j=0
for(i in lista)
{
  tabla = read_delim(i, delim = "\t", col_names = FALSE)
  colnames(tabla) = c("Gene","Reads","RPK","Uniq","Coverage")
  tabla$Sample = i
  
  if(j > 0)
  {
    Full.table = bind_rows(Full.table,tabla)
  }else{
    Full.table = tabla
  }
  j=j+1
}
#### END load data ######################




##### Join AlleleNetwork and Data abundance ############

Full.table = Full.table %>% separate(Sample, c("Sample","DataSet","kk"), sep ="\\.") %>% select(-kk) %>% full_join(.,tmp)

#### Modify to fix the filename of count.reads to your file that contains the number of total reads per sample
reads = read_table("count.reads",col_names = FALSE)
colnames(reads) = c("TotalReads","Sample") 
####

### optional sentences
reads$Sample = gsub(".R1.fastq","",reads$Sample) 
reads$TotalReads = reads$TotalReads/4
###

Full.table = full_join(Full.table, reads)
Full.table = Full.table %>% mutate(FinalClust = ifelse(is.na(Cluster),Gene,paste(DataSet,Cluster,sep="")))
RepresentativesOfCluster = Full.table %>% group_by(FinalClust,Gene) %>% summarise(conn = max(Connections)) %>% group_by(FinalClust) %>% mutate(top = max(conn)) %>% filter(conn == top) %>% group_by(FinalClust) %>% summarise(Representative = first(Gene))

##### END Join AlleleNetwork and Data abundance ############






