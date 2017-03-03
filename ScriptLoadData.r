#####
#
#	R script to load count data using tidyverse package.
#	All count files must be in the R working directory
#
####
library("tidyverse")
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

Full.table = Full.table %>% separate(Sample, c("Sample","DataSet","kk"), sep ="\\.") %>% select(-kk)
