setwd("~/Desktop/Complex-and-Social-Networks/Lab2-Degree-Distribution")

write_summary <- function(language,file) {
   degree_sequence = read.table(file, header = FALSE)
   cat(language,length(degree_sequence$V1),max(degree_sequence$V1),sum(degree_sequence$V1)/length(degree_sequence$V1),length(degree_sequence$V1)/sum(degree_sequence$V1),"\n")
}

source = read.table("list.txt", 
         header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
         as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
        )
for (x in 1:nrow(source)) {
    write_summary(source$language[x], source$file[x])
}


write_AICs <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  cat(language,AICs(degree_sequence$V1),"\n")
}

for (x in 1:nrow(source)) {
  write_AICs(source$language[x], source$file[x])
}



output <- matrix(nrow=length(source$language),ncol=5)

write_AICs <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  
  #cat(language,AICs(degree_sequence$V1),"\n")
  return(AICs(degree_sequence$V1))
}

for (x in 1:nrow(source)) {
  output[x,] <- write_AICs(source$language[x], source$file[x])
}
dimnames(output) <- list(source$language, c("zeta","zeta_2","RT_zeta","geom","poisson"))
output
library(xtable)
xtable(output)

