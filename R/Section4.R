library("Biostrings")
library("dplyr")
found_genes <- data.frame(matrix(nrow=5, ncol=7)) 
colnames(found_genes) <- c('NP', 'VP35', 'VP40', 'GP', 'VP30', 'VP24', 'L')
rownames(found_genes) <- c('reston', 'taiforest', 'sudan', 'zaire', 'bundibugyo')

reston <- readDNAStringSet('../Data/Reston_genome.fasta')
reston_names <- names(reston)
reston_genes <- paste(reston)
reston <- data.frame(reston_names, reston_genes)

marburg <- readDNAStringSet('../Data/Marburg_genes.fasta')
marburg_names <- names(marburg)
marburg_genes <- paste(marburg)
marburg <- data.frame(marburg_names, marburg_genes)

zaire <- readDNAStringSet('../Data/marburg.fasta')
zaire_names <- names(zaire)
zaire_genes <- paste(zaire)
zaire <- data.frame(marburg_names, zaire_genes)

bundibugyo <- readDNAStringSet('../Data/marburg.fasta')
bundibugyo_names <- names(bundibugyo)
bundibugyo_genes <- paste(bundibugyo)
bundibugyo <- data.frame(marburg_names, bundibugyo_genes)

sudan <- readDNAStringSet('../Data/marburg.fasta')
sudan_names <- names(sudan)
sudan_genes <- paste(sudan)
sudan <- data.frame(marburg_names, sudan_genes)

taiforest <- readDNAStringSet('../Data/marburg.fasta')
taiforest_names <- names(taiforest)
taiforest_genes <- paste(taiforest)
taiforest <- data.frame(marburg_names, taiforest_genes)

intervals <- c(0)
for (gene_id in 1:5){
  for(i in 1:7){
    alignment <- pairwiseAlignment(pattern = marburg$marburg[i], subject = reston, type='global-local' )
    found_genes[gene_id, i] <- (reston %>% substr(alignment %>% subject() %>% start(), alignment %>% subject() %>% end()) %>% as.character() )
    if (i == 1){
      intervals <- intervals %>% rbind( alignment %>% subject() %>% start() ) 
      intervals <- intervals %>% rbind( alignment %>% subject() %>% end() )
    }
  }
}

genomes <- c(reston, taiforest, sudan, zaire, bundibugyo)
global_distance <- matrix(0, ncol = 5, nrow = 5)
colnames(global_distance) <- c('reston', 'taiforest', 'sudan', 'zaire', 'bundibugyo')
rownames(global_distance) <- c('reston', 'taiforest', 'sudan', 'zaire', 'bundibugyo')

for(i in 1:5){
  for(j in 1:5){
    global_distance[i, j] <- pairwiseAlignment(pattern = genomes[i], subject = genomes[j], type = 'global') %>% nedit()
  }
}


time_distance <- matrix(0, ncol = 5, nrow = 5)
colnames(time_distance) <- c('reston', 'taiforest', 'sudan', 'zaire', 'bundibugyo')
rownames(time_distance) <- c('reston', 'taiforest', 'sudan', 'zaire', 'bundibugyo')

for(i in 1:5){
  for(j in 1:5){
    alignment <- global_distance[i, j] / ( ( nchar(genomes[i]) + nchar(genomes[j]) ) / 2 )
    time_distance[i, j] <- ( (-3/4) * log(1 - (4/3) * alignment ) ) / (1.9 * 0.001)
  }
}

time_distance %>% View()


intervals <- c(0)

for(i in 1:7){
  alignment <- pairwiseAlignment(pattern = marburg[i], subject = reston, type='global-local' )
  found_genes[1, i] <- (reston %>% substr(alignment %>% subject() %>% start(), alignment %>% subject() %>% end()) %>% as.character() )
  intervals <- intervals %>% rbind( alignment %>% subject() %>% start() ) 
  intervals <- intervals %>% rbind( alignment %>% subject() %>% end() )
}
