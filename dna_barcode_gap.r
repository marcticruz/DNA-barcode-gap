# Load packages    

pkgs <- c("ape", "ggplot2", "reshape2", "gdata", "utils", "spider")   

for(pkg in pkgs){ 
  library(pkg, character.only = TRUE) 
}

# Set directory

setwd("C:\\Users\\SAMSUNG\\Dropbox\\Oligo_barcode_manuscript\\analise2\\sequence_data")

# Import fasta file  

oliseq <- read.dna(file = "coi_oligo_gap3_analise2.fasta", format = "fasta")


# Split names by underline  

oliSpp <- strsplit(dimnames(oliseq)[[1]], split = "_")
oliSpp <- sapply(oliSpp, function(x) x[3])      

oliSpp

# Edit species names

oliSpp <- gsub("sp.nova", "Oligoryzomys sp.n", oliSpp, fixed = TRUE)   

oliSpp <- gsub("flavescens", "O. flavescens", oliSpp, fixed = TRUE)  

oliSpp <- gsub("utiaritensis", "O. utiaritensis", oliSpp, fixe = TRUE)  

oliSpp <- gsub("moojeni", "O. moojeni", oliSpp, fixe = TRUE)  

oliSpp <- gsub("messorius", "O. messorius", oliSpp, fixe = TRUE)  

oliSpp <- gsub("rupestris", "O. rupestris", oliSpp, fixe = TRUE)   

oliSpp <- gsub("stramineus", "O. stramineus", oliSpp, fixe = TRUE)   

oliSpp <- gsub("nigripes", "O. nigripes", oliSpp, fixe = TRUE)  

oliSpp <- gsub("fulvescens", "O. fulvescens", oliSpp, fixe = TRUE)

oliSpp <- gsub("delicatus", "O. delicatus", oliSpp, fixe = TRUE)  

oliSpp <- gsub("mattogrossae", "O. mattogrossae", oliSpp, fixe = TRUE)  

oliSpp <- gsub("microtis", "O. microtis", oliSpp, fixe = TRUE)   

oliSpp <- gsub("vegetus", "O. vegetus", oliSpp, fixe = TRUE)  

oliSpp <- gsub("spodiurus", "O. spodiurus", oliSpp, fixe = TRUE)

# Create a DNA distance matrix

row.names(oliseq) <- oliSpp

oliDist <- dist.dna(oliseq)

olimatrix <- as.matrix(oliDist)    

upperTriangle(olimatrix) <- -999
diag(olimatrix) <- -999      

# Transform the matrix into a data frame

df <- melt(as.matrix(olimatrix), varnames = c("row", "col"), header = TRUE)

df <- subset(df, value!=-999)          


# Subset the intra specific variation from the data frame

disintra <- df[df$col ==  df$row,]

# Subset the inter specific variation from the data frame

disinter <- df[df$col != df$row,]

disinter1 <- disinter[,c("col", "value")]   

colnames(disinter1)[1] <- "col"

disinter2 <- disinter[,c("row", "value")] 

colnames(disinter2)[1] <- "col"

dfdisinter <- rbind(disinter1, disinter2)  

dfdisinter

# Calculate the DNA barcode gap for each species     

names_list  <- unique(oliSpp)  

for (name in names_list){  
  
  maxdistintra <- max(disintra[disintra$col == name, 3])  
  
  mindistinter <- min(dfdisinter[dfdisinter$col == name, 2])
  
  barcodegap <- mindistinter - maxdistintra
  
  cat(name, " ", "barcode gap:", " ", barcodegap*100,"\n")
}


# Plot a histogram to represent the DNA barcode gap for each species


positio <- data.frame(h = rep(7, 14), v = rep(35, 14), gap = c("-", "0.2%", "0.8%", "1.4%", "2.3%", "2.9%", "3.4%", "3.9%", "5.2%", "6.1%", "7.2%", "7.4%", "8.3%", "8.6%"))  

ggplot(df) +
  geom_histogram(aes(df[df$row == df$col,3]), data = disintra, colour = "black", fill = "paleturquoise4", alpha = 0.3, binwidth = 0.01) +  
  geom_histogram(aes(value), data = dfdisinter, colour = "black", fill = "sandybrown", alpha = 0.3, binwidth = 0.01) +    
  facet_wrap(~col, ncol = 3) +
  annotate(geom="text", label= positio$gap, x= 0.04, y= 2400) +   
  theme(strip.text.x = element_text(size = 15, colour = "black", face = "italic")) +  
  labs(x = "Genetic distance", size = 30) + 
  labs(y = "Count", size = 30) + 
  scale_y_sqrt(breaks = c(50, 200, 450, 800, 1300, 2000, 3000))
