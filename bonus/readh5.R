library(rhdf5)

str(h5ls("/home/aliya/Liver/batch1/HG00096.h5"))

myh5 <- h5read("/home/aliya/Liver/batch1/HG00096.h5",
               name = "chr10_100009946_100009948_predictions")

head(myh5)

dim(myh5) #Should be 5313 4
