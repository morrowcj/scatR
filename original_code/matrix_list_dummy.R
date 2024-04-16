Ms <- rep("M", 5)

x = paste(Ms, 1:5, sep = "")

newList = list()
for (i in 1:length(x)){
  newList[[i]] = matrix(c(x[i], 0, 0, x[i]), ncol = 2, byrow = F)
}

newList
