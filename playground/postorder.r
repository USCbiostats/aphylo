rm(list = ls())

# A fake tree
tree <- c("i", "a", "i", "f", "f", "b", "f", "c", "h", "i", "i", "g", "g", "d", "g", "e") #, "j", "i")
# tree <- c(tree, "i", "c")
lvls <- match(tree, letters)
tree <- matrix(lvls, ncol=2, byrow = TRUE)
set.seed(1)
plot(igraph::graph_from_edgelist(tree), layout = igraph::layout_as_tree)
tree <- matrix(lvls, ncol=2, byrow = TRUE)


list_offpring <- function(x, Data) {
  sapply(Data$I, function(y) x[which(x[,1] == y), 2])
}

list_parent <- function(x, Data) {
  sapply(Data$I, function(y) x[which(x[,2] == y), 1L])
}

# Peeling algorithm
prune <- function(i, Data, f) {
  
  if (Data$idx[i] != 0)
    return(NULL)
  
  # Peeling the offspring
  for (j in Data$Off[[i]])
    prune(j, Data, f)

  if (Data$idx[i] != 0)
    return(NULL)
  
  # Increasing counters
  Data$next_i(i)
  
  # Calling the user function
  f(i, Data)

  # Fetching the parent
  for (p in Data$Par[[i]])
    prune(p, Data, f)
  

  return(Data)
}

# Some random function to count the number of steps
# Function to recode the tree
recode_tree <- function(i, D) {
  
  # Initializing the tree
  if (!length(D$new_tree)) {
    D$new_tree <- NULL
    D$Nnode    <- which(sapply(D$Off, length) > 0)
  }
  
  if (length(D$Off[[i]]))
    D$new_tree <- rbind(
      D$new_tree,
      cbind(D$count, D$idx[D$Off[[i]]])
    )
 
}

# Step 1, listing offspring and parent
# n <- 500
# tree <- ape::rtree(n)$edge

#' Creates a peeling object
#' @param tree An integer matrix with two columns. The edgelist
#' @return An environemt to be passed to `prune`:
#' - `off` Current idnumber
newprune <- function(tree) {
  
  n <- max(tree)
  
  # Creating variables
  Data       <- new.env()
  Data$count <- 0L
  Data$I     <- 1L:n
  Data$idx  <- integer(n)
  
  lockBinding("count", Data)
  lockBinding("I", Data)
  lockBinding("idx", Data)
  
  # This function increases the counter and sets the id in the corresponding
  # slot. This is for internal use only
  Data$next_i <- function(i) {
    unlockBinding("count", Data)
    unlockBinding("idx", Data)
    
    Data$count   <- Data$count + 1L
    Data$idx[i] <- Data$count 
    
    lockBinding("count", Data)
    lockBinding("idx", Data)
  }
  lockBinding("next_i", Data)
  
  # Listing offpring and parent
  Data$Off   <- list_offpring(tree, Data)
  Data$Par   <- list_parent(tree, Data)
  
  lockBinding("Off", Data)
  lockBinding("Par", Data)
  
  Data
  
}

ans <- microbenchmark::microbenchmark(
  prune(1, newprune(tree), recode_tree),
  unit = "ms", times = 1000)


ans <- newprune(tree)
prune(8, ans, recode_tree)

oldpar <- par(no.readonly = TRUE)
par(mfrow=c(1, 2))
set.seed(1);plot(igraph::graph_from_edgelist(tree), layout = igraph::layout_as_tree)
set.seed(1);plot(igraph::graph_from_edgelist(ans$new_tree), layout = igraph::layout_as_tree)
par(oldpar)
