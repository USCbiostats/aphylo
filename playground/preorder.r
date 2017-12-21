rm(list = ls())

list_offpring <- function(x, Data) {
  sapply(Data$I, function(y) x[which(x[,1] == y), 2])
}

list_parent <- function(x, Data) {
  sapply(Data$I, function(y) x[which(x[,2] == y), 1L])
}

newprune <- function(tree) {
  
  # Creating variables
  Dat       <- new.env()
  Dat$n     <- max(tree)
  Dat$tree  <- tree
  Dat$count <- 0L
  Dat$I     <- 1L:Dat$n
  Dat$idx   <- integer(Dat$n)
  
  Dat$Off   <- list_offpring(tree, Dat)
  Dat$Par   <- list_parent(tree, Dat)
  
  # Methods
  Dat$nexti <- function(i) {
    Dat$count  <- Dat$count + 1L
    Dat$idx[i] <- Dat$count
  }

  Dat
}

preorder_prune <- function(i, Data, f = function(i, d) {1+1}) {
  
  if (Data$idx[i] != 0)
    return(NULL)
  
  # Take it to the dad
  for (j in Data$Par[[i]]) {
    preorder_prune(j, Data, f)
  }
  
  # Increasing counters
  if (Data$idx[i] != 0)
    return(NULL)
  
  Data$nexti(i)
  
  # Calling the user function
  f(i, Data)
  
  # Fetching the offspring
  for (j in Data$Off[[i]]) 
    preorder_prune(j, Data, f)
  
  
  return(Data)
  
}


# A fake tree
tree <- c("i", "a", "i", "f", "f", "b", "f", "c", "h", "i", "i", "g", "g", "d", "g", "e") #, "j", "i")
# tree <- c(tree, "i", "c")
lvls <- match(tree, letters)
tree <- matrix(lvls, ncol=2, byrow = TRUE)
tree <- matrix(lvls, ncol=2, byrow = TRUE)

# Function to recode the tree
recode_tree <- function(i, D) {
  
  # Initializing the tree
  if (!length(D$new_tree)) {
    D$new_tree <- NULL
    D$Nnode    <- which(sapply(D$Off, length) > 0)
  }
    
  if (length(D$Par[[i]]))
    D$new_tree <- rbind(
      D$new_tree,
      cbind(D$idx[D$Par[[i]]], D$count)
      )
  
}

edgelist_to_ape <- function(x) {
  structure(
    list(
      edge  = x,
      Nnode = 
    )
  )
}

ans <- newprune(tree)
preorder_prune(1, ans, recode_tree)

oldpar <- par(no.readonly = TRUE)
par(mfrow=c(1, 2))
set.seed(1);plot(igraph::graph_from_edgelist(tree), layout = igraph::layout_as_tree)
set.seed(1);plot(igraph::graph_from_edgelist(ans$new_tree), layout = igraph::layout_as_tree)
par(oldpar)
