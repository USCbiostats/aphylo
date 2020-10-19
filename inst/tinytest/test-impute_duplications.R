# Alternative definition of the imputation function, this is slower
# by definition. But it should be good enough to impute duplication.
get_species <- function(i, offspring, species, ans = NULL) {
  
  # Leaf node?
  if (length(offspring[[i]]) == 0) {
    
    ans <- c(ans, species[i])
    
  } else {
    
    for (j in offspring[[i]])
      ans <- get_species(j, offspring, species, ans)
    
  }
  
  return(ans)
  
}

imputate_duplications_alt <- function(tree, species = NULL) {
  
  # Checking the input
  if (!inherits(tree, "phylo"))
    stop("This -tree- should be of class 'phylo'", call. = FALSE)
  
  # Getting the species
  if (!length(tree$tip.label))
    stop("This tree has no labels on its tips.", call. = FALSE)
  
  if (!length(species))
    species <- gsub(pattern = ".+[_]", replacement = "", tree$tip.label)
  
  # Listing offspring
  off   <- aphylo::list_offspring(tree)
  size_ <- ape::Nnode(tree, internal.only = FALSE)
  dpl   <- logical(size_)
  for (p in (ape::Ntip(tree) + 1):size_) {
    
    # Collecting the species
    species_p <- vector("list", length(off[[p]]))
    
    for (i in seq_along(off[[p]])) {
      species_p[[i]] <- get_species(
        i         = off[[p]][i],
        offspring = aphylo::list_offspring(tree),
        species   = species
      )
    }
    
    # Comparing lists
    for (i in combn(1:length(species_p), 2, simplify = FALSE)) {
      if (any(species_p[[i[1]]] %in% species_p[[i[2]]])) {
        dpl[p] <- TRUE   
        break
      }
      
    }
  }
  
  return(dpl)
  
  
}

# Generating the data
set.seed(1231)
tree2 <- aphylo::sim_tree(40)
s     <- sample.int(10, 40, replace=TRUE)

dpl0 <- imputate_duplications(tree2, species = s)
dpl1 <- imputate_duplications_alt(tree2, species = s)

expect_equal(dpl0, dpl1)

