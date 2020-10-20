imput_species <- function(tree, i = 1L, env = NULL, species = NULL) {
  
  if (is.null(env)) {
    
    if (is.null(species))
      stop("The parameter -species- should be specified.", call. = FALSE)
    
    if (!inherits(tree, "phylo"))
      stop("-tree- must be a phylo class object.", call. = FALSE)
    
    env <- list()
    env$species     <- species
    env$offspring   <- aphylo::list_offspring(tree)
    env$checked     <- rep(FALSE, ape::Nnode(tree, internal.only = FALSE))
    env$off_species <- vector("list", length(env$checked))
    env$ntip        <- ape::Ntip(tree)
    env$ntot        <- length(env$checked)
    env$par         <- rbind(tree$edge, c(NA, env$ntip + 1))
    env$par         <- env$par[order(env$par[,2]), 1]
    
    env <- list2env(env)
  }
  
  if (env$checked[i]) {
    return(env)
  }
  env$checked[i] <- TRUE
  
  # Going down
  if (!length(env$offspring[[i]])) { # Is a leaf
    env$off_species[[i]] <- env$species[i]
  } else { # Has offspring
    sapply(env$offspring[[i]], imput_species, tree = tree, env = env)
    # Collecting species
    env$off_species[[i]] <- unique(unlist(env$off_species[env$offspring[[i]]]))
    
  }
  
  # Root node
  if (i == (env$ntip + 1))
    return(as.list(env))
  
  if (is.na(env$par[i])) {
    
    warning("This is fishy")
    
  }
  
  imput_species(tree = tree, i = env$par[i], env = env)
  
}

#' Impute duplication events based on a vector of species
#' 
#' Uses a simple algorithm to impute duplication events based on the
#' terminal genes of the tree. An interior node is a duplication event
#' if a specie has two or more leafs within its clade.
#' 
#' @param tree An object of class [ape::phylo].
#' @param species A character vector of length `ape::Ntip(tree)` (see details).
#' @details 
#' This function will take a vector of species and, based on that, assign
#' duplication events throughout the interior nodes. An interior node is labeled
#' as a duplication event if two or more of the leaves within it are from the
#' same species.
#' @return A logical vector of length `ape::Nnode(tree, internal.only = FALSE)`
#' with `TRUE` to indicate that the corresponding node is a duplication event.
#' The order matches that in the input tree. 
#' @examples 
#' 
#' # Data from PANTHER
#' path <- system.file("tree.tree", package="aphylo")
#' ptree <- read_panther(path)
#' 
#' # Extracting the species
#' sp <- gsub(".+[:]|[|].+", "" , ptree$tree$tip.label)
#' 
#' # Imputing duplications
#' imputate_duplications(ptree$tree, species = sp)
#' @export
imputate_duplications <- function(tree, species) {
  
  # Checking the input
  if (!inherits(tree, "phylo"))
    stop("This -tree- should be of class 'phylo'", call. = FALSE)
  
  if (length(species) != ape::Ntip(tree))
    stop(
      "The -species- parameter does not match the number of tips in the tree.",
      call. = FALSE
      )
  
  # Retrieving the species
  DAT <- imput_species(tree = tree, species = species)
  
  # Vector of answers
  dpl   <- rep(FALSE, length(DAT$par))
  
  for (p in (DAT$ntip + 1):DAT$ntot) {
    
    # Getting the species in p
    species_p <- DAT$off_species[DAT$offspring[[p]]]
    
    # Comparing lists
    for (i in utils::combn(1:length(species_p), 2, simplify = FALSE)) {
      if (any(species_p[[i[1]]] %in% species_p[[i[2]]])) {
        dpl[p] <- TRUE
        break
      }
    }
    
  }
  
  return(dpl)
}

