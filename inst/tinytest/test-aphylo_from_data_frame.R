# Generating a test dataset
set.seed(1371)
x <- raphylo(20)

# Extracting the tree and annotations
tree <- x$tree

anno <- with(x, rbind(tip.annotation, node.annotation))
anno <- data.frame(id = with(tree, c(tip.label, node.label)), anno)

types <- data.frame(id = tree$node.label, x$node.type)

# Creating a aphylo tree without node types
aphylo_from_data_frame(tree, anno)

# Now including types
aphylo_from_data_frame(tree, anno, types)

# Dropping some data
aphylo_from_data_frame(tree, anno[sample.int(nrow(anno), 10),])


# Test 1: Recovering the original dataset --------------------------------------
ans0 <- aphylo_from_data_frame(tree, anno, types)
expect_equal(ans0, x)

if (interactive()) {
  
  op <- par(mfrow = c(1, 2))
  plot(x)
  plot(ans0)
  par(op)
  
}

# Test 2: Shuffling annotations ------------------------------------------------
ans0 <- aphylo_from_data_frame(tree, anno[sample(1:nrow(anno)),], types)
expect_equal(ans0, x)

# Test 3: checking errors ------------------------------------------------------
expect_error(
  aphylo_from_data_frame(tree, rbind(anno, data.frame(id=99, fun0000=0)), types),
  "present"
)

expect_error(
  aphylo_from_data_frame(tree, anno, rbind(types, data.frame(id = 99, x.node.type = 1))),
  "present"
)

tree$node.label <- NULL
expect_error(
  aphylo_from_data_frame(tree, anno, types),
  "labels"
)

# Test 4: Partial annotations
set.seed(13)
tree <- x$tree
ans0 <- aphylo_from_data_frame(tree, anno[sample.int(nrow(anno), 20),], types)

is_not_na <- which(ans0$tip.annotation[,1] != 9)
expect_identical(
  ans0$tip.annotation[is_not_na] ,
  x$tip.annotation[is_not_na]
)

