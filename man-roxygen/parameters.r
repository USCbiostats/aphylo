#' \code{}
#' <%=ifelse(exists(".annotationslab"), "@param annotations Matrix of size \\eqn{n\\times P + 1}{n*P + 1}. Annotation data, where the first column is the nodeid (should coincide with the coding of the tree), the following columns are function annotations coded as \\code{1} if the node has the function, \\code{0} if the node doesn't has the function and \\code{9} or \\code{NA} otherwise (no information).", "") %>
#' <%=ifelse(exists(".annotations"),    "@param annotations Matrix of size \\eqn{n\\times P}{n*P}. Annotation data, each column is a functional annotation coded as \\code{1} if the node has the function, \\code{0} if the node doesn't has the function and \\code{9} or \\code{NA} otherwise (no information).", "") %>
#' <%=ifelse(exists(".tip.annotation"), "@param tip.annotation,node.annotation Annotation data. See [aphylo-class].", "") %>
#' <%=ifelse(exists(".psi"),         "@param psi Numeric vector of length 2. Misclasification probabilities. (see \\code{\\link{LogLike}}).", "") %>
#' <%=ifelse(exists(".mu"),          "@param mu_d,mu_s Numeric vector of length 2. Gain/loss probabilities (see \\code{\\link{LogLike}}).", "") %>
#' <%=ifelse(exists(".eta"),         "@param eta Numeric vector of length 2. Annotation bias probabilities (see \\code{\\link{LogLike}}).", "") %>
#' <%=ifelse(exists(".Pi"),          "@param Pi Numeric scalar. Root node probability of having the function (see \\code{\\link{LogLike}}).", "") %>
#' <%=ifelse(exists(".offspring"),   "@param offspring List of length \\eqn{n}. Offspring of each node (see \\code{\\link{new_aphylo}}).", "") %>
#' <%=ifelse(exists(".S"),           "@param S Integer matrix of size \\eqn{2^P\\times P}{2^P * P}. See \\code{\\link{states}}.", "") %>
#' <%=ifelse(exists(".tree"),       "@param tree An object of class [phylo][ape::read.tree]", "") %>
#' <%=ifelse(exists(".types"),       "@param tip.type,node.type Integer vectors with values {0,1}. 0 denotes duplication node and 1 speciation node. This is used in [LogLike].", "") %>
NULL
