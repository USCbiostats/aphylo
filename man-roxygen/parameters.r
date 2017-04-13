#' \code{}
#' 
#' <%=ifelse(exists("annotations"), "@param annotations Matrix of size \\eqn{n\\times P}{n*P}. Annotation data, one column per function coded as \\code{1} if the node has the function, \\code{0} if the node doesn't has the function and \\code{9} otherwise (no information).", "") %>
#' <%=ifelse(exists("psi"),         "@param psi Numeric vector of length 2. Misclasification probabilities. (see \\code{\\link{LogLike}}).", "") %>
#' <%=ifelse(exists("mu"),          "@param mu Numeric vector of length 2. Gain/loss probabilities (see \\code{\\link{LogLike}}).", "") %>
#' <%=ifelse(exists("Pi"),          "@param Pi Numeric vector of length 2. Root node probabilities (see \\code{\\link{LogLike}}).", "") %>
#' <%=ifelse(exists("noffspring"),  "@param noffspring Integer vector of length \\eqn{n}. Number of offspring each node has (see \\code{\\link{new_aphylo}}).", "") %>
#' <%=ifelse(exists("offspring"),   "@param offspring List of length \\eqn{n}. Offspring of each node (see \\code{\\link{new_aphylo}}).", "") %>
#' <%=ifelse(exists("S"),           "@param S Integer matrix of size \\eqn{2^P\\times P}{2^P * P}. See \\code{\\link{states}}.", "") %>
#' <%=ifelse(exists("edges"),       "@param edges Integer matrix with two columns. An edgelist where each row is in the form of \\code{(parent, offspring)}.", "") %>
