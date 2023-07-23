
profile_ll_B_cup <- function(X,
                             Y,
                             B_cup){
J <- ncol(Y)
p <- ncol(X)
B <- B_from_B_cup(B_cup, J = J, p = p)

return(profile_ll(X,Y,B))

}
