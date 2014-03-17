x <- matrix(1:10, 10, 10, byrow = FALSE)
x[3, 1] = NA
x[1, 10] = NaN
x[10, 10] = NA

A = .clean(x)

stopifnot(A[[1]] == 3, A(identical[[2]], as.integer(c(1, 3, 10))))
