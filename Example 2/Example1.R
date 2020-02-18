# Squaring elements of a given vector

square_for <- function(x){
  # Use the for loop
  y <- x
  for(i in 1 : length(x)) {
    y[i] = x[i]^2
  }
  return(y)
}

square_sapply <- function(x){
  # Use the sapply function
  y <- sapply(x, function(t) t^2)
  return(y)
}

square_vec <- function(x){
  # Use power function in vector form
  y <- x^2
  return(y)
}

# Create a vector x of size 100,000 of normal variables

x <- rnorm(100000)

# Verify that all functions return the same output

y1 <- square_for(x)
y2 <- square_sapply(x)
y3 <- square_vec(x)

sum(abs(y1-y2))
sum(abs(y1-y3))
#sum(abs(y1-y4))

# Use microbenchmark package to compare three functions in terms of speed

library(microbenchmark)

microbenchmark(
  square_for(x),
  square_sapply(x),
  square_vec(x)
)
