# Exercise 1
# Create a directory "exercisesR". Inside the folder create another
# directory named "Es" and set it as working directory for R
# ------------------------------------------------------------------------------------------------------

# Exercise 2
# Open a help window for the mean function
# ------------------------------------------------------------------------------------------------------

# Exercise 3
# Execute the operation sqrt(2*10 + 5), put the result in the variable x and print it to the screen.
# ------------------------------------------------------------------------------------------------------

# Exercise 4
# Define a vector v1 = (1,2,3,4,5), v2 = (23,13,1,33,2). Calculate the vector v3 = v1 + v2. Execute the 
# log10 of each element of v3 and put the result in v4. Print the square of each element of v4 
# ------------------------------------------------------------------------------------------------------

# Exercise 5
# Print the sum of the squares of the first ten numbers using sequences. 
# Hint: the function sum(v) calculate the sum of the elements of the vector v.
# ------------------------------------------------------------------------------------------------------

# Exercise 6
# Print the sum of the first ten even numbers.
# ------------------------------------------------------------------------------------------------------

# Exercise 7
# Create a vector v with ten numbers in the range (1,100). Print how many numbers are
# greater than 50.
# Hint: Use the function sum.
# ------------------------------------------------------------------------------------------------------

# Exercise 8
# Define the following vector: v <- c(NA, NA, NA, NA, 1:100). Create a vector v2 from v with elements greater than 50
# but lower than 80. Exclude NA values.
# ------------------------------------------------------------------------------------------------------

# Exercise 9
# Define a vector v with elements from 1 to 20. Create a vector v2 with the odd elements of v.
# Assign to v2 a names vector with the letters of even position (from variable LETTERS). Print v2. 
# Use the negative indexing.
# ------------------------------------------------------------------------------------------------------

# Exercise 10
# Build a factor of 25 elements (strings) formed by 3 different levels
# ------------------------------------------------------------------------------------------------------

# Exercise 11
# Which data structures could be used to model a data set with a set of patients
# healthy and sick, subjected to 5 different types of clinical analysis?
# ------------------------------------------------------------------------------------------------------

# Exercise 12
# Create a list l with the following elements:
#   - a vector named seq with the first 10 numbers
#   - a vector named letters with the first 10 letters from LETTERS
#   - a vector named bool with some boolean values

# Print the length of l
# Print the elements of the vector "seq" using its name and its numerical index
# Print the second element of vector letters
# Add to l, as last element, the string "last"
# ------------------------------------------------------------------------------------------------------

# Exercise 13
# Given the list li <- list(m = matrix(rnorm(64), nrow = 8), s = c(rep("T", 3), rep("G", 5))) 
# - extract the second column of the matrix
# - the "G"s of the vector s
# Add to the list a vector of 10 random numbers.
# Access in 3 different ways to the second component of the list
# What is the difference between li[1] and li[[1]]?
# With the function diag, print the diagonal of the matrix of li 
# ------------------------------------------------------------------------------------------------------

# Exercise 14
# Define a 5x5 matrix m formed by the numbers from 1 to 25. Assign names to columns and rows with the first 5 letters. 
# Create a new matrix m2, that is m without column 1 and 3 but a last column with all one named 'O'.
# ------------------------------------------------------------------------------------------------------

# Exercise 15
# Define a 5x5 matrix m with all one
# Define a 5x2 matrix m2 with the first 10 numbers.
# Put in m3 the result of matrix multiplication m*m2.
# Use %*% to multiply matrices
# Put m, m2 and m3 in a list l
# Print l
# ------------------------------------------------------------------------------------------------------

# Exercise 16
# Create a 4x4 matrix m with the first 16 letters. Fill the matrix by rows.  
# Define a matrix m2 with all columns of m except the third
# Define a vector v with 4 arbitrary strings
# Create a matrix m3 which is the concatenation of m and v by columns
# ------------------------------------------------------------------------------------------------------

# Exercise 17
# Build 2 matrices M and N with 5 columns. Build (if possible) with rbind a matrix with
# 5 columns and as rows, the rows of M and N. Is it possible build a matrix from M and N
# using cbind?
# ------------------------------------------------------------------------------------------------------

# Exercise 18
# Given 2 square matrices A and B (3x3) formed by random numbers, compute their element-wise
# product and their rows times columns product. Define a vector b with 3 elements and solve
# (function "solve") the linear system Ax = b.
# ------------------------------------------------------------------------------------------------------

# Exercise 19
# Build a 10x10 matrix composed by random numbers using at least two different ways.
# Use the functions matrix and array
# ------------------------------------------------------------------------------------------------------

# Exercise 20
# Build 2 matrices x (3x4) and y (5x3) of characters. Change the matrix x with a unique
# assignment in order to replace its first two columns with the last two rows of y
# ------------------------------------------------------------------------------------------------------

# Exercise 21
# Build a 3x4x5 array formed by the first 60 positive natural numbers. Multiply a with
# an array with the same dimensions caracterized by random numbers. Is it possible
# multiply "a" with a vector with 60 numbers? If yes, what do you get?
# Modify "a" in order to obtain an array 4-D of dimensions 3x5x2x2. 
# ------------------------------------------------------------------------------------------------------

# Exercise 22
# Create a data frame df from an arbitrary 5x5 matrix.
# Add to df a factor with levels "1" and "0".
# ------------------------------------------------------------------------------------------------------

# Exercise 23
# Create a data frame df using three vectors and 1 factor (lengths > 5)
# Create a data frame df2 that is a subset of df binded with a 4x4 matrix
# ------------------------------------------------------------------------------------------------------

# Exercise 24
# Create a data frame df using a list containing 2 matrices with a 
# different number of columns  
# ------------------------------------------------------------------------------------------------------

# Exercise 25
# Create a data frame containing the data frame of the previous exercise plus a vector
# Extract one of the columns in three ways
# ------------------------------------------------------------------------------------------------------

# Exercise 26
# Print the following attributes of the iris dataset:
# - class
# - dimensions
# - number of rows
# - number of columns
# - space occupied by the dataset
# - names attribute
# ------------------------------------------------------------------------------------------------------

# Exercise 27 - Summary - Objects in R
# Define a vector v with the following numbers: 27, 22, 23, 19, 21, 30 
# Print v
# Print the dimensions of v
# Print the length of v
# Print the class of v
# Assign to v the dimension (2,3)
# Print v
# Print the dimensions of v
# Print the class of v
# Assign to v the dimension NULL
# Define a factor s with the values: "M", "F", "M", "F", "F", "F"
# Print s
# Define a vector p with the strings: "Marco", "Monica", "Simone", "Giuliana", "Maria", "Marta"
# Print the class of p
# Use cbind to create a matrix M with columns p, v and s
# Print the class of each column of M
# Define a list l with the elements p, v and s
# Assign to l the names "Name", "Age" and "Gender"
# Print l
# Create a data frame using the list l
# Print the class of df
# Print df
# Print the class of each column of df
# ------------------------------------------------------------------------------------------------------

# Exercise 28
# Define a vector v with 20 numbers between 1 and 100
# Using a for loop, iterate over v and put in s the sum of the numbers greater than 50
# ------------------------------------------------------------------------------------------------------

# Exercise 29
# Define a 10x10 matrix m
# Iterate over m and set each element m[i,j] based on these rules:
# - if i != j then  m[i,j] = i*j
# - otherwise       m[i,j] = 1
# ------------------------------------------------------------------------------------------------------

# Exercise 30
# Define a 5x5 matrix m1 with numbers between 1 and 10 and a 5x2 matrix m2 with numbers from 1 to 10
# Using for loops calculate the matrix multiplication m1*m2. Put the result in m3
# ------------------------------------------------------------------------------------------------------

# Exercise 31
# You're about to write your first function! Just like you would assign a value 
# to a variable with the assignment operator, you assign functions in the following
# way:
#
# function_name <- function(arg1, arg2){
#	# Manipulate arguments in some way
#	# Return a value
# }
#
# The "variable name" you assign will become the name of your function. arg1 and
# arg2 represent the arguments of your function. You can manipulate the arguments
# you specify within the function. After sourcing the function, you can use the 
# function by typing:
# 
# function_name(value1, value2)
#
# Below we will create a function called boring_function. This function takes
# the argument `x` as input, and returns the value of x without modifying it.
# Delete the pound sign in front of the x to make the function work! Be sure to 
# save this script and type submit() in the console after you make your changes.

boring_function <- function(x) {
  x
}

# ------------------------------------------------------------------------------------------------------

# Exercise 32
#Implement the function my_mean
# You're free to implement the function my_mean however you want, as long as it
# returns the average of all of the numbers in `my_vector`.
#
# Hint #1: sum() returns the sum of a vector.
# 	Ex: sum(c(1, 2, 3)) evaluates to 6
#
# Hint #2: length() returns the size of a vector.
# 	Ex: length(c(1, 2, 3)) evaluates to 3
#
# Hint #3: The mean of all the numbers in a vector is equal to the sum of all of
#		   the numbers in the vector divided by the size of the vector.
#
# Note for those of you feeling super clever: Please do not use the mean()
# function while writing this function!!!!!!

# ------------------------------------------------------------------------------------------------------

# Exercise 33
# You're going to write a function called "remainder." remainder() will take
# two arguments: "num" and "divisor" where "num" is divided by "divisor" and
# the remainder is returned. Imagine that you usually want to know the remainder
# when you divide by 2, so set the default value of "divisor" to 2. Please be
# sure that "num" is the first argument and "divisor" is the second argument.
#
# Hint #1: You can use the modulus operator %% to find the remainder.
#   Ex: 7 %% 4 evaluates to 3. 
#
# Remember to set appropriate default values!

remainder <- function(num, divisor) {
  # Write your code here!
  # Remember: the last expression evaluated will be returned! 
}

# ------------------------------------------------------------------------------------------------------

# Exercise 34
# You can pass functions as arguments to other functions just like you can pass
# data to functions. Let's say you define the following functions:
#
# add_two_numbers <- function(num1, num2){
#    num1 + num2
# }
#
# multiply_two_numbers <- function(num1, num2){
#	num1 * num2
# }
#
# some_function <- function(func){
#    func(2, 4)
# }
#
# As you can see we use the argument name "func" like a function inside of 
# "some_function()." By passing functions as arguments 
# some_function(add_two_numbers) willi evaluate to 6, while
# some_function(multiply_two_numbers) will evaluate to 8.
# 
# Finish the function definition below so that if a function is passed into the
# "func" argument and some data (like a vector) is passed into the dat argument
# the evaluate() function will return the result of dat being passed as an
# argument to func.
#
# Hints: This exercise is a little tricky so I'll provide a few example of how
# evaluate() should act:
#    1. evaluate(sum, c(2, 4, 6)) should evaluate to 12
#    2. evaluate(median, c(7, 40, 9)) should evaluate to 9
#    3. evaluate(floor, 11.1) should evaluate to 11

evaluate <- function(){
}
# ------------------------------------------------------------------------------------------------------

# Exercise 35
# The syntax for creating new binary operators in R is unlike anything else in
# R, but it allows you to define a new syntax for your function. I would only
# recommend making your own binary operator if you plan on using it often!
#
# User-defined binary operators have the following syntax:
#      %[whatever]% 
# where [whatever] represents any valid variable name.
# 
# Let's say I wanted to define a binary operator that multiplied two numbers and
# then added one to the product. An implementation of that operator is below:
#
# "%mult_add_one%" <- function(left, right){ # Notice the quotation marks!
#   left * right + 1
# }
#
# I could then use this binary operator like `4 %mult_add_one% 5` which would
# evaluate to 21.
#
# Write your own binary operator below from absolute scratch! Your binary
# operator must be called %p% so that the expression:
#
#       "Good" %p% "job!"
#
# will evaluate to: "Good job!"

"%p%" <- function(){ # Remember to add arguments!
}

#You made your own binary operator! Let's test it out. Paste together the strings: 'I', 'love', 'R!'
#| using your new binary operator.

#> "I" %p% "love" %p% "R!"
#[1] "I love R!"
# ------------------------------------------------------------------------------------------------------

# Exercise 36
# Define a function even_count that takes in input a vector v and return the number of even elements of it.
# ------------------------------------------------------------------------------------------------------

# Exercise 37
# Define a function m_create that takes 2 vectors v1, v2 and a variable n (by default = 10) that returns a matrix with n columns, 
# where each column is the concatenation of v1 and v2
# ------------------------------------------------------------------------------------------------------

# Exercise 38
# The ellipses can be used to pass on arguments to other functions that are
# used within the function you're writing. Usually a function that has the
# ellipses as an argument has the ellipses as the last argument. The usage of
# such a function would look like:
#
# ellipses_func(arg1, arg2 = TRUE, ...)
#
# In the above example arg1 has no default value, so a value must be provided
# for arg1. arg2 has a default value, and other arguments can come after arg2
# depending on how they're defined in the ellipses_func() documentation.
# Interestingly the usage for the paste function is as follows:
#
# paste (..., sep = " ", collapse = NULL)
#
# Notice that the ellipses is the first argument, and all other arguments after
# the ellipses have default values. This is a strict rule in R programming: all
# arguments after an ellipses must have default values. Take a look at the
# simon_says function below:
#
# simon_says <- function(...){
#   paste("Simon says:", ...)
# }
#
# The simon_says function works just like the paste function, except the
# begining of every string is prepended by the string "Simon says:"
#
# Telegrams used to be peppered with the words START and STOP in order to
# demarcate the beginning and end of sentences. Write a function below called 
# telegram that formats sentences for telegrams.
# For example the expression `telegram("Good", "morning")` should evaluate to:
# "START Good morning STOP"
# ------------------------------------------------------------------------------------------------------

# Exercise 39
# Let's explore how to "unpack" arguments from an ellipses when you use the
# ellipses as an argument in a function. Below I have an example function that
# is supposed to add two explicitly named arguments called alpha and beta.
# 
# add_alpha_and_beta <- function(...){
#   # First we must capture the ellipsis inside of a list
#   # and then assign the list to a variable. Let's name this
#   # variable `args`.
#
#   args <- list(...)
#
#   # We're now going to assume that there are two named arguments within args
#   # with the names `alpha` and `beta.` We can extract named arguments from
#   # the args list by using the name of the argument and double brackets. The
#   # `args` variable is just a regular list after all!
#   
#   alpha <- args[["alpha"]]
#   beta  <- args[["beta"]]
#
#   # Then we return the sum of alpha and beta.
#
#   alpha + beta 
# }
#
# Have you ever played Mad Libs before? The function below will construct a
# sentence from parts of speech that you provide as arguments. We'll write most
# of the function, but you'll need to unpack the appropriate arguments from the
# ellipses.

mad_libs <- function(...){
  # Do your argument unpacking here!

  # Don't modify any code below this comment.
  # Notice the variables you'll need to create in order for the code below to
  # be functional!
  paste("News from", place, "today where", adjective, "students took to the streets in protest of the new", noun, "being installed on campus.")
}
# ------------------------------------------------------------------------------------------------------

# Exercise 40
# Define a function sum_vector that takes in input a undetermined number of vectors and returns a
# the sum of all elements. NA values on vectors must be ignored
# Use ellipsis function
# ------------------------------------------------------------------------------------------------------

# Exercise 41
# Define a matrix_by_row function, with the same parameters of function "matrix()", 
# which fills the matrix by rows.
# Use ellipsis function
# ------------------------------------------------------------------------------------------------------

# Exercise 42
# Dataset mtcars
# Show the number of cars by the number of cylinders (cyl)
# Show th average consumption (mpg) by the number of cylinders (cyl)
# Use only one instruction
# ------------------------------------------------------------------------------------------------------

# Exercise 43
# Visualize the first rows of the dataset iris
# Define a function iris_mean that takes as parameter a type of iris specie and returns
# a vector of the means of each column of that specie
# Use the function sapply
# Apply the previous function to all species
# Use for loop
# ------------------------------------------------------------------------------------------------------

# Exercise 44
# Define a function "between" that takes as parameter: v, p and q. The function returns
# the number of elements of v that are between the values p and q (included).
# Define the following list: l <- list(c(23,232,112,2), c(23,33,33,22,32), c(12,34,65,67)).
# Apply the function between to l to obtain a vector of results
# Apply the function between to l to obtain a list of results
# ------------------------------------------------------------------------------------------------------

# Exercise 45
# Create a data frame df with 4 columns. The first with a sample of 20 numbers between 1 and 100,
# the second with the value TRUE if the corresponding element is even (FALSE otherwise), the third with 
# the cube of the x elements. The names of the columns are: x, even, square and cube
# Write the first 10 rows of df (with the function write.table) in a file "powers1.csv". Each element of a row 
# must be separated with a ",". Hold the column names but not the row names.
# Write the last 10 rows of df (with the function write.csv) in a file "powers2.csv".
# Don't hold the row names.
# Merge the contents of files "powers1.csv" and "powers2.csv" in a file "powers3.csv".
# Use the function "read.table" or "read.csv"
# ------------------------------------------------------------------------------------------------------

# Exercise 46
# Use "read.table" to read data from the link "http://data.princeton.edu/wws509/datasets/effort.dat" with following:
# - "a", "b", "c" as columns names
# - number from 1 to 10 as row names
# - rows from 3 to 12
# ------------------------------------------------------------------------------------------------------

# Exercise 47
# Read the first 10 lines of the html page "http://www.r-project.org" in the variable x
# Print the class of x
# Put the elements of x with odd position in y
# Write the elements of y in a file "y_output.txt"
# ------------------------------------------------------------------------------------------------------




