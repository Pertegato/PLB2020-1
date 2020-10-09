#Exercise 48
#Create two numeric vectors x and y of length 3 and 5 respectively. 
#What happens when you execute their product?

#Exercise 49
#(a) Check if the current file is present in your current directory.
#(b) Then, choose another file from your current directory and get some informatin about it. 
#Do point (b) it in a single command.

#Exercise 50
#Write down the R code required to produce the following:
#(a) [1] 1 5 1 5 1 5 1 5 1 5 1 5 1 5 1 5 1 5 1 5 1 5 1 5 1 4 7 10 13 16 19
#(b) [1] 2 4 6 2 4 6 1 1 1 3 3 5 5 5 5

#Exercise 51
#Calculate the sum of i^2 for i=3:177

#Exercise 52
#Use paste to obtain the following
#[1] "Ab-Bc-Cd-De-Ef-Fg-Gh-Hi-Ij-Jk-Kl-Lm-Mn-No-Op-Pq-Qr-Rs-St-Tu-Uv-Vw-Wx-Xy-Yz-Za"
#What is the difference between the parameters sep and collapse?

#Exercise 53
#Sample 50 values from x<-rep(NA,300) and y<-1:300 and put the result in a vector z
#Find the proportion of NAs in z in a single command.

#Exercise 54
#Consider the vector z of the previous exercise. Modify it such that it contains only values greater than -50 
#but lower than 100(in a single command). NA values must be kept.
#What should be added in order to remove NA values?

#Exercise 55
#Create a vector x with named elements. Obtain the first element using three different ways. 
#Then modify x such that it does not contain the third element anymore. Do it by using names.

#Exercise 56
#Create a factor f containing the following levels: "BD Student" "MD Student". The repetitions should be 4 and 5 
#respectively.
#Then add a new level "PhD Student" for which the repetitions are still 0.

#Exercise 57
#Create a list l composed by the two following lists: 
l1<-list(name=c("Alberto","Antonio"),age=c(25,40),work=c("Student","Engineer"))
l2<-list(surname=c("Rossi","Bianchi"),nation=c("Italy","Malta"))
#(a)What happens if you use list to create l?
#(b)What happens when you type l$n and why? What should you type in order to obtain the nation component?

#Exercise 58
#Create a numeric matrix m 4x5.
#(a) Bind to m a numeric 4x2 matrix called n. Do both cbind and rbind work? If no, why?
#(b) Remove from the newly obtained m, the central column. Do this by exploiting the ncol function.
#(c) Multiply the new m matrix by n to obtain the row by columns product. Is it possible in this case?
#If not, how can we solve it?

#Exercise 59
#(a) Create a matrix m 5x5 containing letters by using the sample function such that an element can appear more than once.
#Check the class of one of the elements of m.
#(b) Add a column of numbers. Check the class of one of the new elements of m. Is there something strange about it?
#(c) Which data structure could you use instead of the matrix? Create an appropriate data structure n
#containing the same elements of m.
#Then name your columns C1,...,C5 and your rows R1,...,R5 (do not use c())
#(d) Access to an element of n in 3 different ways. Does it remind you of another data structure?

#Exercise 60
#Create a 3D array a 3x3x4 containing numbers from 1:18.
#(a) Does recycling occur?
#(b) Select only the first row and all but the first column. Leave the third dimension intact. Do not subscribe a.
#(c) Copy a to a2.What happens if you set dim(a2) to NULL?
#(d) Transform a into a matrix.

#Exercise 61
#Analyze the trees dataset.
#(a)Report its class,summary,structure,dimensions,size in bit,colnames and rownames.
#(b) What does R return if you type length(trees)?
#(c)Show the last 7 rows.
#(d) What is the frequency of the different values of Height? Round your result.

#Exercise 62
#Consider the following matrix 
m<-matrix(sample(-20:20,50,replace=T),nrow=5)
#(a) Check if it contains any 0.
#(b) Check if it contains all positive integers
#(c) Which elements are greater than 15? How many are them?

#Exercise 63
#Replace each negative element of m with the corresponding positive one. Replace positive elements
#with their cube.Use control structures.

#Exercise 64
#Build a cycle where each time you sample a new integer between -20 and 30.
#Whenever you encounter an even number (not 0), print the square of the sum of all even numbers sampled until
#that moment.
#When 0 is sampled, exit from the cycle.

#Exercise 65
#Define a function f that given a numeric vector, returns its maximum, minimum, median and mean value.
#Put these values in a vector, and return NULL if the input is not numeric. Test the function.

#Exercise 66
#Define a function euclidean_dist that given a matrix return the euclidean distance between its rows Do not use
#the dist function. Test the function(compare it with dist).

#Exercise 67
#Define a function fact that given a number returns its factorial. Try not to use for loops

#Exercise 68
#Define a function that given a character vector, returns the frequency of the characters.
#HINT: Use the function strsplit to split the string.
#Test it with the following: AAAACCGGTCGACGGGACGG

#Exercise 69
#Define a function with a variable number of input numeric arguments (minimum two) that does the following:
#(((a+b)*c)+d)*e and so on. Example: if the input is 4,5,6,7 what we want is ((4+5)*6)+7=61
#Test it using also length 1 vectors.

#Exercise 70
#Consider the dataset esoph, which contains data from a case-control study of esophageal cancer in France.
#Analyze the dataset (how many features, how many samples, the size of the dataset).
#Provide a summary of it.
#Then
#(a) See how many people are there per each age group
#(b) Find the classes of each feature
#(c) Compute how many controls are in each age group 
#(d) Find how many unique vals exist per each feature. Do it in a single command

#Exercise 71
#What is the difference between using vapply, sapply and lapply?
#When do sapply and lapply return the result in the same format?

#Exercise 72
#Consider the dataset state.x77, which contains statistics about the US. Consider also state.region,
#which indicates the location of each state.
#(a) Find the population for those states with an area greater than 100000 square miles(Single command).
#(b) Compute the mean of life Expectancy for each region(Northeast,South,North,Central,West)

#Exercise 73
#Make sure that that salary.txt is in your working directory.
#(a)Then execute the command below. What can you say about the first row of your file? How would you solve it?
sal<-read.table("salary.txt")
#(b) Select only those rows for which the degree is "doctorate" and put them into a new file called "doc_salary.txt".
#Remove the dg column. When writing the data to your file, specify new column names(sex,rank,year,year_since_deg,salary). 
#Use write.table.
#(c)Consider your doc_salary.txt file. Skip the first 5 lines and read only the first 10 among the remaining 
#ones. Use read.csv.

#Exercise 74
#Use the command sample to generate a vector x of 20 random numbers between 1 and 20. Do it using the fewest number
#of parameters.
#What can you say about the frequency of the numbers in x?
#Now create a vector y of 20 random numbers between 1 and 20 but with the parameter replace to TRUE.
#What is the difference between y and x?

#Exercise 75
#Generate using two different R functions, 100 random and independent events which result can be either success
#or failure, knowing that in 85% of cases you will fail. Encode,if needed, 0 as failure and 1 as success.
#Then verify the proportion of successes and compare it with the expected probability. Is it different?

#Exercise 76
#(a)Generate a vector x of 100 random numbers extracted from a normal distribution with 0 mean and 1 sd.
#Using the function hist() plot the vector. Give a name to your plot.
#(b) Using the same command, generate a vector y of 1000 numbers, but with mean 10 and 30 sd.
#Plot the histogram and compare it with the one you've obtained before.
#(c) Compute the mean, standard deviation and variance of y and compare them with the theoretical values
#(look at the parameters of your command). Are they completely equal?

#Exercise 77
#Find the number of hours that are missing from December 25,2017. Use opportune data formats. 
#Find the weekday of December 25,2017 as well.

#Exercise 78
#Consider the setosa objects of the iris dataset.
#Make a plot with the sepal.width on the x axis, and the sepal.length on the y axis. 
#Color your points dark green, and change their shape to crosses (hint:use the help)
#Assign a title to your plot and your axis.
