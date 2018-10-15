## Homework 2  
Please design questions suitable for an upper level undergraduate course to illustrate the concepts we covered in class.
Then, provide the answer key for your own questions (i.e. answer your own questions). Below, a prompt is given for the 
concept you are to design a question from. When given an option, try your best to choose the one that you're least familiar 
with so you can hone your own understanding.
1. Ask a question that requires a student to understand navigation and manipulation of directories in a filesystem. Your 
question should require an answer using at least the following commands/concepts: cd, ../, mkdir, rmdir  
>**Question:**  Write a code on the terminal that would allow you to make a new directory named *"MyDir"*. After creating an
 empty document called *"First_Document"* which needs to be copied to the parent directory, this newly created directory 
 needs to be deleted.  
**Answer:** 
Option 1 :
<pre><code>$ mkdir MyDir
$ cd ./MyDir
$ touch First_Document
$ cp First_Document ../
$ rm First_Document
$ cd..
$ rmdir MyDir
</code></pre>  
Option 2:
<pre><code>$ mkdir MyDir
$ cd ./MyDir
$ touch First_Document
$ cp First_Document ../
$ cd..
$ rm -r MyDir
</code></pre> 
2. Ask a question that requires a student to understand the difference between accessing a column in a matrix with text 
indices versus accessing a column in a data frame with text indices. Your question should require an answer comparing the 
following: mymatrix[,'col1'] vs. mydf[,'col1'] vs. mydf['col1'] vs. mydf$col1 vs. mydf[['col1']].  
>**Question:**  After creating a dataframe and a matrix of your choice, explain various ways of calling the first column for 
each of those, without using the index value of 1.  
**Answer:** The code given below explains various ways of calling the first column of both the data frame and matrix without 
actually using the index value for the said column.  
There is essentially only a single way for calling a column or a row of a matrix as given below. The command "_mymatrix[,'col1',drop=F]_" gives the first column of the matrix in a matrix form and doesn't drop the dimentions of the matrix 
to that of a vector.  
For a dataframe however, there is more than a single way to do the same as the matrix can be subsetted like a list as well as a 
matrix. The commands "_mydf[,'Name']_", "_mydf[['Name']]_" and "_mydf$Name_" return the same information (in the same format i.e as Factors)  which is the column named "_Name_" as a vector. However, the command "_mydf['Name']_" returns the column as it is, and the information can be considered as another dataframe (sub-dataframe).
<pre><code>
require('RCurl')
mlbweightdfurl <- 'http://goo.gl/rih9v9'
mlb.weight.df <- textConnection(getURL(mlbweightdfurl, followlocation  = TRUE))
mydf<- read.table(mlb.weight.df, header = TRUE)
mydf[,'Name'] 
mydf['Name'] 
mydf$Name 
mydf[['Name']]

mymatrix <- matrix(data = 1:12, nrow = 3, ncol = 4, byrow = TRUE)
dimnames(mymatrix) <- list(c('row1', 'row2', 'row3'), c('col1', 'col2','col3','col4'))
print(mymatrix)
mymatrix[,'col1',drop=F] 
</code></pre> 
3. Ask a question that requires a student to understand how to share access to a directory and a file in that directory 
on a Unix/Linux filesystem from their home directory with a colleague without exposing the user's entire directory. Your 
question should require an answer using chmod {u,g,o}{+,-}{r,w,x} (not using octal permissions).
   * **Hint 1:** in order to view a file, all of its parent directories must be executable
   * **Hint 2:** in order to view a file, the file itself must be readable  
>**Question:**  Write a code to change the permissions for a new directory called "New_Dir" (and a file inside the said directory 
called "New_File") in your parent directory, which would now be easily acessible to other people without using the octal permissions.  
**Answer:**  
<pre><code>$ cd ~
$ mkdir New_Dir
$ ls -l
$ chmod go+x New_Dir
$ cd ./New_Dir
$ touch New_file
$ ls -l
$ chmod go+r New_file
$ ls -l
</code></pre>
