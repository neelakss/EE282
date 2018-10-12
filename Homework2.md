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
$ rmdir My_Dir
</code></pre>  
Option 2:
<pre><code>$ mkdir MyDir
$ cd ./MyDir
$ touch First_Document
$ cp First_Document ../
$ cd..
$ rm -r My_Dir
</code></pre> 
2. Ask a question that requires a student to understand the difference between accessing a column in a matrix with text 
indices versus accessing a column in a data frame with text indices. Your question should require an answer comparing the 
following: mymatrix[,'col1'] vs. mydf[,'col1'] vs. mydf['col1'] vs. mydf$col1 vs. mydf[['col1']].  
>**Question:**  
**Answer:**
3. Ask a question that requires a student to understand how to share access to a directory and a file in that directory 
on a Unix/Linux filesystem from their home directory with a colleague without exposing the user's entire directory. Your 
question should require an answer using chmod {u,g,o}{+,-}{r,w,x} (not using octal permissions).
   * **Hint 1:** in order to view a file, all of its parent directories must be executable
   * **Hint 2:** in order to view a file, the file itself must be readable  
>**Question:**  Write a code to change the permissions for a new directory called "New_Dir" in your parent directory, which
would now be easily acessible to other people without using the octal permissions.  
**Answer:**
