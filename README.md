# Genome Contact Map Explorer - gcMapExplorer

It is platform to visualize and analyze the contact maps that are generated from Hi-C experiments. This package is developed by considering the huge size of contact maps at very fine resolution. It contains
  * Graphical User Interface - Several windows like applications to perform tasks.
  * Command Line Interface - Several commands to perform tasks.
  * Application Programming Interface - It can be used to perform analysis by anyy mathematical operations through programming.


## Features:
  * Support for **huge contact maps** - Use of Disk instead of RAM
    * Matrices/arrays are stored in Disks - mathematical operations by directly reading/writing from/to Disks, **without loading them into RAM**
  * A browser with rich interfaces for **Comparative** and **Interactive** visualization of **two dimensional contact maps** with **one dimensional genomic datasets** such as DNase-seq, ChIP-seq, RNA-seq etc.
  * Contact maps can be **zoomed in/out** from finest resolution to whole chromosome level. 
  * Rich customizations of **color scale for contact maps** visualization
  * Rich customizations of **X- and Y- axis properties**.
  * Publication ready images at one click.
  * Normalization of contact maps by 
    * **Iterative Correction** (IC)
    * **Knight-Ruiz Matrix Balancing** (KR)
    * **Distance-Frequency**
  * A **new file format** for contact map  and genomic datasets:
    * **Protable**, **platform independent** and can be read through C/C++, JAVA, Python and R programming language.
    * **Very fast to read** - fast browsing of contact maps and genomic datasets
  * Another file format for chormosomal contact map - much faster than above format to read/write but not compact
  

# Coming soon !!!
  
