Here is my week1 assignment
I have made a github repository and cloned it.

# Week 1 Assignment

## 1. Samtools Version
```bash
conda activate bioinfo
samtools --version

samtools 1.22.1
Using htslib 1.22.1
Copyright (C) 2025 Genome Research Ltd.



##2. Making nested directories 
(bioinfo) [axd6044@p-sc-2366 Week01]$ mkdir newdir
(bioinfo) [axd6044@p-sc-2366 Week01]$ cd newdir/
(bioinfo) [axd6044@p-sc-2366 newdir]$ mkdir dir1
(bioinfo) [axd6044@p-sc-2366 newdir]$ cd dir1/
(bioinfo) [axd6044@p-sc-2366 dir1]$ mkdir dir2
(bioinfo) [axd6044@p-sc-2366 dir1]$ ls
dir2

##3. Making files in different directories

(bioinfo) [axd6044@p-sc-2366 Week01]$ ls
newdir  README.md
(bioinfo) [axd6044@p-sc-2366 Week01]$ cd newdir/
(bioinfo) [axd6044@p-sc-2366 newdir]$ touch newdir.txt
(bioinfo) [axd6044@p-sc-2366 newdir]$ ls
dir1  newdir.txt

##4. Relative and absolute paths

##relative path 
(bioinfo) [axd6044@p-sc-2366 Bioinformatics-course-projects]$ cd Week01/newdir/dir1/
(bioinfo) [axd6044@p-sc-2366 dir1]$ ls ../../..
README.md  Week01

##absolute path

[axd6044@p-sc-2366 Bioinformatics-course-projects]$ cd ~/scratch/Bioinformatics-course-projects/Week01/newdir/dir1/dir2/'
[axd6044@p-sc-2366 dir2]$ 
