Here is my week1 assignment
I have made a github repository and cloned it.

# Week 1 Assignment

## 1. Samtools Version
### INPUT
```bash
conda activate bioinfo
samtools --version
```
### OUTPUT
```
samtools 1.22.1
Using htslib 1.22.1
Copyright (C) 2025 Genome Research Ltd.
```
## 2. Making nested directories
### INPUT
```bash
mkdir newdir
cd newdir/
mkdir dir1
cd dir1/
mkdir dir2
ls
```

## 3. Making files in different directories
```
ls
newdir  README.md
```
### Input
```
cd newdir/
touch newdir.txt
ls
```
### Output 
```
dir1  newdir.txt
```
## 4. Relative and absolute paths
### relative path 
### Input
```
cd Week01/newdir/dir1/
ls ../../..
```
### Input
```
README.md  Week01
```
## absolute path
### Input
```
cd ~/scratch/Bioinformatics-course-projects/Week01/newdir/dir1/dir2/'
```
### Output
```
$ 
```
