Here is my week1 assignment
I have made a github repository and cloned it.

#1 samtools version


[axd6044@p-sc-2366 Week01]$ conda activate bioinfo
(bioinfo) [axd6044@p-sc-2366 Week01]$  samtools --version
samtools 1.22.1
Using htslib 1.22.1
Copyright (C) 2025 Genome Research Ltd.

Samtools compilation details:
    Features:       build=configure curses=yes 
    CC:             /opt/conda/conda-bld/samtools_1752528053426/_build_env/bin/x86_64-conda-linux-gnu-cc
    CPPFLAGS:       -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /storage/work/axd6044/.conda/envs/bioinfo/include
    CFLAGS:         -Wall -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /storage/work/axd6044/.conda/envs/bioinfo/include -fdebug-prefix-map=/opt/conda/conda-bld/samtools_1752528053426/work=/usr/local/src/conda/samtools-1.22.1 -fdebug-prefix-map=/storage/work/axd6044/.conda/envs/bioinfo=/usr/local/src/conda-prefix
    LDFLAGS:        -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/storage/work/axd6044/.conda/envs/bioinfo/lib -Wl,-rpath-link,/storage/work/axd6044/.conda/envs/bioinfo/lib -L/storage/work/axd6044/.conda/envs/bioinfo/lib
    HTSDIR:         
    LIBS:           
    CURSES_LIB:     -ltinfow -lncursesw

HTSlib compilation details:
    Features:       build=configure libcurl=yes S3=yes GCS=yes libdeflate=yes lzma=yes bzip2=yes plugins=yes plugin-path=/storage/work/axd6044/.conda/envs/bioinfo/libexec/htslib htscodecs=1.6.4
    CC:             /opt/conda/conda-bld/htslib_1752522550715/_build_env/bin/x86_64-conda-linux-gnu-cc
    CPPFLAGS:       -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /storage/work/axd6044/.conda/envs/bioinfo/include
    CFLAGS:         -Wall -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /storage/work/axd6044/.conda/envs/bioinfo/include -fdebug-prefix-map=/opt/conda/conda-bld/htslib_1752522550715/work=/usr/local/src/conda/htslib-1.22.1 -fdebug-prefix-map=/storage/work/axd6044/.conda/envs/bioinfo=/usr/local/src/conda-prefix -fvisibility=hidden
    LDFLAGS:        -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/storage/work/axd6044/.conda/envs/bioinfo/lib -Wl,-rpath-link,/storage/work/axd6044/.conda/envs/bioinfo/lib -L/storage/work/axd6044/.conda/envs/bioinfo/lib -fvisibility=hidden -rdynamic

HTSlib URL scheme handlers present:
    built-in:    file, preload, data
    S3 Multipart Upload:         s3w+https, s3w+http, s3w
    Amazon S3:   s3+https, s3, s3+http
    libcurl:     gophers, smtp, wss, smb, rtsp, tftp, pop3, smbs, imaps, pop3s, ws, ftps, https, ftp, gopher, sftp, imap, http, smtps, scp, dict, mqtt, telnet
    Google Cloud Storage:        gs+http, gs+https, gs
    crypt4gh-needed:     crypt4gh
    mem:         mem
(bioinfo) [axd6044@p-sc-2366 Week01]$ 

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
