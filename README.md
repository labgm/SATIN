# SATIN - mini and  microSATellite IdentificatioN Tool
- Satin is a command line tool for Mini and Microsatellite identification with many new features to optimize SSR prospection. It can be executed in Linux Terminal. It is very intuitive to use and only requires a fasta or a gbk to run.

# How to Install
## 1. Download SATIN’s folder. 
The executable main.py is the file the user should run in order to make SSR analysis, but it requires the supporting files of the folder to run the codes.

1.1. The download of the folder can be done by clicking in “Code > Download Zip” at SATIN’s github page.

![image](https://github.com/labgm/SATIN_c/assets/101668229/faae9662-6651-4e5c-885a-e684230493ba)

 
1.2. The download of the folder can also be downloaded through Linux Terminal by using git clone:
```sh
git clone https://github.com/labgm/SATIN_c.git
```


## 2. Verify if the dependencies are installed are your computer. 
Satin depends on those dependencies to run.
   
  2.1 Dependencies:

  Recommended installation of dependecies (conda or mamba)

```sh
mamba env create -f SATIN.yml
```

  To activate the environment

```sh
mamba activate SATIN
```
  Alternative installation
   
* pip3
```sh
sudo apt install python3-pip
```

* Pandas 
```sh
pip3 install pandas
```
* Numpy 
```sh
pip3 install numpy
```
* Biopython
```sh
pip3 install biopython
```
* Collections 
```sh
pip3 install collection
```
* Zlib 
```sh
sudo apt-get install libz-dev
```
* Bzlib.h
```sh
sudo apt-get install libbz2-dev
```

2.2 Optional R dependencies
SATIN uses R to calculate abundance and statistics. It’s adviseable to install R in your computer, and install the R packages dplyr and tidyr.
*Install R
```sh
sudo apt-get install r-base
```

* Install dplyr for R
```sh
install.packages('dplyr')
```

* Install tidyr for R 
```sh
install.packages('tidyr')
```

3.3 In some cases when using the Linux OS, you may need to install GLIBC. When required, an error message may occur when trying to run the application.
> Error message: "**version GLIBC_{VERSION} not found..**"

* In this case, you need to install GLIBC.
```sh
wget http://ftp.gnu.org/gnu/libc/glibc-{VERSION}.tar.gz
tar -xvf glibc-{VERSION}.tar.gz
cd glibc-{VERSION}
mkdir build 
mkdir glibc-{VERSION}-install
cd build
~/glibc/glibc-{VERSION}/configure --prefix=$HOME/glibc/glibc-{VERSION}-install
make -j
make install
```

# To Execute SATIN 
Satin is very intuitive to use. The user needs to execute the main executable file (main.py) and write in the command line the specifications of input (“-f”, “-mf”, “-fs”, “-g” or “-gs”)  and output (-o), what files are to be analyzed and Where the output should be placed after the analysis is done. 
Satin can use as input Genbank or Fasta Files. If Genbank files are used as input, SATIN generates additional outputs.
Multiple files can be used at input at the same run, as long as they are of the same type.
 
## 1.	How to run using Genbank files as Input
* Single GBK file
```sh
python3 main.py -g input_folder/filename.gbk -o output/
```
* Multiple GBK files
```sh
python3 main.py -gs input_folder/ -o output/
```
## 2.	How to run using Fasta files as Input
* Single FASTA file
```sh
python3 main.py -f input_folder/filename.fna -o output/
```
* MultiFASTA file
```sh
python3 main.py -mf input_folder/filename.fna -o output/
```
* Multiple FASTA/MultiFASTA files
```sh
python3 main.py -fs input_folder/ -o output/
```