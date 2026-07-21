### pathPhynder Installation

#### Conda installation

Create an environment (must have bioconda channel set up, see https://bioconda.github.io/index.html).
Note this comes with phynder and all other dependencies pre-installed.

``` bash
# set up bioconda (note: channel_priority set to flexible)
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible

# create environment
conda create -n pathPhynder -c bioconda pathphynder samtools=1.20
```

Once installed, activate the environment.

```bash
conda activate pathPhynder
pathPhynder -h
```

Once finished, deactivate the environment.

```bash
conda deactivate
```

#### Manual installation

Prerequisites:

 - R language (https://www.r-project.org/)
 - python 3 (https://www.python.org/)
 - samtools (http://www.htslib.org/)


Install the following R packages:
```
install.packages('optparse')
install.packages('phytools')
install.packages('scales')
```

1) Download pathPhynder to your computer
```
git clone https://github.com/ruidlpm/pathPhynder.git
```

2) Add the following line to your ~/.bash_profile (create one if necessary). Replace <path_to_pathPhynder_folder> with the location of the pathPhynder folder in your system.
```
alias pathPhynder="Rscript <path_to_pathPhynder_folder/pathPhynder.R>"
```
For example, if you have downloaded the folder to your ~/software/ directory, then you would add the following lines to ~/.bash_profile.
```
alias pathPhynder="Rscript ~/software/pathPhynder/pathPhynder.R"
```
and then:
```
source ~/.bash_profile
```

3) Test the installation.

```
pathPhynder -h
```

### phynder Installation

For more details, you can visit the original repository:
https://github.com/richarddurbin/phynder

```
git clone https://github.com/samtools/htslib.git
cd htslib
make lib-static htslib_static.mk
cd ..
git clone https://github.com/richarddurbin/phynder.git
cd phynder
make
make install
```

(Note: If on macOS, you may need conda install clangxx_osx-64)
