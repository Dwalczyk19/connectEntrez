# connectEntrez
Obtain Gene Sequences straight from Entrez IDs.  
* Works in two Steps:
  - Get data from Python specified by the user.
  - R parses it out and returns to you a single column file with gene names and their respective sequence
* Please set up an NCBI account at https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
* Enter your email prior to using.
* Input the specified genes dataset or small list of genes prior to starting the package.
* Use of virtual enviornment is recommended, but you don't have to.
  - This is due to problems with protobuf package 4.23.0 (if there are problems specify this package be installed in the requirments file)
* Download Anaconda here: https://www.anaconda.com/download

### Python
```
conda create -c conda-forge -n <ENV_NAME> spyder #using spyder as my default ide for venv

conda activate <ENV_NAME>

pip install -r requirements.txt
```

### R
```
install.packages(c("BiocManager", "seqinr", "stringr", "dplyr"))
BiocManager::install(("Biostrings"), force = TRUE)
```
* If installing from Ubuntu and having problems installing RCurl dependency, use command
```
sudo apt -y install libcurl4-openssl-dev
```


# Repo works in truncated steps: 
- NCBI has a limit on gene ids per request
- Find optimal truncated steps [10,50,100,etc.]
- Right now steps work in <= 1000. For e.g if parsing through 10000 genes you'll get 10 files.
![image](https://github.com/Dwalczyk19/connectEntrez/assets/92831596/2131bc75-1c89-4b41-b30f-b13862e6b9b4)
