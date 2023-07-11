# connectEntrez
Obtain Gene Sequences straight from Entrez IDs.  

* Please set up an NCBI account at https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
* Enter your email prior to using.
* Input the specified genes dataset or small list of genes prior to starting the package.

### Python
```
pip install -r requirements.txt
```

### R
```
install.packages(c("seqinr", "stringr", "dplyr))
```


# Repo works in truncated steps: 
- NCBI has a limit on gene ids per request
- Find optimal truncated steps [10,50,100,etc.]
![image](https://github.com/Dwalczyk19/connectEntrez/assets/92831596/2131bc75-1c89-4b41-b30f-b13862e6b9b4)
