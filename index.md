# PhlashyNAMe Tutorial

## Introduction 

PhlashyNAMe is a command line tool for downstream analysis Proteomics and Phosphoproteomics data. The tool can be run on linux, Mac and Windows operating systems. 

## Dependencies 
* Java 7
* R
  * ggplot 
  * ... 

## Setting up 

1. Clone this repo (as below) or create a new directory and place the provided scripts 

```
git clone https://github.com/HannahHuckstep/Db_Compare.git
```
2. If you do not have the required depndencies the following `conda` command will create an environment called `PhlashyNAMe` with all the dependencies installed. 

```
conda create --name DbCompareConda \
  --channel conda-forge \
  --channel bioconda \
  r=3.6 \
  r-upsetr \
  r-sna \
  r-plotrix \
  r-ggplot2 \
  bioconductor-clusterprofiler \
  bioconductor-org.hs.eg.db \
  bioconductor-org.hs.eg.db \
  pandoc \
  openjdk=8
```

3. [OPTIONAL] download the most current version of Reactome and PhosphositePlus. There are currently pre-made integrated and non-integrated databases which can be found in the databases file in this repo for both human and mouse. 

4. Map your data as shown below 

5. Data can then be analysed with useing the number of functions provided below with detailed examples. 

## Tool options and commands 
To start, navigate into the repo directory and type the following command to view all of the tool options: 

```
java -jar jars/ReactoSitePlus.jar -h
```
The first 2 arguments are named -h (or --help) to access the above again (This will be useful later). While the second named option -m (or --mode) is used to specify which function you would like to perform. the current options are listed in the '{}' and then again below in a list form. As an example, the first option is: 
* "CreateDB", takes an  OWL  file  [-iof],  an  output  path  [-op],  an  optional  update  boolean  [-u]  (can  be  T  or  F,  default  is  T),  and the
                         species of graph you'd like to make [-s] (can be human (h) or mouse(m))

So if you would like to create a database you would put the `CreateDB` option after the -m when performing the command. However, this function requires a few parametes to be specified which are shown in the '\[\]' 

For this option you will need to specify:
* the name of the **I**nput **O**wl **F**ile \[-iof\] of [Reactome](https://reactome.org/download/current/biopax.zip) that you would like to build. 
* the **o**ut***p***ut directory where you'd like your graph to reside \[-op\]
* optionally you can specify if you'd like your database to be **u**pdated \[-u\] which can only be set to 'T' (true) or 'F' (false). (updating the database will update secondary UniProt accessions (reffered to here as UniProt ID's) to current UniProt ID's, it will also annotate any deleted or non-human UniProt ID's as deleted) 
* you must also specify which **s**pecies of database you'd like to create \[-s\], which can be set to Human with 'h', 'human', '9606', or Mouse 'm', 'mouse', '10090'

Thus, the final command would look something like this: 
```
java -jar ./path/to/jars/RactoSitePlus.jar -m CreateDB -iof ./path/to/file/Reactome.owl -op ./path/to/graph/ -u T -s h
```

**All comands will start with** ```java -jar ./path/to/jars/RactoSitePlus.jar```

Finally, below all the commands are all of the parameters required to run each command, and an explanation of that they are used for. e.g.,
* --input_owl_file [INPUT_OWL_FILE], -iof [INPUT_OWL_FILE]
                         The OWL file to input

## Mapping Phosphoproteomic Data 

Now we should be ready to map our data! Data can be mapped to pre-built neo4j databases located in the databases folder, or you can build your own up-to-date database as shown above (and again with more detail below). 
This command is shown in the help guide as: 
* "MapPeptides",  takes  an  input  database  [-idb]  and  an  output  path  [-op],   a   file  to  map  onto  the  database  [-idf],  and  the  optional
                         Abundance Score mapping method preferred [-as] ("HighestSupport" is defalut)

Which needs the following parameters specified: 
* the input database directory \[-idb\]
* the output directory for the mapping report \[-op\]
* the file of data you would like to map onto the database \[-idf\]
* The mapping option you'd prefer. There are 5 options 'HighestSupport', 'Max', 'Mean', 'Median', and 'Extreme'



To map the data you can run the following command
``` 
java -jar ./path/to/jars/RactoSitePlus.jar -m MapPeptides -idb ./path/to/graph/ -op ./path/to/output/ -idf ./path/to/data.txt -as HighestSupport
```
                         
In order to map your data you will need an input file with the following 4 columns; 
1. A column where each cell contains a single UniProt ID corresponding to a peptide. *The name of the column cannot contain spaces*. 
    1. e.g., leading_razor_protein 
2. A column specifying the modified peptide sequence. *The name of the column cannot contain spaces*. 
    1. The modified peptide may be one of 2 different formats. Either \_(ac)AAAITDM(ox)ADLEELSRLS(ph)PLPPGS(ph)PGSAAR\_  or AAAITDMADLEELSRLpSPLPPGpSPGSAAR
3. A column containing the value you would like mapped. *The name of the column cannot contain spaces*. 
    1. e.g., the column that holds: p-value, SILAC_ratio, log2_intensity
4. A column you would like the mapped values associated with, such as the experiment name or time. *The name of the column cannot contain spaces*. 
    1. e.g., the column that holds:  stimulated, unstimulated, time_point_1




