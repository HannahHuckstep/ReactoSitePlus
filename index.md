# PhlashyNAMe Tutorial

## Introduction 

PhlashyNAMe is a command line tool for downstream analysis Proteomics and Phosphoproteomics data. The tool can be run on linux, Mac and Windows operating systems. What makes PhlashyNAMe unique is it's ability to use phospho-peptide level information to analyse your data. Reactome annotates proteins in specific phosphorylation states as separate entities, therefore these phosphoproteins can take part in different signalling cascades depending on their phosphorylation state. Proteins in these differing modification states are hereby referred to as proteoforms. Taking advantage of this structure, PhlashyNAMe takes an input file of phosphopeptides and maps the phosphorylations onto proteins and complexes according to a set of guidelines explained in detail in the figure below. It works by assigning a confidence score and an abundance score to each mapped protein and complex. 
#### Confidence Score 
This confidence score summarises support for a given phosphorylation state of a protein over its other phosphorylation states based on observations from the data. 

![Support Score](https://user-images.githubusercontent.com/9949832/120878866-98292680-c602-11eb-9e33-aaf8e3549ee1.png)

In Reactome, each modified protein is treated as a separate entity to the unmodified protein. In this document these entities are referred to as proteoforms, as they are different versions of the same protein and participate in different parts of the network accordingly. Proteoforms are highlighted in green throughout this figure. Across all proteoforms (green), all recorded phosphorylations are aggregated into one group called Points of Interest (orange) as shown in above. The phosphorylations recorded in the network in this example (PA,PB, and PC) are shown in purple while, phosphorylations found in the data (and not in the network) are shown in red. 
Looking at the rules developed to score each proteoform, we first score each peptide and take the average of the peptide scores to get the proteoform Confidence score. In rule 1 of the peptide score, a match means at the amino acid in the proteoform, and the peptide have the same status (both phosphorylated or both unphosphorylated). We give a score of 0.5 for a mismatch because that peptide supports the existence of that proteoform, but we know that peptide contains a P.O.I. and will support another proteoform better. Next, looking at Rule 2, we take the average of all possible matches. In Rule 3, we multiply by 0.9 to reflect that the database takes precedent over the data. However, we also reduce the weight of the score to reflect the extra unknown phosphorylation. Finally, in Rule 4, because a peptide mapping perfectly to a modified proteoform is uncommon and of interest we multiply the score by 1.5 to highlight the match. However, we do not highlight this if the proteoform is unmodified as perfect matches would be unmodified peptides (which mostly occur from unspecific binding). 


#### Abundance Score 
Concurrently but separately, abundances (e.g., intensity, SILAC ratio, p-value, log fold change) are mapped onto all proteins and complexes across the network. The process of computing this score is illustrated in the figure below. The default abundance score mapping takes the abundance from the peptide with the highest confidence score (if there are multiple peptides with the same confidence score the average of the abundances are taken). The user can also choose to compute the mean or median abundance for each phosphopeptide mapped to a protein. 

![Abundance score](https://user-images.githubusercontent.com/9949832/120879240-5b126380-c605-11eb-8d3f-6b691f454cfd.png)

When the confidence score is mapped, the abundance score is mapped simultaneously, although the scores remain separate. The user has 4 options for mapping the abundance score as described above. 

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
conda create --name PhlashyNAMe \
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
* the name of the **i**nput **O**WL **f**ile \[-iof\] of [Reactome](https://reactome.org/download/current/biopax.zip) that you would like to build. 
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
    * Each option is explained above in the abundance score section. The default option is 'HighestSupport'.  



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


## Analysing Mapped Network 

Now that the network has the data mapped to it, there are a number of ways to analyse the mapped network. 
* Traversal analysis 
* Neighbourhood analysis 
* Find the shortest path between two proteins 
* Minimal Connection Network
* Manual investigation via cytosscape or Neo4j

### Traversal Analysis 


### Neighbourhood Analysis 

### Shortest Path 

### Minimal Connection Network 

### Visualization 
#### Writing the mapped database to a SIF 
After the database is mapped to you can write it to a Simple Interaction Format (SIF) file to import into [cytoscape](https://cytoscape.org/) along with an attribute file. The SIF file will be names SIF.sif, which can be renamed but must keep the .sif extension. An example of a SIF file looks like: 

```
21	INPUT	20
23	PHOSPHORYLATION	21
20	OUTPUT	24
```
Which we would read as the node with the id ```21``` is an ```INPUT``` to the reaction node with the id ```20```, or the node with the id ```23``` is a ```PHOSPHORYLATION``` on the node ```21```.

The attribute file will contain all attributes associated with each node. An example of an attribute file looks like: 

Node_ID | Database_ID | Display_Name | Type | Database_Link | Location | Status | Kinase | Transcription_Factor | Cell_Surface_Receptor | UniProt_Gene_Name | Integrated | ABUNDANCE_SCORE_wt | SUPPORT_SCORE_wt | ABUNDANCE_SCORE_stim | SUPPORT_SCORE_stim
--------|-------------|--------------|------|---------------|----------|--------|--------|----------------------|-----------------------|-------------------|------------|--------------------|------------------|----------------------|-------------------
21 | Protein2090 | p-S568-MLXIPL | Protein | http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=R-HSA-163687.1 | nucleoplasm | | | TRANSCRIPTION_FACTOR | | MLXIPL | | 3.2 | 0.9 | -1.5 | 0.5
24 | SmallMolecule848 | Pi | SmallMolecule | http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=R-ALL-113550.4 | nucleoplasm | | | | | | | | | | 
20 | BiochemicalReaction1310 | [Dephosphorylation of pChREBP (Ser 568) by PP2A] | BiochemicalReaction | http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=R-HSA-164056.2 | | | | | | | | | | | 
23 | Protein2090_p-S_568 | p-S_568 | p_S | | 568 | | | | | | | | | | | 


Where the attributes are associated with the node ids found in the SIF file. The first line contains all of the attributes that are found in that database.

Attribute name | Description
---------------|-----------
Node_ID | The unique node number used to link nodes to attributes in cytoscape
Database_ID | The unique node ID given by Reactome 
Display_Name | The name of the node 
Type | The type of node 
Database_Link | The link to the database that has more information on the node. Can link back to Reactome, PhosphoSitePlus or UniProt
Location | The cellular location of the node 
Status | The status of UniProt nodes, can be nothing (if the database was not updated), 'Current', 'Updated', or 'Deleted?' - the question mark indicates further review is needed. The UniProt Id could be old and outdated or not consistent with the species in the database.  
Kinase | If the node is a Kinase this column with have the values 'KINASE' otherwise it will be left blank. These are the lists of kinases for [Human](https://www.uniprot.org/uniprot/?query=keyword:%22Kinase%20[KW-0418]%22&format=list&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22), and  [Mouse](https://www.uniprot.org/uniprot/?query=keyword:%22Kinase%20[KW-0418]%22&format=list&fil=organism:%22Mus%20musculus%20(Mouse)%20[10090]%22).
Transcription_Factor | If the node is a Transcription Factor this column with have the values 'TRANSCRIPTION_FACTOR' otherwise it will be left blank. These are the lists of Transcription factors for [Human](https://www.uniprot.org/uniprot/?query=goa:(%22DNA-binding%20transcription%20factor%20activity%20[3700]%22)%20(reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22)&format=list), and [Mouse](https://www.uniprot.org/uniprot/?query=goa:(%22DNA-binding%20transcription%20factor%20activity%20[3700]%22)%20(reviewed:yes%20organism:%22Mus%20musculus%20(Mouse)%20[10090]%22)&format=list).
Cell_Surface_Receptor | If the node is a Cell Surface Receptor this column with have the values 'CELL_SURFACE_RECEPTOR' otherwise it will be left blank. These are the lists of Cell surface receptors for [Human](https://www.uniprot.org/uniprot/?query=goa:(%22cell%20surface%20receptor%20signaling%20pathway%20involved%20in%20cell-cell%20signaling%20[1905114]%22)%20(reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22)&format=list) and [Mouse](https://www.uniprot.org/uniprot/?query=goa:(%22cell%20surface%20receptor%20signaling%20pathway%20involved%20in%20cell-cell%20signaling%20[1905114]%22)%20(reviewed:yes%20organism:%22Mus%20musculus%20(Mouse)%20[10090]%22)&format=list).
UniProt_Gene_Name | The gene name associated with that UniProt ID from UniProt.
Integrated | If the node is from PhosphoSitePlus the value will be True. 
ABUNDANCE_SCORE_ | The abundance score mapped from the data for a node. There may be multiple Abundance scores in a single database derived from multiple experiments. Each row labelled ABUNDANCE_SCORE_ will contain the experiment name in the suffix after the last '_'. 
SUPPORT_SCORE_ | The support score mapped from the data for a node. There may be multiple Support scores in a single database derived from multiple experiments. Each row labelled SUPPORT_SCORE_ will contain the experiment name in the suffix after the last '_'. 

To write the SIF and attribute files the java command requires the input database directory and the path to the direcory you'd like the output created in. The command is as follows: 
``` 
java -jar ./path/to/jars/RactoSitePlus.jar -m WriteDBtoSIF -idb ./path/to/graph/ -op ./path/to/output/
```

To then load the SIF file into cytoscape, first open the cytoscape application. Next click the network button (highlighted in red in the figure below).

![Cytoscape Network Button](https://user-images.githubusercontent.com/9949832/120894578-81172280-c65c-11eb-9e3a-d2da0ccb2ff5.png)

Next, navigate to your SIF file and open it. I do not recommend a network view is made at this point as the network is incredibly large and cytoscape often crashes while attempting to make a network view this large. Following this, click the attribute file button (highlighted in red in the figure below). 

![Cytoscape Attribute Button](https://user-images.githubusercontent.com/9949832/120894665-fb47a700-c65c-11eb-8576-4179716395da.png)

Navigate to your attibute file and open it. As in the figure below, make sure the Node_ID column is chosen as the key column. You can also change the value type of the attributes. One change you may like to make is to ensure the ABUNDANCE_SCORE_ and SUPPPORT_SCORE_ columns are set to numeric values. 

![Cytoscape attributes](https://user-images.githubusercontent.com/9949832/120894779-950f5400-c65d-11eb-871a-9d3e8a0c3597.png)

Now your database should be ready to be explored. You may select nodes from the node table (make sure to left click the highlighted nodes and choose 'Select nodes from selected rows') and create a newtork view using the button highlighted in the figure below. 

![Cytoscape new network button](https://user-images.githubusercontent.com/9949832/120895262-8164ed00-c65f-11eb-9c49-383c636fe57f.png)






