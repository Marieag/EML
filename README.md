### Data import instructions

This analytical pipeline imports data into R and performs exploratory statistics and a suite of differential analyses, such as DeSeq2 and constructing heat trees. A range of statistical tools is available, with the option of running ANOVAs or PERMANOVAs with post hoc tests.
The pipeline is based on R packages Phyloseq (1), vegan (2) and metacodeR (3). 

Input is required in the form of an S4 object, often colloquially named a phyloseq or physeq object. 
This is an R object consisting of an array of three separate matrices: 

**A table of counts** (in the case of microbial data, usually an OTU or ASV table)

![Example table of compound area counts, equivalent to a OTU table](https://github.com/Marieag/EML/blob/main/COMPOUNDS.jpg)

**A taxonomy table** describing the classification levels of each identified species or compound,

![Example table of chemical compound classifications, equivalent to a taxonomy table](https://github.com/Marieag/EML/blob/main/CLASS.jpg)

And **a metadata table**, detailing all extra information relevant to the study such as site, time points, temperature or pH values.

![Example table of sample metadata](https://github.com/Marieag/EML/blob/main/METADATA.jpg)

 
Fitting a chemical dataset to a tool developed for microbial ecology becomes a matter of processing the data to fit the requirements of the packages. 
When performing chemical statistical analysis, the data values used for analysis are extracted from the selected data points representing a single compound, given as the area under the peak in the HPLC chromatogram. These values are presented as decimal numbers; essentially values congruent with read counts in microbial diversity analysis. As such, a table of compound areas per sample can be used as the chemical equivalent of an OTU table. 
The taxonomy table equivalent presents a greater challenge – a classification system based on levels. In a taxonomy table, a common set of levels consist of classification of species into genera, families, orders, classes, phylae, and kingdoms, and this is of course not applicable to chemical compounds. Instead, compounds may be classified by their structure or their use, and by extension origin. This way, an equivalent of a taxonomy table – a compound classification table – can be constructed, and each identified compound be designated a source, as discussed in the accompanying paper (https://doi.org/10.1016/j.envres.2024.119242). This approach allows grouping of compound classes. If a classification table is constructed this way, and combined with the compound area table and a table of metadata, the result is an object that can be analysed using the numerous validated tools found in this optimised analysis pipeline, allowing for a comprehensive analysis of chemical data with the option of adding further phyloseq-compatible analyses or packages of interest.  



The data to be analyzed can be extracted from the same file (i.e. a curated excel file or similar containing all information, like “tutorial_dataset.csv”) or constructed from three separate files (“otu_table.csv”, “tax_table.csv” and “metadata.csv”).  For ease of syntax, the original terms from the Phyloseq packages are used in this instruction. As such, the area table is referred to as the “OTU table” (OTU’s – operational taxonomic units, in layman's terms single bacterial species – here being treated as data points equivalent to single compounds). 
The chemical classification table is referred to as the taxonomy (or tax) table. This allows for grouping of compounds and construction of cladograms – even if the classification is only two or three levels of chemical classes rather than taxonomic levels.

When constructing an S4 object (physeq object), keep in mind that sample id and compound/otu id MUST be identical across taxonomy table (chemical classification table), otu table (compound area table) and metadata.   

 
## Data import 
Load relevant .csv files, and check the files. This tutorial dataset should be a matrix of 38 compounds, with all acquired data in 76 rows. 
```
DATA <- read.csv2(file = "Tutorial_dataset.csv", header = TRUE, sep=",")
dim(DATA)
head(DATA)
```

If data rownames are numbers, they must be converted to character vectors to avoid issues later on. Consider streamlining both rownames and colnames this way, as is relevant for the data.

```
rowname <- as.character(as.matrix(DATA[1])); rowname[is.na(rowname)] <- "unidentified"
rowname <- paste("Compound", rowname, sep = "_")
```

__All rownames MUST be unique:__ 

```
rowname <- make.unique(rowname)
```

Consider backing up raw data to ease of access in case of errors. 

```
rawdata <- DATA
```

From here, the relevant data is extracted into separate objects: 



**Compound table (OTU table equivalent):**

Select relevant columns in the dataset, and save as numeric data. Set any NAs to 0. 

```
COMP <- as.data.frame(DATA[29:43])
COMP <- apply(COMP, 2, as.numeric); COMP[is.na(COMP)] = 0
Set rownames – compound and sample IDs MUST match across the S4 object. 
rownames(COMP) = rowname

head(COMP); dim(COMP) 
38 rows, 15 columns – 9 samples, 3 procedural blanks and 3 QC pools. 
```


**Classification table (Tax table equivalent):**

Same issues apply – select relevant columns and check row and column names. 
```
CLASS <- DATA[5:8]
rownames(CLASS) = rowname
head(CLASS)
```

**Metadata:** 
Metadata can be imported from a .csv file, or created directly in R. 

Metadata rownames MUST be identical to sample ID in the compound table. 

Creating tutorial metadata: 
```
temp <- as.data.frame(colnames(COMP))
site = c("A","A","A","A","D","D", "D","D","QC", "QC","QC","S","S","S","S")
type = c("Sample","Sample","Sample","Procedural_Blank","Sample", "Sample","Sample",
   "Procedural_Blank","QC","QC","QC", "Sample","Sample","Sample","Procedural_Blank")
PE = c(350000,350000,350000,350000,400000,400000,400000,400000,NA,NA,NA,12000,12000,12000,12000)
metadata <- cbind(temp,site,type,PE)
colnames(metadata) = c("Sample_ID","Site","Sample_Type","PE")
metadata
METADATA <- sample_data(metadata)
rownames(METADATA) = METADATA$Sample_ID
```

These three files can also be created in excel or similar, and imported directly as .csv files: 

```
COMP <- read.csv2(file = "compound_table.csv", header = TRUE, sep=",")
CLASS <- read.csv2(file = "classification_table.csv", header = TRUE, sep=",")
metadata <- read.csv2(file = " metadata.csv", header = TRUE, sep=",")
```

Following generation of the three data files, these can now be converted into the correct matrix formats and merged into a phyloseq (physeq) object, and any relevant elements from the analysis pipeline may be run. 

```
COMP_phy = otu_table(COMP, taxa_are_rows = TRUE)
CLASS_phy = tax_table(as.matrix(CLASS))
METADATA <- sample_data(metadata)
```

Create preliminary physeq object: 
```
physeq = phyloseq(COMP_phy, CLASS_phy)
```

Merge this object with the metadata:
```
physeq <- merge_phyloseq(physeq, METADATA)
```

Final physeq object, ready for the pipeline and other downstream analysis: 
```
physeq
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 38 taxa and 15 samples ]
tax_table()   Taxonomy Table:    [ 38 taxa by 4 taxonomic ranks ]
```

Check the Phyloseq and metacodeR tutorial pages for further info and help: 

[https://joey711.github.io/phyloseq/]

[https://github.com/grunwaldlab/metacoder]

1.   P. J. McMurdie and S. Holmes, “phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data,” PLoS One, vol. 8, no. 4, p. e61217, Apr. 2013, doi: 10.1371/journal.pone.0061217.

2.	 J. Oksanen et al., “vegan: Community Ecology Package. R package   version 2.0-10,” http://CRAN.R-project.org/package=vegan. 2013.

3.	 Z. S. L. Foster, T. Sharpton, and N. J. Grunwald, “MetacodeR: An R package for manipulation and heat tree visualization of community taxonomic data from metabarcoding,” bioRxiv, vol. 13, no. 2, p. 071019, 2017, doi: 10.1101/071019.


