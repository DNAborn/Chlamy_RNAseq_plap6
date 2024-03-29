---
  title: "RNAseq plab6"
author: "Kelterborn"
date: "2024-03-15"
output:
  md_document:
  variant: gfm
toc: true
always_allow_html: true
editor_options: 
  chunk_output_type: console
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = T,
                      error=FALSE,
                      warning=FALSE,
                      message=FALSE,
                      dpi=300)
```

# 0. Prepare System
## Install System
BiocManager::install()

## Load System
```{r library}
library(stringr)
library(R.utils)
library(RColorBrewer)
library(sessioninfo)
library(data.table)
library(plyr)
library(tidyverse)
library(tximeta)
library(tximport)
library(curl)
library(AnnotationHub)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(writexl)
library(biomaRt)
library(ape)
library(kableExtra)
library(knitr)

library(stringr)
library(R.utils)
library(RColorBrewer)

library(sessioninfo)
library(data.table)
library(plyr)
library(tidyverse)
library(tximeta)
library(tximport)
library(curl)

library(SummarizedExperiment)
library(GenomicRanges)
library(ape)

library(viridis)
library(patchwork)
library("ggpubr")

# library(wget)

```
## Download data
Username: uFZjGr6t
Password: nYGTAdUw
Link https://filetransfer.mdc-berlin.de/?u=uFZjGr6t&p=nYGTAdUw

You can find detailed instruction accessing the folder in different settings here: https://t1p.de/gdmhowtocrush 

```{r folders}
getwd()
# setwd("S:/AG/AG-Scholz-NGS/Daten/Simon/P3043")
# wget("https://filetransfer.mdc-berlin.de/?u=uFZjGr6t&p=nYGTAdUw")

ifelse(Sys.info()["sysname"]== "Linux",
       s <- "/mnt/s",
       s <- "S:")
dir <- paste(s,"AG/AG-Scholz-NGS/Daten/Simon/P3043",sep="/")
gitdir <- paste(dir,"git_Chlamy_RNAseq_plap6",sep="/")
datadir <- paste(dir,"data",sep="/")
outdir <- gitdir
quantdir <- paste(dirname(datadir),"quants2",sep="/")

dir1 <- paste(dirname(dir),"linux-ngs/salmon/salmon_index_chlamy",sep="/")

indexDir <- file.path(dirname(dir),"linux-ngs/salmon/salmon_index_chlamy/chlamy_index_v6.1")
fastaPath <- file.path(dir1,"Phytozone_v6.1/CreinhardtiiCC_4532_707_v6.0.hardmasked.fa.gz")
# list.files(dirname(fastaPath))
# head(readLines(fastaPath,n=50))
gtfPath <- file.path(dir1,"Phytozone_v6.1/CreinhardtiiCC_4532_707_v6.1.repeatmasked_assembly_v6.0.gff3.gz")
# readLines(gtfPath,n=20)
gtfPath <- file.path(dir1,"Phytozone_v6.1/CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3.gz")
# readLines(gtfPath,n=20)
gtfPath <- file.path(dir1,"Phytozone_v6.1/CreinhardtiiCC_4532_707_v6.1.gene.gff3.gz")
# readLines(gtfPath,n=20)

file.exists(indexDir, fastaPath, gtfPath)

setwd(gitdir)

```
setwd("S:/AG/AG-Scholz-NGS/Daten/Simon/P3043")

# 1. Process Data
## Sample names
```{r sample_table}

f <- list.files(path = quantdir)
files <- as.factor(str_remove(f,pattern ="_quant"))
# short form
sample <- str_remove(files,pattern ="P3043_") 
sample %>% head() %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")

sampleid <- {}
strain  <- {}
lane  <- {}
condition <- {}
media  <- {}
replicate <- {}
trial <- {}

for (i in list.files(path = quantdir)){
  # print(i)
  si <- substring(i,14,18)
  sampleid <-c(sampleid,str_split(i, pattern = "_", simplify = TRUE)[,5])
  condition.i <- str_split(i, pattern = "_", simplify = TRUE)[,4]
  strain <-c(strain,str_split(condition.i, pattern = "-", simplify = TRUE)[,1])
  media <-c(media,str_split(condition.i, pattern = "-", simplify = TRUE)[,2])
  replicate <-c(replicate,str_split(condition.i, pattern = "-", simplify = TRUE)[,3])
  lane <-c(lane,str_split(i, pattern = "_", simplify = TRUE)[,6])
}
condition <- paste(strain,media,sep="_") %>% factor(levels = c("WT_TAP","fib_TAP","WT_HSM","fib_HSM"))
strain  <- strain %>% factor() %>% relevel("WT")
lane  <- lane %>% factor()
media  <- media %>% factor()  %>% relevel("TAP")
replicate <- replicate %>% factor()

trial <- c(rep(1,6),rep(2,8),rep(1,2),rep(2,4),rep(1,8),rep(1,12),rep(2,12),rep(1,12)) %>% factor() %>% relevel("1")

# Make table
sample.table <- data.frame(sample,sampleid,strain,media,condition,replicate,trial,lane)
sample.table %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
sample.table$strain %>% levels()

```


## Mapping Rates
https://stackoverflow.com/questions/14261776/extracting-data-from-text-files

working dir:
  "C:/Users/kelterbs/OneDrive - Charité - Universitätsmedizin Berlin/RNA-seq/Chlamy RNA-Seq P3043"

data dir:
  "S:/AG/AG-Scholz-NGS/Daten/Simon/P3043/quants2"


```{r mapping_rates}

# list.files(path = outdir)
# list.files(path = datadir)
# list.files(path = quantdir)

samplelist <- {}
mappingrates <- {}
for (i in list.files(path = quantdir)){
  # print(i)
  si <- paste("",str_sub(i,11,-7),sep = "")
  samplelist <- c(samplelist,si)
  f <- readLines(paste(quantdir,i,"logs/salmon_quant.log", sep="/"))
  line <- grep("Mapping rate = ",f,value=TRUE)
  sl <- str_length(line)
  notime <- substring(line,30,sl)
  manual <- substring(line,sl-7,sl-1)
  val <- as.numeric(str_extract(notime,"[0-9.]+"))
  valr<-round(val, digits=2)
  # print(paste("Mapping rate of ",si," is: ",valr," %"))
  mappingrates <- c(mappingrates,valr)
}

# Make table

m.table <- data.frame(sample.table,mappingrates)
m.table %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

### Plot mapping rates
```{r plot mappingrates, fig.height=6, out.width="100%"}

# Colours

# Plot
# plot(m.table$mappingrates)
# -> boring

# increase margin for longer names
# par(mar=c(4,15,4,4)+.1)
# barplot(height=mappingrates, names=m.table$sample, horiz=T, las=1)

par(mar=c(4,15,4,4)+.1)
xx <- barplot(height=mappingrates, names=m.table$sample, horiz=T, las=1, xlim=c(0,100)) #col=col
text(x = mappingrates, y = xx, label = mappingrates, pos = 4, cex = 0.8, col = "red")


```

## Tximeta

```{r tximeta}

example.quant <- paste(quantdir,list.files(quantdir)[1],"quant.sf",sep="/")

list.files(quantdir)[1]
file.exists(example.quant)

# first file:
example.table <- read.table(example.quant, header=T)
head(example.table)
basename(example.quant)
dirname(example.quant)
basename(dirname(example.quant))

# generate file list & prepare variables
files={}
for (i in list.dirs(path = dirname(dirname(example.quant)), full.names = TRUE, recursive = FALSE)) {
  files <- c(files,paste(i,"/quant.sf",sep=""))
  # print(basename(i))
  # print(head(read.table(files[length(files)], header=T)))
  # print (files[length(files)])
}
# head(files)

dir <- dirname(files)
# head(dir)
# head(basename(dir))
run <- basename(dir)
dirdir <- dirname(dir)
# head(dirdir)

tximeta_files <- file.path(dirdir, run, "quant.sf") 
file.exists(files,tximeta_files) %>% summary()

# names as working sample names
# see 'm.table' at mapping rates
coldata <- data.frame(files, names = run, stringsAsFactors=FALSE)
coldata <- data.frame(coldata,m.table)
coldata %>% head(10) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
colnames(coldata)
condition %>% levels()

# load tximeta
# with linked Transcriptome
## Chlamy

file.exists(indexDir, fastaPath, gtfPath)

makeLinkedTxome(indexDir=indexDir,
                source="Phytozome",
                organism="Chlamydomonas reinhardtii",
                release="v6.1",
                genome="CreinhardtiiCC_4532_707",
                fasta=fastaPath,
                gtf=gtfPath,
                write=FALSE)

gtfdata <- readLines(gtfPath)
# head(gtfdata, n=20)


# gene_name info sind in gff3 abder nicht in se?

se <- tximeta(coldata, useHub=TRUE)
se

colData(se) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")

# meta infos sind da (genotype, media,...)

# rename Samples
rownames(colData(se)) <- colData(se)$sample.id

rowRanges(se)
# genome info
seqinfo(se)
# ?
head(assays(se)[["counts"]])
# counts. THE DATA

# Mapping infos:
# names(metadata(se)[["quantInfo"]])
# str(metadata(se)[["quantInfo"]]) # Infos from Salmon Mapping
metadata(se)[["quantInfo"]]$percent_mapped
metadata(se)[["quantInfo"]]$num_processed

par(mfrow=c(1,2),
    mar=c(4,4,4,4)+.1)
barplot(metadata(se)[["quantInfo"]]$percent_mapped, main="Mapping Rate")
barplot(metadata(se)[["quantInfo"]]$num_processed/1000000, main="Mio. Reads")

edb <- retrieveDb(se)
class(edb)
genes(edb) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
columns(edb)


se.exons <- addExons(se)
rowRanges(se.exons)[[1]] %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
gse <- summarizeToGene(se)
rowRanges(gse) %>% head() %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")

head(assays(gse)[["counts"]])[1:5,1:5] %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
head(assays(gse)[["abundance"]])[1:5,1:5] %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
head(assays(gse)[["length"]])[1:5,1:5] %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")

#### add gene symbol
columns(retrieveDb(se))
edb
metadata(se)$txomeInfo$source

TXNAME <- as.character(mapIds(edb,keys = mcols(gse)$gene_id, column = "TXNAME", keytype = "GENEID", multiVals="first"))
# select()' returned 1:many mapping between keys and columns = 1 gene has many transcripts...
TXNAME.list <- (mapIds(edb,keys = mcols(gse)$gene_id, column = "TXNAME", keytype = "GENEID", multiVals="list"))
head(TXNAME.list)
# info is already there
CDSID <- as.character(mapIds(edb,keys = mcols(gse)$gene_id, column = "CDSID", keytype = "GENEID", multiVals="first"))
head(CDSID)
# ?
CDSNAME <- as.character(mapIds(edb,keys = mcols(gse)$gene_id, column = "CDSNAME", keytype = "GENEID", multiVals="first"))
head(CDSNAME)
# same as gene_id
EXONNAME <- as.character(mapIds(edb,keys = mcols(gse)$gene_id, column = "EXONNAME", keytype = "GENEID", multiVals="first"))
head(EXONNAME)
TXTYPE <- as.character(mapIds(edb,keys = mcols(gse)$gene_id, column = "TXTYPE", keytype = "GENEID", multiVals="first"))
head(TXTYPE)
# no useful info... no gene Symbol!!

# add CDSID
mcols(gse)$CDSID <- as.character(mapIds(edb,keys = mcols(gse)$gene_id, column = "CDSID", keytype = "GENEID", multiVals="first"))


colnames(mcols(gse))
head(rownames(gse))

getwd()
datadir
saveRDS(gse, file = paste(datadir,"gse.RDS",sep="/"), compress = FALSE)
gse <- {}
gse <- readRDS(paste(datadir,"gse.RDS",sep="/"))
gse

```

## Get Gene Names from GFF3 file
### a-I) Prepare gene.gff3 data
Use data from transcriptome used fore mapping: CreinhardtiiCC_4532_707_v6.1.gene.gff3.gz

```{r gene.gff3}
readLines(gtfPath,n=10)
gff <- read.delim(gtfPath, header=F, comment.char="#")
head(gff)

getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}

gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end","score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}

Chlamy_gff <- gffRead(gtfPath)
Chlamy_gff
Chlamy_gff$Name <- getAttributeField(Chlamy_gff$attributes, "Name")
Chlamy_gff$ID <- getAttributeField(Chlamy_gff$attributes, "ID")
# Chlamy_gff$Parent <- getAttributeField(Chlamy_gff$attributes, "Parent")
Chlamy_gff$pacid <- getAttributeField(Chlamy_gff$attributes, "pacid")
# Chlamy_gff$tID <- getAttributeField(Chlamy_gff$attributes, "tID")
Chlamy_gff$geneName <- getAttributeField(Chlamy_gff$attributes, "geneName")
# Chlamy_gff$longest <- getAttributeField(Chlamy_gff$attributes, "longest")
dim(Chlamy_gff)

# all IDs, some Gene symbols but no description etc...

# get only genes with Gene symbol
Chlamy_genes <- subset(Chlamy_gff, !geneName =="")
Chlamy_genes$gene_id <- paste(str_split(Chlamy_genes$Name, pattern = "_", simplify = TRUE)[,1])
gene2symbol <- data.frame(gene_id=Chlamy_genes$gene_id,
                          Symbol=Chlamy_genes$geneName)
dim(gene2symbol)
gene2symbol[duplicated(gene2symbol$gene_id),]
gene2symbol <- gene2symbol[!duplicated(gene2symbol$gene_id),]
dim(gene2symbol)
summary(mcols(gse)$gene_id %in% gene2symbol$gene_id)
# all FALSE (Cre01.g050800_4532)
summary(gene2symbol$gene_id %in% mcols(gse)$gene_id)
# 1 outliner

gene2symbol[!gene2symbol$gene_id %in% mcols(gse)$gene_id,] %>% dim()
# CreCp.g802322	icreI	

```


### a-II) add gene.gff3 data into gse
```{r add_gff3}
gse <- readRDS(paste(datadir,"gse.RDS",sep="/"))

# Remove "_4532"
mcols(gse)$gene_id2 <- mcols(gse)[,"gene_id"]
mcols(gse)$gene_id <- str_remove(mcols(gse)[,"gene_id2"], pattern = "_4532")
rownames(gse) <- mcols(gse)$gene_id

gsetable <- as.data.frame(mcols(gse)[,c("gene_id","gene_id2")])

dim(gsetable)
dim(gene2symbol)

summary(mcols(gse)$gene_id %in% gene2symbol$gene_id)
summary(gene2symbol$gene_id %in% mcols(gse)$gene_id)

new <- plyr::join(gsetable,gene2symbol)
dim(new)
head(new)
new2 <- left_join(gsetable,gene2symbol,by = "gene_id")
head(new2)

summary(new == new2)
# both methods worked!

mcols(gse)$symbol <- new$Symbol
mcols(gse)

getwd()
saveRDS(gse, file = paste(datadir,"gse.RDS",sep="/"), compress = FALSE)
gse <- readRDS(paste(datadir,"gse.RDS",sep="/"))
gse

```


### b-I) More gene data: master_annotation_table.tsv
Use other files from Phytozome database: CreinhardtiiCC_4532_707_v6.1.master_annotation_table.tsv

```{r master_annotation}

# more gene info from master_annotation_table
annopath <- file.path(dir1,"Phytozone_v6.1/CreinhardtiiCC_4532_707_v6.1.master_annotation_table.tsv")

anno <- read.delim2(annopath, header = T, sep = "\t",)
head(anno) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")

# search for genes
anno[str_detect(anno[["locusName_4532"]],"Cre01.g000050"),]
anno[str_detect(anno[["geneSymbol"]],"HKR"),]
anno[str_detect(anno[["previousIdentifiers"]],"HKR"),]
anno[str_detect(anno[["geneSymbol"]],"COP"),]
anno[str_detect(anno[["previousIdentifiers"]],"COP1"),]
anno[str_detect(anno[["geneSymbol"]],"cop"),]
anno[str_detect(anno[["geneSymbol"]],"PLAP6"),]
anno[str_detect(anno[["geneSymbol"]],"COQ3"),]
anno[str_detect(anno[["geneSymbol"]],"PCRY"),]
anno[str_detect(anno[["geneSymbol"]],"CRY"),]
anno[str_detect(anno[["geneSymbol"]],"CRY-DASH"),]
anno[str_detect(anno[["previousIdentifiers"]],"CRY-DASH"),]

dim(anno)

# Remove "_4532"
anno$gene_id <- str_remove(anno$locusName_4532,pattern = "_4532")
rownames(anno) <- anno$gene_id

# convert previousIdentifiers to list
anno$previousIdentifiers_list <- strsplit(anno$previousIdentifiers, "#")
head(anno) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
anno$previousIdentifiers_list[1]

# polish "previousIdentifiers"
# remove Cre & g-numbers
head(anno$previousIdentifiers)
tail(anno$previousIdentifiers)
prev.symbols <-  
  str_replace_all(anno$previousIdentifiers,pattern="Cre\\d+.g\\d+.t\\d+.\\d+", replacement="?") %>% 
  str_replace_all(pattern="Cre...g\\d+", replacement="?") %>% 
  str_replace_all(pattern="g\\d+.t\\d+", replacement="?") %>% 
  str_remove_all(pattern="\\?") %>%
  str_remove_all("^#|#$| ") %>% str_remove_all("^#|#$") %>% str_remove_all("^#|#$|^0") %>% 
  str_replace_all(pattern="##","#") %>% str_replace_all(pattern="##","#")
prev.symbols[str_detect(prev.symbols,"#")] %>% head(n=20)
prev.symbols[str_detect(prev.symbols,"0")] %>% head(n=20)
prev.symbols[str_detect(prev.symbols,"Cre")] %>% head(n=20)
tail(anno$previousIdentifiers[str_detect(anno$previousIdentifiers,"Cre")])


anno$prev.symbols <- prev.symbols
head(anno, n=10) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")

# Make colum with gene symbol or prev.symbol or Gene id
anno$combname <- anno$geneSymbol

summary(anno$combname)
summary(anno$combname=="")
summary(anno$combname==" ")
summary(is.na(anno$combname))
summary(anno$combname=="0")

# replace "" with prev.symbol
summary(anno$combname=="")
length(anno$combname[anno$combname==""])
length(str_split(
  anno$prev.symbols[anno$combname==""],pattern ="#", simplify = T)[,1])
anno$combname[anno$combname==""] <- str_split(
  anno$prev.symbols[anno$combname==""],pattern ="#", simplify = T)[,1]
summary(anno$combname=="")

anno$combname[anno$combname==""] <- anno$gene_id[anno$combname==""]
summary(anno$combname=="")
subset(anno,anno$combname=="")

summary(duplicated(anno$combname))
# some duplicates

colnames(anno)[which(names(anno) == "combname")] <- "id.symbol"


summary(anno$geneSymbol!="")
# 4408 gene Symbols
summary(anno$previousIdentifiers!="")
# all (!) 17712 genes have "previousIdentifiers"
summary(anno$prev.symbols!="")
# 7236 have alternative/old symbols
summary(anno$id.symbol)

dim(anno)

saveRDS(anno, file = paste(datadir,"anno.RDS",sep="/"), compress = FALSE)
anno <- readRDS(paste(datadir,"anno.RDS",sep="/"))

```


### b-II) add anno data into gse ####
```{r add_anno}
gse <- readRDS(paste(datadir,"gse.RDS",sep="/"))
anno <- readRDS(paste(datadir,"anno.RDS",sep="/"))

summary(anno$gene_id %in% mcols(gse)$gene_id)
summary(mcols(gse)$gene_id %in% anno$gene_id)

# choose info
anno_join <- anno[,c("gene_id", "geneSymbol", "id.symbol","prev.symbols", "previousIdentifiers", "previousIdentifiers_list", "Description", "Comments", "TargetP", "Predalgo")]
anno_join

class(gsetable)
class(anno_join)

dim(gsetable)
dim(anno_join)

# any duplicated:
summary(duplicated(anno$gene_id))

new <- plyr::join(gsetable,anno_join, by = "gene_id", type = "left")
dim(new)
head(new)
new2 <- left_join(gsetable,anno_join,by = "gene_id")
dim(new2)
head(new2)

mcols(gse) <- left_join(as.data.frame(mcols(gse)),anno_join,by = "gene_id")
mcols(gse)

getwd()

saveRDS(gse, file = paste(datadir,"gse.RDS",sep="/"), compress = FALSE)
gse <- readRDS(paste(datadir,"gse.RDS",sep="/"))
gse

```



## DESeq2 Analysis

```{r desq2,fig.show="hold"}

gse <- readRDS(paste(datadir,"gse.RDS",sep="/"))

colData(gse) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
head(mcols(gse)) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
gse$condition %>% levels()

# dds <- DESeqDataSet(gse, ~media+strain+media:genotype)
dds <- DESeqDataSet(gse, ~condition)
colData(dds)$names %>% head()

##########################################\n
## filter all rows with rowsum = 0 ####

hist(log(counts(dds)), breaks=100, ylim = c(0,100000))
hist(counts(dds), breaks=1000000, ylim = c(0,100000), xlim = c(0,100))
# -> dip at log(counts(dds)=2)

# remove 0 counts
remove0 <- rowSums(counts(dds) == 0) == length(colData(dds)[, 1])
summary(remove0)
dds <- dds[!remove0,]

# min 16 counts in 16 samples
sample.number <- nrow(colData(dds)) / colData(gse)$condition %>% levels() %>% length()
keep.sn <- rowSums(counts(dds) >= sample.number) >= sample.number
summary(keep.sn)
dds <- dds[keep.sn,]

dds
head(assays(dds)[["counts"]])[1:5,1:5] %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")

length(rownames(counts(dds)))

colData(dds) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")

## combine technical replicates #####
dds <- collapseReplicates(dds, dds$sampleid,dds$lane)
colData(dds) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
dds$runsCollapsed
colnames(dds)
matchFirstLevel <- dds$sampleid == levels(factor(dds$sampleid))[1]
all(rowSums(counts(dds[,matchFirstLevel])) == counts(dds[,1]))

# dds <- estimateSizeFactors(dds)
colData(dds)$genotype <- colData(dds)$strain


###########
# run DESeq
dds <- DESeq(dds)
plotMA(dds)


plotCounts(dds, gene = "Cre01.g000150", intgroup = "condition",col="trail")


###############
# save dds #####

saveRDS(dds, file = paste(datadir,"dds.RDS",sep="/"), compress = FALSE)
dds <- readRDS(paste(datadir,"dds.RDS",sep="/"))


```


## Quality I
### - Data transformations
```{r pre_trans}
dds <- readRDS(paste(datadir,"dds.RDS",sep="/"))
vsd <- vst(dds, blind=FALSE) #Variance stabilized transformation
ntd <- normTransform(dds)
```

```{r pre_trans_rlog, eval=FALSE}
rld <- rlog(dds, blind=FALSE) #regularized logarithm
save(rld,file=paste(data,"rlog.rld", sep="/"))
rld <- 1
load(file=paste(data,"rlog.rld", sep="/"))
rld
```

```{r pre_trans_fig, figures-side, fig.show="hold", out.width="33%"}
load(file=paste(data,"rlog.rld", sep="/"))
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

```

### - Check sample distance
```{r pre_sample_dist, fig.height=12, out.width="100%"}
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$names
# colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
anno_col <- as.data.frame(colData(vsd)[,c("experiment","treatment","genotype")])
anno_colors <- list(media = c("lightcoral","skyblue1"),
                    genotype = c("grey","seagreen3","turquoise3","tan2"),
                    trial = viridis(2, option="plasma"))

names(anno_colors$treatment) <- levels(anno_col$treatment)
names(anno_colors$genotype) <- levels(anno_col$genotype)
names(anno_colors$experiment) <- levels(anno_col$experiment)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation_col=anno_col,
         annotation_colors = anno_colors,
         show_colnames     = FALSE,
         col=viridis(20),
         cutree_rows = 8,
         cutree_cols = 8)


```




### Rename gene.ids 

```{r eval=FALSE, include=FALSE}
################
# outliner: S15, S16, S17, S18, S19, S26
outliner <- c("S15", "S16", "S17", "S18", "S19", "S26")
rmo <- colnames(dds)[str_detect(colnames(dds), paste(outliner, collapse = "|"))]
keep <- colnames(dds)[!str_detect(colnames(dds), paste(outliner, collapse = "|"))]
keep

dds[,keep]
dds <- dds[,keep] 

############



###########
# run DESeq
dds <- DESeq(dds)
plotMA(dds)

plotCounts(dds, gene = "Cre01.g000150", intgroup = "condition")

head(mcols(dds)$symbol)

anno[str_detect(anno[["geneSymbol"]],"CRY-DASH"),]
# no hit
g <- anno[str_detect(anno[["previousIdentifiers"]],"CRY-DASH1"),"gene_id"]
g
plotCounts(dds, gene = g, intgroup = "condition")
g <- anno[str_detect(anno[["geneSymbol"]],"PCRY"),"gene_id"]
plotCounts(dds, gene = g, intgroup = "condition", col=colData(dds)$genotype, main =anno[g,"geneSymbol"])
g <- anno[str_detect(anno[["geneSymbol"]],"ACRY"),"gene_id"]
plotCounts(dds, gene = g, intgroup = "condition", col=colData(dds)$genotype, main =anno[g,"geneSymbol"])
g <- anno[str_detect(anno[["geneSymbol"]],"ROC15"),"gene_id"]
plotCounts(dds, gene = g, intgroup = "condition", col=colData(dds)$genotype, main =anno[g,"geneSymbol"])
g <- anno[str_detect(anno[["geneSymbol"]],"ROC40"),"gene_id"]
plotCounts(dds, gene = g, intgroup = "condition", col=colData(dds)$genotype, main =anno[g,"geneSymbol"])


# search genes
anno[str_detect(anno[["geneSymbol"]],"CRY"),]
anno[str_detect(anno[["geneSymbol"]],"ROC"),]
anno[str_detect(anno[["previousIdentifiers"]],"CRY"),]


###############
# save dds #####

saveRDS(dds, file = paste(datadir,"dds.RDS",sep="/"), compress = FALSE)
dds <- readRDS(paste(datadir,"dds.RDS",sep="/"))
dds

```




### Rename gene.ids 

```{r eval=FALSE, include=FALSE}
# aleady done
# head(rownames(dds))
# gene_ids <- rownames(dds) %>% str_remove(pattern = "_4532")
# head(gene_ids)
# length(gene_ids)
# length(rownames(dds))
# rownames(dds) <- gene_ids
# plotCounts(dds, gene = "Cre01.g000150", intgroup = "condition")
# PCRY <- "Cre06.g295200"
# plotCounts(dds, gene = PCRY, intgroup = "condition")
# 
# saveRDS(dds, file = paste(datadir,"dds.RDS",sep="/"), compress = FALSE)
# dds <- readRDS(paste(datadir,"dds.RDS",sep="/"))
# dds

```



## Quality

### PCA
```{r}
datadir
dds <- readRDS(paste(datadir,"dds.RDS",sep="/"))

colData(dds)
vsd <- vst(dds, blind=FALSE)
z <- plotPCA(vsd, intgroup = c("condition"))
z$data$name <- colData(dds)$sample 
z + geom_label(aes(label = name))

mappingrates <- colData(dds)[,c("sampleid","mappingrates")]
plot(mappingrates)
names(colData(vsd))
pcaData <- plotPCA(vsd, intgroup=colnames(colData(vsd)), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

g1 <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape = trial)) +
  geom_point(size=3) +
  scale_color_viridis_d() +
  geom_text_repel(aes(label = name), size = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) # + 
# coord_fixed()

g2 <- ggplot(pcaData, aes(PC1, PC2, color=mappingrates, shape = condition)) +
  geom_point(size=3) +
  scale_color_viridis_c() +
  geom_text_repel(aes(label = name), size = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

g3 <- ggplot(pcaData, aes(PC1, PC2, color=trial, shape = condition)) +
  geom_point(size=5) +
  scale_color_viridis_d() +
  geom_text_repel(aes(label = name), size = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

g1+g2+g3+plot_layout(guides = "collect",axis = "collect")

g4 <- ggplot(pcaData, aes(PC1, PC2, color=mappingrates, shape = genotype)) +
  geom_point(size=3) +
  scale_color_viridis_c() +
  geom_text_repel(aes(label = name), size = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

# run all lines together
pdf(file="graphs3/PCA.pdf", width=5, height=5)
par(mar=c(4,4,4,4)+.1)
plot(g)
dev.off()

ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))



```






## Make results 

```{r}
resultsNames(dds)

# TAP = fib_TAP - WT_TAP
# HSM = fib_HSM - WT_HSM
# fib = fib_HSM - fib_TAP
# WT = WT_HSM - WT_TAP

anno

res_fib.TAP.vs.WT.TAP <-
  results(dds, contrast = c("condition","fib_TAP","WT_TAP"))
res_fib.TAP.vs.WT.TAP$Symbol <- mcols(dds)$id.symbol

res_fib.HSM.vs.WT.HSM <-
  results(dds, contrast = c("condition","fib_HSM","WT_HSM"))
res_fib.HSM.vs.WT.HSM$Symbol <- mcols(dds)$id.symbol

res_fib.HSM.vs.fib.TAP <-
  results(dds, contrast = c("condition","fib_HSM","fib_TAP"))
res_fib.HSM.vs.fib.TAP$Symbol <- mcols(dds)$id.symbol

res_WT.HSM.vs.WT.TAP <-
  results(dds, contrast = c("condition","WT_HSM","WT_TAP"))
res_WT.HSM.vs.WT.TAP$Symbol <- mcols(dds)$id.symbol

res1 <- res_fib.TAP.vs.WT.TAP
res2 <- res_fib.HSM.vs.WT.HSM
res3 <- res_fib.HSM.vs.fib.TAP
res4 <- res_WT.HSM.vs.WT.TAP

summary(res1)
top_res1 <- subset(res1, padj < 0.01 & baseMean > 50 &
                     (log2FoldChange < -1 | log2FoldChange > 1))
top_res1 <- top_res1[order(top_res1$log2FoldChange, decreasing = T),]
dim(top_res1)

summary(res2)
top_res2 <- subset(res2, padj < 0.01 & baseMean > 50 &
                     (log2FoldChange < -1 | log2FoldChange > 1))
top_res2 <- top_res2[order(top_res2$log2FoldChange, decreasing = T),]
dim(top_res2)

summary(res3)
top_res3 <- subset(res3, padj < 0.01 & baseMean > 50 &
                     (log2FoldChange < -2 | log2FoldChange > 2))
top_res3 <- top_res3[order(top_res3$log2FoldChange, decreasing = T),]
dim(top_res3)

summary(res4)
top_res4 <- subset(res4, padj < 0.01 & baseMean > 50 &
                     (log2FoldChange < -2 | log2FoldChange > 2))
top_res4 <- top_res4[order(top_res4$log2FoldChange, decreasing = T),]
dim(top_res4)

resLFC_TAP <- lfcShrink(dds, contrast = c("condition","fib_TAP","WT_TAP"), type="ashr")
resLFC_TAP$Symbol <- mcols(dds)$id.symbol
resLFC_HSM <- lfcShrink(dds, contrast = c("condition","fib_HSM","WT_HSM"), type="ashr")
resLFC_HSM$Symbol <- mcols(dds)$id.symbol
resLFC_fib <- lfcShrink(dds, contrast = c("condition","fib_HSM","fib_TAP"), type="ashr")
resLFC_fib$Symbol <- mcols(dds)$id.symbol
resLFC_WT <- lfcShrink(dds, contrast = c("condition","WT_HSM","WT_TAP"), type="ashr")
resLFC_WT$Symbol <- mcols(dds)$id.symbol

res1["Cre01.g000150",]
plotMA(resLFC_TAP, main = "fib vs. WT in TAP", ylim=c(-4,4))
summary(res1)

pdf(file="graphs3/MA-fib.vs.WT(TAP).pdf", width=8, height=5)
plotMA(resLFC_TAP, main = "fib vs. WT in TAP", ylim=c(-4,4))
dev.off()

res2["Cre01.g000150",]

```


# 2.) Data Dive

## Get gen names

```{r}
# Load data
outdir <- "C:/Users/kelterbs/OneDrive - Charité - Universitätsmedizin Berlin/NGS (RNA-Seq & CRISPR)/Chlamy RNA-Seq P3043"
setwd(outdir)

datadir <- "S:/AG/AG-Scholz-NGS/Daten/Simon/P3043/"
dds <- readRDS(file=paste(datadir,"dds.RDS", sep="/"))
anno <- readRDS(file=paste(datadir,"anno.RDS", sep="/"))
list.files(outdir)

colData(dds)$genotype

COQ3 <- "Cre10.g423750"
PHO1 <- "Cre08.g359300"

# search genes
anno[str_detect(anno[["geneSymbol"]],"COQ3"),]
# get only gene IDs
anno[str_detect(anno[["geneSymbol"]],"CRY"),"gene_id"]

# gene ID
# g245450
anno[str_detect(anno[["gene_id"]],"g245450"),]


# More searches:
anno[str_detect(anno[["geneSymbol"]],"ROC"),]
anno[str_detect(anno[["prev.symbols"]],"CRY"),]
anno[str_detect(anno[["Description"]],paste(c('Opsin', 'opsin','photoreceptor'), collapse="|")),]
anno[str_detect(anno[["Comments"]],paste(c('phototaxis'), collapse="|")),]

# all Chloroplast
anno %>% filter_at(.vars = vars(Predalgo, TargetP),
                   .vars_predicate = any_vars(str_detect
                                              (. , paste0("^(", paste("Chloroplast", collapse = "|"), ")"))))
# all Flagellar
anno[str_detect(anno[["Flagellar_Proteome"]],"Total Peptides:"),]

# all transmembrane with 7 helices
anno[str_detect(anno[["TMHMM_transmembrane"]],"TMHMM: 7 helices"),]

# fibrillin related
coqs <- anno[str_detect(anno[["geneSymbol"]],"COQ"),]
anno[str_detect(anno[["Comments"]],paste(c('Coenzyme Q'), collapse="|")),]
goi <- anno[c(PHO1,coqs[,"gene_id"]),]

```

## Plot Counts
```{r}


gene <- COQ3
plotCounts(dds, gene = gene, intgroup = "condition", col=colData(dds)$genotype, main =anno[gene,"geneSymbol"])

# with ggplot
d <- plotCounts(dds, gene = gene, intgroup = c("condition","media","genotype"), col=colData(dds)$genotype, main =anno[gene,"geneSymbol"], returnData=TRUE)
d
# extract data and save as table

g1 <- ggplot(d, aes(x = condition, y = count, color = genotype)) + 
  geom_point() +
  ggtitle(anno[gene,"geneSymbol"])
plot(g1)

# more advanced:
g1 <- ggplot(d, aes(x = media, y = count, color = genotype)) + 
  geom_boxplot(aes(group = condition, colour = genotype)) +
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  # geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle(anno[gene,"geneSymbol"]) +
  theme(plot.title = element_text(hjust = 0.5))
plot(g1)

ggsave(paste("graphs3/Counts_",anno[gene,"geneSymbol"],".pdf"),
       width = 5,
       height = 5)

## multiple genes
# all ROC genes
goi

n <- {0}
for (i in goi$gene_id){
  n <- n+1
  print(n)
  print(i)
  s <- print(anno[i,"geneSymbol"])
  d <-  plotCounts(dds, gene=i, intgroup=c("condition","media","genotype"),main=s,returnData=TRUE)
  g <- ggplot(d, aes(x = condition, y = count, color = genotype)) + 
    geom_boxplot(aes(group = condition, colour = genotype)) +
    geom_point() + #position=position_jitter(w = 0.1,h = 0)
    # geom_text_repel(aes(label = rownames(d))) + 
    theme_bw() +
    ggtitle(s) +
    theme(plot.title = element_text(hjust = 0.5))
  assign(x = paste("g",n, sep=""),g)
  plot(g)
}
length(goi)
ga <- ggarrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,ncol = 3, nrow = 4)
ga
plot(ga)
ggexport(ga, filename = "graphs3/Counts_COQs.pdf",width = 12, height = 12)


```


## Volcano Plot

```{r}
library(EnhancedVolcano)

plot(res1$log2FoldChange,res1$baseMean )

# colnames(mcols(dds))
# mcols(dds) <- left_join(as.data.frame(mcols(dds)),anno[,c("gene_id","prev.symbols","id.symbol")],by = "gene_id")

summary(res1)

EnhancedVolcano(res1,
                x = 'log2FoldChange',
                y = 'pvalue',
                lab = mcols(dds)$id.symbol,
                # boxedLabels = TRUE,
                labSize = 2.0,
                drawConnectors = TRUE,
                max.overlaps = 20,
                # xlim = c(-10, 10),
                # ylim = c(0, 50),
                ylab = "Padj (-log10)",
                title = 'fib vs. WT in TAP',
                subtitle = "DE genes: ?? (red)",
                # sub = "SVF",
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 2.0,
                legendLabels=c('Not sig.','|L2F| > 2','p-adj < 0.01',
                               'p-adj & L2F'),
                legendPosition = 'right',
)

ggsave("graphs3/EnhancedVolcano.pdf",
       width = 12,
       height = 10)

# with ggplot
ggplot(as.data.frame(res1),aes(log2FoldChange,-log(padj))) +
  geom_point() +
  geom_point(data = as.data.frame(res1[rownames(top_res1),]), colour = "blue") + 
  geom_point(data = as.data.frame(res1[goi$gene_id,]), colour = "red") +
  geom_text_repel(data = as.data.frame(res1[goi$gene_id,]), aes(label=anno[goi$gene_id,"id.symbol"]),colour = "red",hjust=-0.5, vjust=-1) + xlim(-5,5)

getwd()
ggsave("graphs3/Vulcano COQs-2.pdf",
       width = 10,
       height = 6)

```


## Heatmap
```{r}
library("pheatmap")
ntd <- normTransform(dds)

# select top 100 highest expressed genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("condition","media","genotype")])
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

# all top genes
select <- unique(
  rownames(top_res1),
  rownames(top_res2)) %>% 
  unique(rownames(top_res3)) %>%
  unique(rownames(top_res4)) %>%
  unique(rownames(goi))

length(select)
df <- assay(ntd)[select,]
rownames(df) <- mcols(dds)[select,"id.symbol"]

anno_col <- as.data.frame(colData(dds)[,c("condition","media","genotype")])
xx <- pheatmap(df, cluster_rows=TRUE, show_rownames=TRUE,
               cluster_cols=TRUE, annotation_col=anno_col)
ggsave("graphs3/Heatmap_top.pdf",plot=xx,
       width = 10,
       height = 30)

# all fib genes
select <- unique(c(
  rownames(top_res1),
  rownames(top_res2),
  rownames(goi)))

select <- unique(c(head(rownames(top_res1),n=50),
                   tail(rownames(top_res1),n=50),
                   head(rownames(top_res2),n=50),
                   tail(rownames(top_res2),n=50),
                   rownames(goi)))

length(select)

df <- assay(ntd)[select,]
rownames(df) <- mcols(dds)[select,"id.symbol"]

anno_col <- as.data.frame(colData(dds)[,c("condition","media","genotype")])
xx <- pheatmap(df, cluster_rows=TRUE, show_rownames=TRUE,
               cluster_cols=TRUE, annotation_col=anno_col)
ggsave("graphs3/Heatmap_fib.vs.WT.pdf",plot=xx,
       width = 10,
       height = 30)

```
## Genes Yousef 
```{r}
#
tocopherol <- c("Cre09.g414000", "Cre01.g013801", "Cre14.g624350", "Cre09.g393400")
plastoquinone <- c("Cre09.g398993", "Cre06.g283750", "Cre14.g625450", "Cre12.g503550")
phylloquinone <- c("Cre10.g455950", "Cre16.g671000", "Cre04.g219787", "Cre09.g416500")
goi_y <- c(tocopherol,plastoquinone, phylloquinone)
goi_y

```

### Plot counts 
```{r}
library("ggpubr")

## multiple genes
goi <- anno[goi_y,]
dim(goi)

n <- {0}
for (i in goi$gene_id){
  n <- n+1
  print(n)
  print(i)
  s <- print(anno[i,"id.symbol"])
  print
  d <-  plotCounts(dds, gene=i, intgroup=c("condition","media","genotype"),main=s,returnData=TRUE)
  g <- ggplot(d, aes(x = condition, y = count, color = genotype)) + 
    geom_boxplot(aes(group = condition, colour = genotype)) +
    geom_point() + #position=position_jitter(w = 0.1,h = 0)
    # geom_text_repel(aes(label = rownames(d))) + 
    theme_bw() +
    ggtitle(s) +
    theme(plot.title = element_text(hjust = 0.5))
  assign(x = paste("g",n, sep=""),g)
  plot(g)
}
length(goi)
ga <- cowplot::plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12, byrow = FALSE, nrow = 4, ncol = 3)

ga
plot(ga)
ggexport(ga, filename = "graphs3/Counts_YYK.pdf",width = 12, height = 12)



```
### Heatmap

```{r}
library("pheatmap")
ntd <- normTransform(dds)

df <- assay(ntd)[goi$gene_id,]
rownames(df) <- mcols(dds)[goi$gene_id,"id.symbol"]

anno_col <- as.data.frame(colData(dds)[,c("condition","media","genotype")])
xx <- pheatmap(df, cluster_rows=TRUE, show_rownames=TRUE,
               cluster_cols=TRUE, annotation_col=anno_col)

ggsave("graphs3/Heatmap_YYK.pdf",plot=xx,
       width = 10,
       height = 10)

```



# Export

```{r}
library(writexl)

# combine results

colnames(anno)

res_exp <- data.frame(
  # general
  "gene_id" = rownames(dds),
  "gene_name" = anno[rownames(dds),"geneSymbol"],
  "Alias" = anno[rownames(dds),"prev.symbols"],
  "id.symbol" = anno[rownames(dds),"id.symbol"],
  "Description" = anno[rownames(dds),"Description"],
  "Comments" = anno[rownames(dds),"Comments"],
  "TargetP" = anno[rownames(dds),"TargetP"],
  "Predalgo" = anno[rownames(dds),"Predalgo"],
  "baseMean" = res1$baseMean,
  # results 1
  "L2FC.fib_TAP.vs.WT_TAP" = res1$log2FoldChange,
  "pvalue.fib_TAP.vs.WT_TAP" = res1$pvalue,
  "padj.fib_TAP.vs.WT_TAP" = res1$padj,
  # results 2
  "L2FC.fib_HSM.vs.WT_HSM" = res2$log2FoldChange,
  "pvalue.fib_HSM.vs.WT_HSM" = res2$pvalue,
  "padj.fib_HSM.vs.WT_HSM" = res2$padj,
  # results 3
  "L2FC.WT_HSM.vs.WT_TAP" = res3$log2FoldChange,
  "pvalue.WT_HSM.vs.WT_TAP" = res3$pvalue,
  "padj.WT_HSM.vs.WT_TAP" = res3$padj,
  # results 4
  "L2FC.fib_HSM.vs.fib_TAP" = res4$log2FoldChange,
  "pvalue.fib_HSM.vs.fib_TAP" = res4$pvalue,
  "padj.fib_HSM.vs.fib_TAP" = res4$padj,
  counts(dds, normalized = TRUE))
res_exp

outdir <- "C:/Users/kelterbs/OneDrive - Charité - Universitätsmedizin Berlin/NGS (RNA-Seq & CRISPR)/Chlamy RNA-Seq P3043"

write_xlsx(data.frame(res_exp),
           paste(outdir,"2023_08 P3043 results.xlsx",sep="/"))

write_xlsx(data.frame(colData(dds)),
           paste(outdir,"2023_08 P3043 samples.xlsx",sep="/"))


```


# End
```{r}

```



