
cd /data/sn/HErik/testosteron_220403/FASTQ_Generation_2022-04-03_20_27_02Z-548161623

mv *ds*/*.gz .
rm -rf *ds*

export PATH=$PATH:/data/tools/FastQC

mkdir qc01

for f in *fastq.gz 
do 
  fastqc -t 38 $f -o qc01
done 

multiqc qc01 . 


for f in *_L001_R1_001.fastq.gz
do
  zcat ${f/'_L001_R1_001.fastq.gz'/'*'} > ../${f/'_L001_R1_001.fastq.gz'/'.fastq'}
done

TRIM='/data/tools/Trimmomatic-0.38/trimmomatic-0.38.jar'

for f in *.fastq
do 
  o=${f/'.fastq'/'_trimmed.fastq'}
  java -jar $TRIM SE -threads 38 \
    $f $o \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:30 \
    MINLEN:50
done

export PATH=/data/tools/STAR-2.7.9a/bin/Linux_x86_64
export PATH=/data/tools/subread-2.0.2-Linux-x86_64/bin:$PATH


idx='/data/dbs/STAR/GRCm39_104'

for f in *_trimmed.fastq
do
  root=${f/'.fastq'/''}
  STAR -- genomeDir $idx \
       -- readFilesIn $f \
       -- outFileNamePrefix 'GRCm39_104_'$root \
       -- outFilterMultimapNmax 1 \
       -- outReadsUnmapped Fastx \
       -- outSAMtype BAM SortedByCoordinate \
       -- twopassMode Basic \
       -- runThreadN 14 
done
 
featureCounts -O \
   -a $idx/Mus_musculus.GRCm39.104.gtf \
   -o featureCounts_GRCm39.104_O.txt \
   GRCm39_104*bam 


# R

library(org.Mm.eg.db)
library(GO.db)
library(KEGGREST)
library(edgeR)
library(DESeq2)
library(openxlsx)
library(tidyverse)

setwd('/data/sn/HErik/set04/KOKI_210701-276856580')

org.db = org.Mm.eg.db

uniprot.db_sel = read_tsv('MOUSE_10090_idmapping_selected.tab', col_names=c('UniProtKB_AC', 'UniProtKB_ID', 'GeneID_EntrezGene', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI_taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL', 'EMBL_CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional_PubMed')) 

PTHR = read_tsv(
    'PTHR16.0_mouse', 
    col_names=c('GeneIdentifier', 'ProteinID', 'SFID', 'FamilyName', 'SubfamilyName', 'MolecularFunction', 'BiologicalProcess', 'CellularComponents', 'ProteinClass', 'Pathway'), 
    quote='') %>% 
  separate(GeneIdentifier, c('GeneIdentifier', 'UniProt'), 'UniProtKB=')
    
ens2uniprot = read_tsv('Mus_musculus.GRCm39.104.uniprot.tsv', col_names=T) %>%
  dplyr::rename(UniProt=xref, ens=gene_stable_id) %>%
  dplyr::select(ens, UniProt) %>% 
  unique()

go = read_tsv('mgi.gaf', skip=43, col_names=F, 
  col_types = cols(.default = "c")) %>%
  dplyr::select(2,5,7,9,10,11) %>%
  mutate(X9=case_when(X9=='C' ~ 'CC', X9=='P' ~ 'BP', X9=='F' ~ 'MF')) %>%
  dplyr::rename(UniProt=1, GOid=2, evidence=3, ontology=4, GeneName=5, symbol=6)  

ens2entrez = read_tsv('Mus_musculus.GRCm39.104.entrez.tsv', 
  col_types = cols(.default = "c")) 

 
setwd('/data/sn/HErik/testosteron_220403')
  
read.counts = read_delim('featureCounts_GRCm39.104_O.txt', delim='\t', col_names=T, skip=1) %>% 
  rename_at(vars(contains('GRCm39_104_')), list( ~ gsub('GRCm39_104_', '', .))) %>%  
  rename_at(vars(contains('_trimmedAligned.sortedByCoord.out.bam')), list( ~ gsub('_trimmedAligned.sortedByCoord.out.bam', '', .))) %>%  
  rename_at(vars(contains('-')), list( ~ gsub('-', '_', .)))

readcounts = read.counts %>% 
  dplyr::select(grep('_', tolower(colnames(.)))) %>% 
  as.data.frame()
rownames(readcounts) = read.counts$Geneid

meta = data.frame(smpl=colnames(readcounts))
rownames(meta) = meta$smpl
# meta$rep = 1:12
# meta$grp = 'a'
     
bres = tibble(ens=rownames(readcounts))

dds = DESeqDataSetFromMatrix(
    countData = readcounts, 
    colData = meta, 
    design = ~ 1 
)

dds = estimateSizeFactors(dds)

counts_raw = readcounts
colnames(counts_raw) = paste0('raw_', colnames(counts_raw))
counts_raw$ens = rownames(counts_raw)
counts_raw = as_tibble(counts_raw)

counts_normalized = as.data.frame(counts(dds, normalized=T))
colnames(counts_normalized) = paste0('norm_', colnames(counts_normalized))
counts_normalized$ens = rownames(counts_normalized)
counts_normalized = as_tibble(counts_normalized)

m = as.matrix(readcounts)
counts_cpm = as.data.frame(cpm(m))
colnames(counts_cpm) = paste0('cpm_', colnames(counts_cpm))
counts_cpm$ens = rownames(counts_cpm)
counts_cpm = as_tibble(counts_cpm)

# ids = rownames(readcounts)
# n = 1
# tmp = tibble(.rows=0, ens='', UniProt='')
# for(ens in ids){
#   tmp = rbind(tmp, 
#     tibble(ens, UniProt= uniprot.db_sel %>% 
#     filter(str_detect(Ensembl, ens)) %>%
#     pull(UniProtKB_AC)  
#     )
#   )
#   n=n+1
#   print(n)
# }
# save(tmp, file='ens_tmp_104.RData')

load('ens_tmp_104.RData')

ens2uniprot = rbind(ens2uniprot, tmp) %>% 
  unique()

ens_uniprot_pthr = inner_join(ens2uniprot, PTHR)

tib = left_join(tibble(ens=rownames(readcounts)), ens_uniprot_pthr) %>% 
  select(ens, UniProt)
  
tib = inner_join(tib, counts_raw)
tib = inner_join(tib, counts_normalized)
tib = inner_join(tib, counts_cpm)
  
annot = left_join(
  tib,
  PTHR %>% 
    select(UniProt, FamilyName, SubfamilyName, ProteinClass, BiologicalProcess, CellularComponents, MolecularFunction)
)

evs = sort(unique(go$evidence))
onts = sort(unique(go$ontology))

for(ont in onts){
  for(ev in evs){
    tmp = go %>% 
      filter(ontology==ont, evidence==ev) %>% 
      select(UniProt, GeneName)
    if(dim(tmp)[1]>0){ 
      lst = split(tmp$GeneName, tmp$UniProt) %>% 
        lapply(unique) %>% 
        lapply(sort) %>%  
        lapply(paste, collapse='\n')   
      annot = left_join(
        annot, 
        tibble(UniProt=names(lst), tmp=as.character(lst)) %>%
          rename_at(vars(tmp), ~ paste(ont, ev, sep='_'))
        )
    }
  }
}

lst = keggList('pathway', 'mmu') 
paths = tibble(
  PATH=gsub('path:mmu', '', names(lst)), 
  pathway=gsub(' - Mus musculus \\(mouse\\)', '', as.character(lst))
) 

i = 1
query = keggGet(paste0('mmu', paths$PATH[i]))
kegg = as_tibble(matrix(query[[1]]$GENE, nc=2, byrow=T)) %>% 
  rename(GeneID=1, descr=2) %>%
  mutate(PATH=paths$PATH[i])
  
for(i in 2:dim(paths)[1]){
  query = keggGet(paste0('mmu', paths$PATH[i]))
  if(!is.null(query[[1]]$GENE)){
    kegg = rbind(kegg, 
      as_tibble(matrix(query[[1]]$GENE, nc=2, byrow=T)) %>% 
      rename(GeneID=1, descr=2) %>%
      mutate(PATH=paths$PATH[i])
    )
  }
}

  
tmp = inner_join(
  inner_join(kegg,  
    inner_join(
      tibble(gene_stable_id=rownames(readcounts)),
      ens2entrez
    ) %>% 
    select(gene_stable_id, xref) %>% 
    unique() %>% 
    rename(ens=1, GeneID=2)
  ) %>% 
  select(ens, PATH),
  paths
)

lst = split(tmp$pathway, tmp$ens) %>% 
  lapply(unique) %>% 
  lapply(sort) %>%  
  lapply(paste, collapse='\n') 
kegg_res = tibble(ens=names(lst), KEGG=as.character(lst))  

annot = left_join(annot, kegg_res) %>% 
  rename(Ensembl=1)
  
# # https://www.biostars.org/p/69647/  
# library(annotate)
# d = getSYMBOL('118567800', data='org.Mm.eg')
# wget http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt
mgi = read_tsv('/data/sn/HErik/testosteron_220118-319026709/FASTQ_Generation_2022-01-19_01_14_16Z-515250736/MGI_Gene_Model_Coord.rpt')
# mgi %>% select(`11. Ensembl gene id`, `3. marker symbol`) %>% filter(`11. Ensembl gene id`=='ENSMUSG00000018326')

annot = left_join(annot, mgi %>% select(`11. Ensembl gene id`, `3. marker symbol`) %>% rename(Ensembl=1, symbol=2))
  
cs = createStyle(wrapText=T)
wb = createWorkbook() 
addWorksheet(wb, 'with_annotation')
writeData(wb, 1, annot)  
addStyle(wb, 1, style=cs, rows=-1, cols=-1)
saveWorkbook(wb, 'Data_annotation_GRCm39_104.xlsx', overwrite = TRUE)


   
   
