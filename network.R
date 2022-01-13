library(STRINGdb)
library(biomaRt)
library(igraph)
library(RJSONIO)

#databases
ciona <- STRINGdb$new(version='11',species=7719,score_threshold=200, input_directory="")
mouse <- STRINGdb$new(version='11',species=10090,score_threshold=200, input_directory="")

mmProt <- mouse$get_proteins()
mmProt$protein_external_id <- sub('10090\\.','',mmProt$protein_external_id)

#genes of interest
genes <- c('Arhgef8','Gna12.13','Farp','Depdc','Tyrosinase','Eph','Foxf','Gata','Rho',
	   'Dync2h','NDE','EPB4.1','GNA.L.S','Rab5','Paxillin','Rabep','Rock','Col9a1',
	   'Ddr' ,'Daam','ADF.Cofilin')

khid <- c('KH.C8.840','KH.L37.10','KH.C8.441','KH.L108.56','KH.C12.469','KH.C1.404',
	  'KH.C3.170','KH.L20.1','KH.C2.651','KH.C14.432','KH.C3.451','KH.C9.562',
	  'KH.L121.2','KH.L22.59','KH.L108.1','KH.C11.355','KH.C7.324','KH.C8.248',
	  'KH.C9.371','KH.C10.209','KH.L37.33')

khid <- paste0("KH2012:",khid)

khname <- do.call(rbind,fromJSON('all_genes.json'))
row.names(khname) <- khname[,1]

genes <- khname[khid,]

genes <- do.call(rbind,apply(
	genes,1,
	function(x) cbind(x[1],unlist(strsplit(x[2],'; ')))
))

#Gene names differ from KH2013 genome
genes <- data.frame(khid=paste0("KH2013:",khid),name=genes)
gene.names <- read.delim("gene_name.txt", row.names=1, stringsAsFactors=F)
genes <- merge(genes,gene.names,by.x=1,by.y=0)

#generate regex from gene names to search mouse db
regex <- paste0('^',gsub('\\(\\|','\\(',gsub(
  '([A-Z])|','\\1',
  gsub(
    '/','|',
    gsub(
      "([0-9/]+)","\\(\\1\\)a?",
      gsub('([0-9])([0-9][A-Za-z])','\\1\\.?\\2',genes[,2]#genes$UniqueNAME
    ))))),'$')

#select proteins matching regex
prots <- lapply(
	regex,#setNames(regex,genes$UniqueNAME),
	function(x) mmProt[grep(x,mmProt$preferred_name,T),]
)
sapply(prots,dim)

#merge with KHIDs
prots <- do.call(rbind,mapply(
      function(x,y) {
	      if(nrow(y)>0){
		      y$UniqueNAME <- x
		      merge(genes,y)
	      }
      },names(prots),prots
))

#ENSEMBL IDs for KH2013 genome
ensembl <- read.delim('KH-ENS.blast',stringsAsFactors = F,header=F)

#get ENSEMBL geneIDs from transcript IDs
mart <- useMart('ensembl','cintestinalis_gene_ensembl')
transcriptToGene <- select(mart,ensembl[,1],c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),"ensembl_transcript_id")
ensembl.gene <- data.frame(ensembl_transcript_id=ensembl[,1],GeneID=sub('KH2012:(.*).v.*','KH2013:\\1',ensembl[,2]))
ensembl.gene <- merge(ensembl.gene,transcriptToGene,'ensembl_transcript_id')
ensembl.gene <- ensembl.gene[!duplicated(ensembl.gene[,-1]),-1]

genes <- merge(genes,ensembl.gene,by=1,all.x=T)
sel <- sapply(khid,grep,ensembl[,2])

#look up orthologs for genes that don't share a name
query <- setdiff(genes$ensembl_gene_id,prots$ensembl_gene_id)

genes.mm <- select(mart,query,c("ensembl_gene_id","mmusculus_homolog_ensembl_peptide"),'ensembl_gene_id')
genes.mm <- merge(genes,genes.mm,'ensembl_gene_id')

genes.mmProt <- merge(mmProt,genes.mm,by.y='mmusculus_homolog_ensembl_peptide',by.x='protein_external_id')

prots <- genes.mmProt[,c(6,4,1:3,5,7:9)]
names(prots) <- names(genes.mmProt)
prots <- rbind(prots,genes.mmProt)
prots <- prots[!duplicated(prots[,c(1,3:9)]),]

pdf('mmusculus_homologs.pdf')
mouse$plot_network(prots$protein_external_id)
dev.off()
