#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

bx<-commandArgs(trailingOnly=TRUE)[2]
message(sprintf("Processing %s",bx))

require(data.table)
require(seqinr)

tabi<-fread(bx,header=FALSE,sep="\t")
colnames(tabi)<-c("qseqid", "qlen", "sseqid", "slen", "pident", "qstart", "qend", "sstart", "send", "length", "mismatch", "gaps", "evalue", "bitscore")
tabi$SampleID<-gsub(".RE8.RE13.bla","",basename(bx))

tabi$MaxScore<-ave(tabi$bitscore,tabi$qseqid,FUN=max)
tabi$Orientation<-"FWD"
tabi$Orientation[tabi$sstart>tabi$send]<-"REV"

bREl<-vector("list",length(unique(tabi$sseqid)))
names(bREl)<-unique(tabi$sseqid)

fx<-commandArgs(trailingOnly=TRUE)[3]
fai<-read.fasta(file=fx,seqtype="DNA",as.string=TRUE, forceDNAtolower = FALSE)
fai.rc<-fai
z<-lapply(fai[names(fai) %in% tabi$qseqid[tabi$Orientation %in% "REV"]],function(X)toupper(paste0(rev(comp(s2c(X))),collapse="")))
fai.rc[names(fai) %in% tabi$qseqid[tabi$Orientation %in% "REV"]]<-z

for(i in seq_along(bREl)){
    bREl[[i]]<-tabi[tabi$sseqid %in% unique(tabi$sseqid)[i]]

    faREx<-fai.rc[names(fai.rc) %in% bREl[[i]]$qseqid] 
    write.fasta(sequences=faREx,names=names(faREx),file.out=paste0(gsub(".RE8.RE13.bla","",basename(bx)),".",unique(tabi$sseqid)[i],".hits.fasta"),as.string=TRUE)


    print(paste0(i,"_processed"))

}


