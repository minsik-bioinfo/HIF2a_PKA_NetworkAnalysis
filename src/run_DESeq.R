# load DESeq library
library(DESeq2)
library(optparse)

# make options
option_list <- list( 
    make_option(c("-s", "--sample_list"), action="store"),
    make_option(c("-r", "--result_dir"), action="store"),
    make_option(c("-o", "--deseq_result"), action="store"),
    make_option(c("-l", "--sample_num_list"), action="store"),
    make_option(c("-c", "--count_table"), action="store")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# parsing arguments
result_dir <- opt$result_dir
sample_list <- unlist(strsplit(opt$sample_list, ";"))
sample_num_list <- unlist(strsplit(opt$sample_num_list, ";"))
count_table_file <- opt$count_table
final_result_file <- opt$deseq_result

output_prefix <- paste(sample_list, collapse="_")


print("**** arguments ****")
result_dir
sample_list
count_table_file
final_result_file
sample_num_list
output_prefix
print("**** arguments ****")

# read couting table
mycounts <- read.table(count_table_file, header=TRUE, row.names=1)
#print (colnames(mycounts),rownames(mycounts))
# set sample information
# temp <- c(test[3:length(test)]) handle input argument like this to make below
		# for DESeq condition list	

#samples <- data.frame(row.names=c("Nip_2D_0_1_head100000.fastq","Nip_2D_6_1_head100000.fastq"), condition=as.factor(c(rep("Nip_2D_0_1_head100000.fastq",1),rep("Nip_2D_6_1_head100000.fastq"))))
samples <- data.frame(row.names=sample_list, condition=as.factor(sample_num_list))
# create DESeq dataset
colnames(mycounts) <-NULL
dds <- DESeqDataSetFromMatrix(countData = mycounts, colData=samples, design=~condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts<-counts(dds,normalized=TRUE)
#write.table(normalized_counts,file=paste(result_dir,"/normalized.txt"),sep="\t",quote=F,col.names=NA)
#exit(0)
# run DESeq analysis
dds <- DESeq(dds)

# store result
#res <- results(dds, contrasts=c("conditions", sample_num_list))
res <- results(dds)

#res <- results(dds)
res <- res[order(res$padj),]


res

write.table(res,file=final_result_file, quote= FALSE , sep="\t")

