################################################################################
# Purpose: Pulling Gene Lists and Variant Annotations                          #
# Author: Alexander Han                                                        #
################################################################################

# Variant Annotations
SYN <- c("synonymous_variant")
MISSENSE <- c("missense_variant", "inframe_deletion", "inframe_insertion")
LOF <- c("frameshift_variant", "stop_gained", "start_lost", 
         "splice_acceptor_variant", "splice_donor_variant",
         "stop_gained&frameshift_variant", "splice_donor_variant&coding_sequence_variant&intron_variant", 
         "splice_donor_variant&intron_variant", "stop_gained&inframe_insertion", "frameshift_variant&stop_lost", 
         "frameshift_variant&start_lost", "start_lost&splice_region_variant", 
         "stop_gained&protein_altering_variant", "stop_gained&frameshift_variant&splice_region_variant", 
         "stop_gained&inframe_deletion", "frameshift_variant&stop_retained_variant")

# Gene Lists
setwd(paste0(
    "https://github.com/alhanster/ancestral_diversity/tree/b7655122ec234847c11dba025f2563ac69057292/data/genelist"))
dee_monoallelic <- fread("dee_monoallelic.tsv", header = F)
asd_monoallelic <- fread("asd_monoallelic.tsv", header = F)
dd_monoallelic <- fread("dd_monoallelic.tsv", header = F)
clingen_HI <- fread("clingen_HI.tsv", header = F)
mgi_essential <- fread("mgi_essential.tsv", header = F)