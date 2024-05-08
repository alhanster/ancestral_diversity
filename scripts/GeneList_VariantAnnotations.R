################################################################################
# Purpose: Pulling Gene Lists and Variant Annotations                          #
# Author: Alexander Han                                                        #
################################################################################

library(tidyverse)
library(data.table)
library(ggplot2)
library(R.utils)

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
dee_monoallelic <- fread("data/genelist/dee_monoallelic.tsv", header = F)
asd_monoallelic <- fread("data/genelist/asd_monoallelic.tsv", header = F)
dd_monoallelic <- fread("data/genelist/dd_monoallelic.tsv", header = F)
clingen_HI <- fread("data/genelist/clingen_HI.tsv", header = F)
mgi_essential <- fread("data/genelist/mgi_essential.tsv", header = F)
