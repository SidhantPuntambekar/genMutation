library(dplyr)




fail.set <- dplyr::filter(syn.vars, freq > 0, ac==0 | filter!="PASS" | has_star | lcr | rf_label=="FP" | allele_type != "snv")
fail.labels <- paste0(fail.set$chrom, "_", fail.set$pos)
syn.vars$label <- paste0(syn.vars$chrom, "_", syn.vars$pos)
syn.vars <- dplyr::filter(syn.vars, ! label %in% fail.labels)
syn.vars <- dplyr::filter(syn.vars, !(freq==0 & interp_dist==0 & allele_type != "snv"))
syn.vars <- dplyr::filter(syn.vars, qual=="high", over_20>0.8,
                          !( freq > 0 & ( ac==0 | filter!="PASS" | has_star | lcr | rf_label=="FP" | allele_type != "snv")),
                          gene != "ENSG00000242265")