### Plots information:

**leopard_17ind.subset.frq.png** : distribution of allele frequency per variant
    - You usually don't want extreme MAF sites because they can cause bias in the downstream analyses.
    
**leopard_17ind.subset.het.f.png** : inbreeding coefficient F per individual
**leopard_17ind.subset.het.homExp.png** : distribution of homozygosity expected per individual
**leopard_17ind.subset.het.homObs.png** : distribution of homozygosity observed per individual
    - Filtering out outliers individuals: very negative F, different from the population pattern, because it could be an antefact (ex. allelic drop out due to really low depth)
    - Filtering out outliers individuals: very positive F, different from the population pattern, because it could be an antefact (ex. contamination with other individual or sequencing errors)

**leopard_17ind.subset.idepth.png** : mean depth per individual
**leopard_17ind.subset.ldepth.mean.png** : distribution of mean depth per site (averaged across all individuals)
    - It's commom to exclude individuals and sites with very low (it can call SNPs wrong) and very high depth (it can be from paralogous and/or repetitive regions, which is very hard to be sure that it had mapped to the right place).

**leopard_17ind.subset.imiss.png** : proportion of missing data per individual
**leopard_17ind.subset.lmiss.png** : distribution of missing data proportion per site
    - The missing filter and thresholds can be more relevant for one downstream analyses than the others and can vary with your data.

**leopard_17ind.subset.lqual.png** : distribution of quality score per site (as found in the QUAL column of the VCF file)
    - High values means high confidence that it's a variant site. The best threshold is usually 30, but above 20 the chance of the site have been called wrong is very low.