# installation of libsbml
Ref: https://synonym.caltech.edu/software/libsbml/libsbml-docs/installation/
pip install python-libsbml

/opt/R/4.1.2/bin/R CMD INSTALL /home/jihun/Downloads/libSBML_5.19.0.tar.gz

# Installation of rsbml
## Method 1
BiocManager::install("rsbml") # it may not work.
## MEthod 2
Download 


# INstallation of piano
BiocManager::install("piano", dependencies=TRUE)
library(piano)
?piano
?runGSA


GEM_liver_hepatocytes = loadGSC(liver_hepatocytes.xml)
load('/R/default_working_directory/Resources/reporter metabolite assay/hmrGsc.RData')
mmrGsc대신 hmrGsc

#HCC(male tumor w/o Hepatitis vs female tumor w/o Hepatitis)
gsa.HCC_noHepatitis = runGSA(input_reporterMetabolite_padj, input_reporterMetabolite_lfc, gsc=hmrGsc, geneSetStat="reporter", signifMethod="nullDist", nPerm=1000, gsSizeLim=c(5,100))
gsa.HCC_noHepatitis.summary = GSAsummaryTable(gsa.HCC_noHepatitis)
##network.HCC_noHepatitis <- networkPlot(gsa.HCC_noHepatitis, class="distinct", direction="both", significance=0.005, label="numbers")

save(gsa.HCC_noHepatitis, file='gsa.HCCnoH__2021.06.12.rda')
save(gsa.HCC_noHepatitis.summary, file= 'gsa.gsa.HCCnoH_2021.06.12.rda')
write.table(gsa.HCC_noHepatitis.summary, file = "gsa.HCCnoH_summary_2021.06.12.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote=FALSE, append=FALSE)
