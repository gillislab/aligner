# Code and scripts for: 
## The fractured landscape of RNA-seq alignment: The default in our STARs

Many tools are available for RNA-seq alignment and expression quantification, with comparative value being hard to establish. Benchmarking assessments often show high performance, with occasional outliers, but are often focused on either model data or fail to explain variation in performance in detail. This leaves us to ask, what is the most meaningful way to assess different alignment choices? And importantly, where is there room for progress? In this work, we explore the answers to these two questions by performing an exhaustive assessment of the STAR aligner. We assess STAR’s performance across a range of alignment parameters using common metrics, and then on biologically focused tasks. We find technical metrics such as fraction mapping and correlations uninformative, capturing properties unlikely to have any role in biological discovery. Surprisingly, we find that changes in alignment parameters within a wide range have very little impact on likely biological and not just technical performance. Yet, when performance finally does break, it happens in difficult regions, such as X Y paralogs and MHC genes. We believe improved reporting by developers will help establish where results are likely to be robust or fragile, providing a better baseline to establish where methodological progress can still occur. 

# The aligner: 
## STAR 
Available here: 
https://github.com/alexdobin/STAR

# The assesor: 
## AuPairWise
Available here: 
https://github.com/gillislab/aupairwise


![summary](https://github.com/sarbal/aligner/blob/master/imgs/icefloe.png "summary")
