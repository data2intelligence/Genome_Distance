# Genome_Distance

Two steps to identify distance between two genes:  
1, Identify the nearest probe with a confidence score = exp(-a*x). If we assign a half-confidence distance as D, a = ln(2)/D  
2, convert from pixel distance to physical distance as sqrt((103nM * (x1-x2))^2 + (103nM * (y1-y2))^2 + (250nM * (z1-z2))^2), then convert nM to uM.


Two steps to insert new data:  
1, Format data into the same format as Long Cai paper's supplementary materials  
2, run the following command lines:  
python manage.py build_data insert_geneprobe your_gene_probe_file  
python manage.py build_data insert_loci your_loci_file  
