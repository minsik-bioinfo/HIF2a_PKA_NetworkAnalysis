# Adipocyte HIF2α functions as a thermostat via PKA Cα regulation in beige adipocytes [[link]](https://www.nature.com/articles/s41467-022-30925-0)
<p align="center">
  <img src="sample_plot.png"/>
</p>


## About data

The RNA-seq data ("raw" and "processed") of this study has been deposited in the Gene Expression Omnibus (GEO) database, [[GSE179385]](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179385)

## Usage
```
bash src/start.sh
```
---

##### *About path*
*The path is detected automatically in the script code, but if there is a problem, you have to manually modify the "base" variable to your own path.*

---

## Result files

result_dir/deseq_result/coldKO_coldwt.deseq : DEG analysis result between coldKO and coldwt samples.
result_dir/deseq_result/coldKO_coldwt.deseq_0.2 : DEG analysis result between coldKO and coldwt samples with the significance thresholds were |log2FC | > 0.2 and P < 0.15.
result_dir/deseq_result/coldKO_coldwt.deseq_0.2_up : Up-regulated gene list in coldKO samples with above thresholds.

result_dir/network_result/

## Contact
oominsik@gmail.com
