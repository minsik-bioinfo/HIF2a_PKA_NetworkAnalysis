base=`pwd | sed -e "s/\/[^\/]*$//"`
echo $base
res_dir=$base"/resource/"
src_dir=$base"/src/"
data_dir=$base"/processed_rnaseq/"
result_dir=$base"/result_dir/"
htseq_matrix=$data_dir"/htseq-count-matrix-allsample"
htseq_matrix_coldKO_coldwt=$data_dir"/htseq-count-matrix-coldKO-coldwt"

<<"COMMENT"
############ DEG ############
cut -f 1-7 $htseq_matrix > $htseq_matrix_coldKO_coldwt
declare -a sample=("cold-2KO-2" "cold-2KO-4" "cold-2KO-7" "cold-wt-3" "cold-wt-5" "cold-wt-7")
Rscript $src_dir/run_DESeq.R -s "cold-2KO-2;cold-2KO-4;cold-2KO-7;cold-wt-3;cold-wt-5;cold-wt-7" -r $result_dir/deseq_result/ -o $result_dir/deseq_result/coldKO_coldwt.deseq -l "1;1;1;2;2;2" -c $htseq_matrix_coldKO_coldwt

#############################
COMMENT

############ PPI Network Mapping #####
ppi_cut=550
fc_cut=0.2
deseq=$result_dir"/deseq_result/coldKO_coldwt.deseq"
network=$res_dir"/ppi.convert_"$ppi_cut
deg=$result_dir/deseq_result/coldKO-coldwt.deseq_$fc_cut
degup=$result_dir/deseq_result/coldKO-coldwt.deseq_up_$fc_cut

awk -F"\t" -v cut=$fc_cut '{if(($3>cut||$3<cut*-1)&&($6<0.15)){print}}' $deseq > $result_dir/deseq_result/coldKO-coldwt.deseq_$fc_cut
awk -F"\t" -v cut=$fc_cut '{if(($3<cut*-1)&&($6<0.15)){print}}' $deseq > $result_dir/deseq_result/coldKO-coldwt.deseq_up_$fc_cut


awk -F"\t" 'NR==FNR{h[$1]=$1; next} {if(h[$1]!=""&&h[$2]!=""){print}}' $deg $network | awk -F"\t" '{if($2>$1){print $2"\t"$1} else{print $1"\t"$2}}' | sort | uniq > $result_dir/network_result/$fc_cut''_up_down_string_network_deseq_$ppi_cut
awk -F"\t" 'NR==FNR{h[$1]=$1; next} {if(h[$1]!=""&&h[$2]!=""){print}}' $degup $network | awk -F"\t"  '{if($2>$1){print $2"\t"$1} else{print $1"\t"$2}}' | sort | uniq > $result_dir/network_result/$fc_cut''_up_string_network_deseq_$ppi_cut


awk -F"\t" 'NR==FNR {h[$1]=$1;next} {if(h[$1]!=""&&h[$2]!=""){print}}' <(cat $deg $res_dir/Trust_tg | grep -v "Epas1") $network  | cat - $res_dir/epas_Trust_tg | awk -F"\t" '{if($2>$1){print $2"\t"$1} else{print $1"\t"$2}}' | sort | uniq > $result_dir/network_result/$fc_cut''_up_down_Trust_string_network_deseq_v2_ppi_$ppi_cut
awk -F"\t" 'NR==FNR {h[$1]=$1;next} {if(h[$1]!=""&&h[$2]!=""){print}}' <(cat $degup $res_dir/Trust_tg | grep -v "Epas1") $network  | cat - $res_dir/epas_Trust_tg | awk -F"\t" '{if($2>$1){print $2"\t"$1} else{print $1"\t"$2}}' | sort | uniq > $result_dir/network_result/$fc_cut''_up_Trust_string_network_deseq_v2_ppi_$ppi_cut


python $src_dir/network_analysis.py -input_graphs $result_dir/network_result/$fc_cut''_up_down_string_network_deseq_$ppi_cut -seed $res_dir/seed -o $result_dir/network_result/result_0.2_up_down_string_deseq_restart_0.05_$ppi_cut -e 0.05


cat $result_dir/network_result/result_0.2_up_down_string_deseq_restart_0.05_$ppi_cut | sort -r -g -k2 | grep -w -f $res_dir/MGI/filter_gene - | awk -F"\t" '{print $0"\t"NR}' | tee $result_dir/network_result/result_up_down_0.05_$ppi_cut | grep Prkaca 


awk -F "\t" 'NR==FNR{h[$1]=$3;next} {if(h[$1]!=""){print $1"\t"$2"\t"h[$1]"\t"$3}}' $deseq $result_dir/network_result/result_up_down_0.05_$ppi_cut > $result_dir/network_result/result_up_down_0.05_$ppi_cut''_fc


awk -F "\t" '{printf "%s\t%.6f\t%.6f\n",$1,$2,-$3}' $result_dir/network_result/result_up_down_0.05_$ppi_cut''_fc | cat <(echo -e "Gene\tNPscore\tlog2FC") - > $result_dir/network_result/NP_FC.tsv

python $src_dir/between_cen.py $result_dir/network_result/$fc_cut''_up_down_string_network_deseq_$ppi_cut $result_dir/network_result/ 

awk -F"\t" 'NR==FNR{h[$2] = $3; next} {if(h[$1]!=""){printf "%s\t%.6f\t%.6f\t%.6f\n",$1,$2,h[$1],-$3} else{printf "%s\t%.6f\t%.6f\t0\n",$1,$2,h[$1]}}' $result_dir/network_result/Betweenness_Centrality.tsv $result_dir/network_result/result_up_down_0.05_550_fc | sort -r -g -k3 | cat <(echo -e "Gene\tNPscore\tCentrality\tlog2FC") - > $result_dir/network_result/NP_Centrality_FC.tsv

python $src_dir/plot.py $result_dir/network_result/NP_Centrality_FC.tsv $result_dir/network_result/
