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
awk -F"\t" 'NR==FNR{h[$1]=$1; next} {if(h[$1]!=""&&h[$2]!=""){print}}' $degup $net | awk -F"\t"  '{if($2>$1){print $2"\t"$1} else{print $1"\t"$2}}' | sort | uniq > $result_dir/network_result/$fc_cut''_up_string_network_deseq_$ppi_cut


