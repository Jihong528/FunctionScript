
outputpath=/Users/zhengyiyi/Desktop/projects/Ob/SmartSeq2_Ob_P7
path_normal=$outputpath/P7_control_male_filtering_exp.csv
path_treat=$outputpath/P7_ob_male_filtering_exp.csv

Run_Seurat_function.R -o $outputpath --path_normal $path_normal --path_treat $path_treat -t "Smartseq2"  --project_name "Ob_female_P7"


## Smartseq2: Ob_male_P7

outputpath=/home/zhengjh/projects/Ob/Results/Smartseq2_102620_P7_seurat
cd $outputpath
path_normal=$outputpath/P7_control_male_filtering_exp.csv
path_treat=$outputpath/P7_ob_male_filtering_exp.csv

nohup Run_Seurat_function.R -o $outputpath --path_normal $path_normal --path_treat $path_treat -t "Smartseq2" --project_name "Ob_male_P7" \
                            --folder_name "Filter_gene_500_mt_40_rp_20"  --min_cell 20  --min_gene 500  > myout.file 2>&1 &


## 10X: Ob_female_13wk_2
cd /home/zhengjh/projects/Ob/Results/20201207_ob_female_13wk_2

outputpath=/home/zhengjh/projects/Ob/Results/20201207_ob_female_13wk_2
control=/home/zhengjh/projects/Ob/Data/10X_113020_Ob_female_13wk/Ob_13wk_control-all/outs/filtered_feature_bc_matrix
obesity=/home/zhengjh/projects/Ob/Data/10X_113020_Ob_female_13wk/Ob_13wk_ob-all/outs/filtered_feature_bc_matrix

nohup Run_Seurat_function.R -o $outputpath --path_normal $control --path_treat $obesity --project_name "Ob_female_13wk_2" \
                            --folder_name "Filter_gene_500_mt_40_rp_20_3"  --min_cell 20  --min_gene 500  > myout.file 2>&1 &


