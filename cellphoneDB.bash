###### 1.Download cellphonedb
# Create python=>3.6 environment
conda create -n cpdb python=3.7
# Activate environment
source activate cpdb
# Install
pip install cellphonedb


###### 2.Use
source activate cpdb
cellphonedb method statistical_analysis  female_wt_cellphonedb_meta.txt  female_wt_cellphonedb_count.txt      --counts-data=gene_name  # 如果是直接写出基因名的加这个参数，转化为基因ID的话不用加。

# plot
#配体受体对点图 
#--rows 一列名字 选择需要的互作对 --columns 选择要看的细胞类型对
cellphonedb plot dot_plot --output-name dotplot.pdf
cellphonedb plot dot_plot --output-name dotplot_gpcr.pdf --rows rows.txt 

#细胞类型网络热图
cellphonedb plot heatmap_plot Cellphonedb_meta.txt 