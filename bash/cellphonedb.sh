conda create --name python3.6.7 python=3.6.7
conda activate python3.6.7
conda install -c conda-forge r-base r-ggplot2 r-pheatmap r-dplyr r-tidyr
pip install cellphonedb

vi /home/yah2014/anaconda3/envs/python3.6.7/lib/python3.6/site-packages/cellphonedb/src/core/database/sqlalchemy_repository/InteractionRepository.py

change
"multidata_expanded: pd.DataFrame = self.database_manager.get_repository('multidata').get_all_expanded(include_gene)"
to
"multidata_expanded = self.database_manager.get_repository('multidata').get_all_expanded(include_gene)"

#conda activate cellphonedb
cd Yang/Lung_30/Cell_Phone_DB/

for con in proximal distal terminal COPD:
    mkdir $con
    cd $con
    done
    cellphonedb method statistical_analysis ../Lung_30-${con}_meta.data_cellphone.txt ../Lung_30-${con}_counts.txt
    cellphonedb plot dot_plot
    cellphonedb plot heatmap_plot ../Lung_30-${con}_meta.data_cellphone.txt
