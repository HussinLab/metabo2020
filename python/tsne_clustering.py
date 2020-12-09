
import pandas as pd
import numpy as np
import pickle
import os
from bokeh.plotting import figure
from bokeh.io import show, output_notebook, push_notebook
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import LinearColorMapper, CategoricalColorMapper, CDSView, GroupFilter, Button, CustomJS
from bokeh.palettes import Category20
from bokeh.transform import factor_cmap
from bokeh.layouts import gridplot
from bokeh.plotting import figure, output_file, save

from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

from definitions import *

import argparse


def compute_PCA(X):
    pca = PCA(n_components=2)
    pca.fit(X)
    print(pca.explained_variance_ratio_)
    
    p = pca.transform(X) 
    
    xx = [a for (a, b) in p]
    yy = [b for (a, b) in p]
    
    return xx, yy

def compute_tsne(X, ppl=30, n_components=2, n_iter=1000, lr=200, random_state=None):
    X_embedded = TSNE(n_components=n_components, perplexity=ppl, n_iter=n_iter, random_state=random_state,
                      learning_rate=lr).fit_transform(X)
    xx = [a for (a, b) in X_embedded]
    yy = [b for (a, b) in X_embedded]
    
    return xx, yy

def recompute_clusters(data, eps, min_samples, use_clusters_of=None):
    """
    in case you want to update the clustering parameters
    """
    result = data
    
    if use_clusters_of is None:
        print("No clusters passed, each graph will be clustered independently")
        for i, split in enumerate(result):
            tsne_xs = split['tsne_xs']
            tsne_ys = split['tsne_ys']

            clustering = DBSCAN(eps=eps,
                            min_samples=min_samples).fit(np.array([[x, y] for x, y in zip(tsne_xs, tsne_ys)]))

            clusters = clustering.labels_.astype(str)
            result[i]['clusters'] = clusters
    else:
        print(f"using the clustering of {use_clusters_of}")
        found= False
        for i, category in enumerate(result):
            if category['name'] == use_clusters_of:
                found = True
                break
        if not found:
            raise NotFoundException
            
        tsne_xs = category['tsne_xs']
        tsne_ys = category['tsne_ys']

        clustering = DBSCAN(eps=eps,
                        min_samples=min_samples).fit(np.array([[x, y] for x, y in zip(tsne_xs, tsne_ys)]))

        clusters = clustering.labels_.astype(str)
        
        for i in range(len(result)):
            result[i]['clusters'] = clusters
        
    return result


def parse_args():
    """
    Parse the arguments from the command line
    All arguments have default values
    """
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--data_path',  # Specified in the config file
                        help='Folder containing all the datasets',
                        default='../data',
                        type=str)
    
    parser.add_argument('--f1name',  # Specified in the config file
                        help='Name of the normalized dataset file',
                        default='data_normalized_all_BL_DC.csv',
                        type=str)
    
    parser.add_argument('--f2name',  # Specified in the config file
                        help='Name of the metadata file',
                        default='patient_IDs_metadata.csv',
                        type=str)
    
    parser.add_argument('--load_precomputed',  # Specified in the config file
                        help='If the file has already been computed, just reuse it',
                        default=True,
                        type=bool)
    
    parser.add_argument('--random_seed',  # Specified in the config file
                        help='random seed. It can help with reproducibility',
                        default=888,
                        type=int)
    
    parser.add_argument('--tsne_ppl',  # Specified in the config file
                        help='random seed. It can help with reproducibility',
                        default=7,
                        type=int)
    
    parser.add_argument('--tsne_lr',  # Specified in the config file
                        help='random seed. It can help with reproducibility',
                        default=200,
                        type=int)
    
    parser.add_argument('--dbscan_eps',  # Specified in the config file
                        help='random seed. It can help with reproducibility',
                        default=8,
                        type=int)
    
    parser.add_argument('--dbscan_min_samples',  # Specified in the config file
                        help='dbscan_ seed. It can help with reproducibility',
                        default=4,
                        type=int)
    return parser
if __name__ == "__main__":
    parser = parse_args()
    parser = parser.parse_args()
    
    data_path = parser.data_path #'../Data/original'
    file_1 = parser.f1name#"data_normalized_all_BL_DC.csv"
    file_2 = parser.f2name#"patient_IDs_metadata.csv"
    
    compare_all_relationship_method = "kendall"
    
    compare_all_tsne_ppl = parser.tsne_ppl
    compare_all_tsne_n_iter = 5000
    compare_all_tsne_lr = parser.tsne_lr #200
    compare_all_tsne_random_state = parser.random_seed
        
    compare_all_dbscan_eps = parser.dbscan_eps #2
    compare_all_dbscan_min_samples = parser.dbscan_min_samples # 7
    
    #First time run, set it to false, it will take a bit of time
    load_precomputed = parser.load_precomputed#True
    
    use_special_metab_shapes = True
    use_clusters_of = 'BL Only'
    
    first_analysis_df = pd.read_csv(os.path.join(data_path,file_1), sep='\t')
    meta_1st_analysis_df = pd.read_csv(os.path.join(data_path,file_2), sep='\t')

    #Rename the 1st column to be samples
    first_analysis_df['samples'] = first_analysis_df['metabolites']
    meta_1st_analysis_df['CONDITION'] = meta_1st_analysis_df['TIME']
    
    filtered_df = first_analysis_df.copy()
    filtered_df = pd.merge(filtered_df, meta_1st_analysis_df[['SAMPLE_NAME','AGE', 'SEX', 'CONDITION', 'TREATMENT']], left_on='samples', right_on='SAMPLE_NAME')
    filtered_df.drop(['metabolites'], axis=1, inplace=True)
    
    bl_true_condition = filtered_df['CONDITION'] == 'BL'
    dc_true_condition = ~bl_true_condition
    
    female_condition = filtered_df['SEX'] == 'Female'
    male_condition = ~female_condition
    
    young_age_condition = filtered_df['AGE'] < young_age_thresh
    old_age_condition = ~young_age_condition
    
    clo_condition = filtered_df['TREATMENT'] == 'Clo'
    tig_condition = ~clo_condition
    
    filtered_df.drop(['samples', 'SAMPLE_NAME', 'AGE', 'SEX', 'CONDITION'], axis=1, inplace=True)
    
    df_bl = filtered_df[bl_true_condition]
    df_dc = filtered_df[dc_true_condition]
    
    df_female = filtered_df[female_condition]
    df_male =   filtered_df[male_condition]
    
    df_young = filtered_df[young_age_condition]
    df_old =   filtered_df[old_age_condition]
    
    df_bl_female = filtered_df[bl_true_condition & female_condition]
    df_bl_male =   filtered_df[bl_true_condition & male_condition]
    
    df_dc_female = filtered_df[dc_true_condition & female_condition ]
    df_dc_male =   filtered_df[dc_true_condition & male_condition]
    
    df_bl_young_female = filtered_df[bl_true_condition & young_age_condition & female_condition]
    df_bl_young_male =   filtered_df[bl_true_condition & young_age_condition & male_condition]
    
    df_dc_young_female = filtered_df[dc_true_condition & young_age_condition & female_condition]
    df_dc_young_male =   filtered_df[dc_true_condition & young_age_condition & male_condition]
    
    df_bl_old_female = filtered_df[bl_true_condition & old_age_condition & female_condition ]
    df_bl_old_male =   filtered_df[bl_true_condition & old_age_condition & male_condition ]
    
    df_dc_old_female = filtered_df[dc_true_condition & old_age_condition & female_condition]
    df_dc_old_male =   filtered_df[dc_true_condition & old_age_condition & male_condition]
    
    df_bl_young = filtered_df[bl_true_condition & young_age_condition ]
    df_bl_old =   filtered_df[bl_true_condition &  old_age_condition]
    
    df_dc_young = filtered_df[dc_true_condition & young_age_condition]
    df_dc_old =   filtered_df[dc_true_condition & old_age_condition]
    
    #Code used to inser dc and clo
    df_bl_clo = filtered_df[bl_true_condition & clo_condition ]
    df_dc_clo = filtered_df[dc_true_condition & clo_condition ]
    
    df_bl_tig = filtered_df[bl_true_condition & tig_condition ]
    df_dc_tig = filtered_df[dc_true_condition & tig_condition ]
    
    df_splits = [{"DC and Clo": df_dc_clo}, {"DC and Tig": df_dc_tig}, {"BL and Clo": df_bl_tig}, {"BL and Tig": df_bl_tig}]
    
    df_splits = [{"BL Only": df_bl}, {"DC Only": df_dc}, 
                {"BL and Clo": df_bl_clo}, {"DC and Clo": df_dc_clo}, 
                {"BL and Tig": df_bl_tig}, {"DC and Tig": df_dc_tig},
                {"Female Only": df_female}, {"Male Only": df_male}, 
                {"Young Only": df_young}, {"Old Only": df_old}, 
                {"BL and Female": df_bl_female}, {"BL and Male": df_bl_male}, 
                {"DC and Female": df_dc_female}, {"DC and Male": df_dc_male}, 
                {"BL and Young and Female": df_bl_young_female}, 
                {"BL and Young and Male": df_bl_young_male}, 
                {"DC and Young and Female": df_dc_young_female}, 
                {"DC and Young and Male": df_dc_young_male}, 
                {"BL and Old and Female": df_bl_old_female}, 
                {"BL and Old and Male": df_bl_old_male}, {"DC and Old and Female": df_bl_old_female}, 
                {"DC and Old and Male": df_bl_old_male}, {"BL and Young": df_bl_young}, {"BL and Old": df_bl_old}, 
                {"DC and Young": df_dc_young}, {"DC and Old": df_dc_old}]
    
    if load_precomputed:
        with open('../data/compare_all_data_output_allbls.pickle', 'rb') as f:
            compare_all_data = pickle.load(f)
            
    else:
        compare_all_data = []

        for i, split in enumerate(df_splits):
            k,v = list(*split.items())
            current_data_dict = {}

            corrs = v.corr(compare_all_relationship_method)
            tsne_xs, tsne_ys = compute_tsne(corrs, ppl= compare_all_tsne_ppl, 
                                            n_iter=compare_all_tsne_n_iter, 
                                            random_state=compare_all_tsne_random_state,
                                            lr=compare_all_tsne_lr)
            clustering = DBSCAN(eps=compare_all_dbscan_eps,
                            min_samples=compare_all_dbscan_min_samples).fit(np.array([[x, y] for x, y in zip(tsne_xs, tsne_ys)]))

            clusters = clustering.labels_.astype(str)

            current_data_dict['name'] = k
            current_data_dict['corrs'] = corrs
            current_data_dict['tsne_xs'] = tsne_xs
            current_data_dict['tsne_ys'] = tsne_ys
            current_data_dict['clusters'] = clusters

            compare_all_data.append(current_data_dict)

            print(f'{i+1}- {k} computation is done.')
            
        with open('../data/compare_all_data_output_allbls.pickle', 'wb') as f:
            pickle.dump(compare_all_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    compare_all_data = recompute_clusters(compare_all_data, eps=7.1, min_samples=4, use_clusters_of=use_clusters_of)
    
    #Remove redundunt entries from list
    compare_all_data_unfiltered = compare_all_data
    compare_all_data = []
    seen_before = []
    for i, e in enumerate(compare_all_data_unfiltered):
        d = e['name']
        if d not in seen_before:
            compare_all_data.append(e)
            seen_before.append(d)
        else: 
            continue
        
    data_dict = {}
    ordered_names = []
    data_dict['metabol'] = list(compare_all_data[0]['corrs'].columns)
    for i, split in enumerate(compare_all_data):
        ordered_names.append(split['name'])
        data_dict['name_'+str(i)] = [split['name']]*len(split['tsne_xs'])
        data_dict['x_'+str(i)] = split['tsne_xs']
        data_dict['y_'+str(i)] = split['tsne_ys']
        data_dict['cluster_'+str(i)] = list(split['clusters'].astype(str))
        
        shape_list = list(['circle']*len(split['tsne_xs']))
        if use_special_metab_shapes:
            for shape, metabol_to_shape in shapes_dict.items():
                for j, metabol in enumerate(data_dict['metabol']):
                    if metabol in metabol_to_shape:
                        shape_list[j] = shape
        data_dict['shape_'+str(i)] = shape_list
        
        
    with open('data_dict_independent_clustering_allbls.pickle', 'wb') as f:
        pickle.dump(data_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    compare_all_TOOLS = "box_select,lasso_select,help,save"
    compare_all_gridplot = []
    compare_all_single_plot_dims = (350,350)

    source = ColumnDataSource(data=pd.DataFrame(data_dict))
    
    factors = np.array(list(range(-1, 19))).astype(str)

    for i, split in enumerate(compare_all_data):
        circles = CDSView(source=source, filters=[GroupFilter(column_name=f'shape_{i}', group='circle')])
        triangles = CDSView(source=source, filters=[GroupFilter(column_name=f'shape_{i}', group='triangle')])
        squares = CDSView(source=source, filters=[GroupFilter(column_name=f'shape_{i}', group='square')])
        diamonds = CDSView(source=source, filters=[GroupFilter(column_name=f'shape_{i}', group='diamond')])
        star = CDSView(source=source, filters=[GroupFilter(column_name=f'shape_{i}', group='*')]  )

        if i%2 == 0:
            current_plot_line = []
            compare_all_gridplot.append(current_plot_line)
        
        TOOLTIPS = [
                ("metabo_name", "@metabol"),
                ("cluster", f"@cluster_{i}"),
                ('name', f'@name_{i}')
            ]
        
        fig=figure(title=split['name'],
                plot_width=compare_all_single_plot_dims[0], 
                plot_height=compare_all_single_plot_dims[1], 
                tooltips=TOOLTIPS, 
                tools=compare_all_TOOLS)
        
        c = fig.circle(f'x_{i}', f'y_{i}', size=5, source=source, view=circles, alpha=0.6,
                    color=factor_cmap(f'cluster_{i}', palette=Category20[20], factors=factors),
                    selection_color=factor_cmap(f'cluster_{i}', palette=Category20[20], factors=factors),
                    nonselection_alpha=0.0 )
        
        fig.triangle(f'x_{i}', f'y_{i}', size=7, source=source, view=triangles,alpha=0.6, line_color='black', 
                    color=factor_cmap(f'cluster_{i}', palette=Category20[20], factors=factors),
                    selection_color=factor_cmap(f'cluster_{i}', palette=Category20[20], factors=factors),
                    nonselection_alpha=0.01 )
        
        fig.square(f'x_{i}', f'y_{i}', size=7, source=source, view=squares,alpha=0.6, line_color='black',
                    color=factor_cmap(f'cluster_{i}', palette=Category20[20], factors=factors),
                    selection_color=factor_cmap(f'cluster_{i}', palette=Category20[20], factors=factors),
                    nonselection_alpha=0.01 )
        
        fig.diamond(f'x_{i}', f'y_{i}', size=7, source=source, view=diamonds,alpha=0.6, line_color='black',
                    color=factor_cmap(f'cluster_{i}', palette=Category20[20], factors=factors),
                    selection_color=factor_cmap(f'cluster_{i}', palette=Category20[20], factors=factors),
                    nonselection_alpha=0.01 )
                        
        fig.circle_cross(f'x_{i}', f'y_{i}', size=7, source=source, view=star,alpha=0.6, line_color='black',
                    color=factor_cmap(f'cluster_{i}', palette=Category20[20], factors=factors),
                    selection_color=factor_cmap(f'cluster_{i}', palette=Category20[20], factors=factors),
                    nonselection_alpha=0.01 )
        
        current_plot_line.append(fig)

    comparisons_plot = gridplot(compare_all_gridplot)
    selected_mets = []
    
    source.selected.js_on_change('indices', CustomJS(args=dict(s1=source), 
                                                     code="""
                                                            var inds = cb_obj.indices;
                                                            console.log(inds)
                                                            var selected_mets = "["
                                                            for (var i = 0; i < inds.length; i++) 
                                                            {
                                                                console.log(s1.data.metabol[inds[i]])
                                                                selected_mets = selected_mets + "'" + s1.data.metabol[inds[i]] + "',"
                                                            }
                                                            selected_mets = selected_mets + "'0']"
                                                            console.log("selected_mets="+selected_mets);
                                                        """))
    
    output_file("tsne.html")
    save(comparisons_plot)
    
