import scanpy as sc
import scvelo as scv
import anndata as ad
import numpy as np
import scanpy.external as sce
import matplotlib.pyplot as plt
import pertpy as pt
import seaborn as sns
import pandas as pd

def preprocess_adata(adata, n_top_genes=2000):
    # QC
    adata.var["mt"] = adata.var_names.str.startswith("mt:")
    adata.var["ribo"] = adata.var_names.str.startswith("RpS", "RpL")
    # calculate the metrics 
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace = True, percent_top = None, log1p=False)
    # Filter (keep high count ones (we will remove doublets later)
    sc.pp.filter_cells(adata, min_genes=250) # filter cells that have less than 250 genes
    sc.pp.filter_cells(adata, min_counts=500) # filter cells that have less than 500 counts
    sc.pp.filter_genes(adata, min_cells=5) # filter genes that are expressed in less than 5 cells
    adata = adata[(adata.obs.pct_counts_mt>0.5)&(adata.obs.pct_counts_mt<30)].copy()
    # keep filtered counts in new layer
    adata.layers['counts'] = adata.X.copy()
    # normalize
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    """# ? why does gene filter after normalization? """
    # ? sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    return adata

def merge_a_v_datas(adata, vdata):
    # make both datasets unique
    """make cell barcode names the same in adata_velocyto as in adata (remove the first 7 characters from 
    the cell barcodes in the velocyto data 
    and then add the sample_name to the end of the cell barcodes with "-" in between)"""
    ## rename in the obs sample_id to sample_name in the velocyto data
    vdata.obs.rename(columns={"sample_id": "sample_name"}, inplace=True)
    ## remove sample names in the beginning of the cell barcode names 
    vdata.obs.index = vdata.obs.index.map(lambda x: x[7:])
    ##  convert the column sample_name to object because I wont be able to add sample name to the end of the cell barcode if it is a category
    vdata.obs["sample_name"] = vdata.obs["sample_name"].astype("object") # convert the column to object
    ## add sample name to the end of the cell barcode name
    vdata.obs.index = vdata.obs.index + "-" + vdata.obs["sample_name"] # add sample name to the end of the cell barcode name
    ## rename the index column in the velocyto data
    vdata.obs.index.name = "cell_barcode"
    # check the cells which overlap
    shared_cells = adata.obs_names.intersection(vdata.obs_names)
    shared_cells
    
    # select cells which are overlaping in both datasets
    ## assign the cells from adata which are in shared cells
    adata_shared = adata[adata.obs_names.isin(shared_cells)] 
    print(adata_shared)
    # select cells which intersect for adata_velocyto
    vdata_shared = vdata[vdata.obs_names.isin(shared_cells)]
    print(vdata_shared)
    print(f" the number of cells which in the same order are {sum(adata_shared.obs_names == vdata_shared.obs_names)}")

    """The cells are not in the same order in both datasets as well as genes
    and the number of genes differ"""
    # sort cell barcodes in the same order in both datasets

    ## sort the cell barcodes and genes in the adata
    # assign to list common genes and cells
    common_genes = list(set(adata.var.index) & set(vdata.var.index))
    common_cells = list(set(adata.obs.index) & set(vdata.obs.index))

    ## align the data
    aligned_adata = adata[common_cells, common_genes] # .copy is important!!
    aligned_vdata = vdata[common_cells, common_genes]

    # check if the cells are in the same order
    if sum(aligned_adata.obs_names == aligned_vdata.obs_names) == aligned_adata.shape[0]: # if the number of cells in the same order equal to the number of cell in the anndata objects
        print("The cell barcodes are in the same order in both datasets")

    # check if the genes are in the same order
    print(f"The number of genes which are in the same order: {sum(aligned_adata.var_names == aligned_vdata.var_names)}")

    # transfer layers from aligned_vdata to aligned_adata object
    aligned_adata.layers["spliced"] = aligned_vdata.layers["spliced"].copy()
    aligned_adata.layers["unspliced"] = aligned_vdata.layers["unspliced"].copy()    
    aligned_adata.layers["ambiguous"] = aligned_vdata.layers["ambiguous"].copy()

    # change the column type to category for plotting
    aligned_adata.obs["sample_name"] = aligned_adata.obs["sample_name"].astype("category")
    aligned_vdata.obs["sample_name"] = aligned_vdata.obs["sample_name"].astype("category")
    vdata.obs["sample_name"] = vdata.obs["sample_name"].astype("category")

    # save the results to h5ad 
    """The data is merged correctly"""
    aligned_adata.write("/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/data/a_v_datas_merged_DECODE_v1_local_machine.h5ad")


# function for RNA velocity estimation
def velocity_inferring(adata, min_shared_counts=None, min_shared_cells=None, n_top_genes=None, mode = "stochastic", condition=None,
                       annotation = None, moments_approach=None, groupby=None, groups=None, groups_for_fit=None, perc=[5,95], min_r2=0.01):
    if min_shared_counts is None:
        scv.pp.filter_genes(adata, min_shared_cells=min_shared_cells) 
    else:
        print("\n Applying min_shared_counts filter \n")
        scv.pp.filter_genes(adata, min_shared_counts=min_shared_counts)
        
    scv.pp.normalize_per_cell(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    sc.pp.log1p(adata)

    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    if "X_diffmap" in adata.obsm:
        del adata.obsm["X_diffmap"]
    
    print("'\n Computing PCA, knn and UMAP \n ")
    # compute the PCA
    scv.pp.pca(adata, n_comps=50) # 50 is the default value
    # compute the knn
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15) # knn 15 is the default value
    # compute the UMAP
    scv.tl.umap(adata, n_components=2, min_dist=0.8) # 2 is the default value
    print("\n Computing RNA velocity \n ")
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    if mode == "stochastic":
        print("\n Applying stochastic mode \n")
        scv.tl.velocity(adata, mode=mode, groupby=groupby, groups_for_fit=groups_for_fit, perc=perc, min_r2=min_r2)
        scv.tl.velocity_graph(adata, n_jobs=1, backend="threading")
    elif mode == "dynamical":
        print("\n Applying dynamical mode \n")
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode=mode, groupby=groupby, groups=groups, groups_for_fit=groups_for_fit, perc=perc, min_r2=min_r2)
        scv.tl.velocity_graph(adata)
    else: 
        print("\n Applying deterministic mode \n")
        scv.tl.velocity(adata, mode="deterministic", groupby=groupby, groups_for_fit=groups_for_fit, perc=perc, min_r2=min_r2)


    #####
    """Move the diffmap coordinates to the dataset"""
    diffmap_csv = pd.read_csv("/g/huber/users/kharuk/DECODE/adata_notch_control_3p_diffmap.csv", index_col=0)
    # let’s say Y is a pandas Df with diffusion components and the cell names as index.
    adata.obsm["X_diffmap"] = diffmap_csv.loc[adata.obs_names].values
    #####

    if moments_approach == "split":
     # plot some plots
    # compute and plot diffmap for notch + control samples
        sc.pl.diffmap(adata, color=annotation, components="2,3", title=f"Dmap for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments=split", ncols=1, frameon=True)
        scv.pl.velocity_embedding_stream(adata, components="1,2", basis="diffmap", color=annotation, dpi=140, legend_loc="right margin", title=f"Stream for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments split, condition={condition}")
        scv.pl.velocity_embedding_grid(adata, components="1,2", basis="diffmap", color=annotation, dpi=140, legend_loc="right margin", title=f"Stream for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments split, condition={condition}", min_mass=17, arrow_size=10, scale=.7)
    else:
     # plot some plots
    # compute and plot diffmap for notch + control samples
        sc.pl.diffmap(adata, color=annotation, components="2,3", title=f"Dmap for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments=default", ncols=1, frameon=True)
        scv.pl.velocity_embedding_stream(adata, components="1,2", basis="diffmap", color=annotation, dpi=140, legend_loc="right margin", title=f"Stream for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments default, condition={condition}")
        scv.pl.velocity_embedding_grid(adata, components="1,2", basis="diffmap", color=annotation, dpi=140, legend_loc="right margin", title=f"Stream for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments default, condition={condition}", min_mass=17, arrow_size=10, scale=.7)

    ### ranking genes 
    scv.tl.rank_velocity_genes(adata, groupby=annotation, min_corr=0.3)
    # # get the table of the top genes
    top_ranked_genes = scv.get_df(adata.uns["rank_velocity_genes"]["names"])
    
    return adata, top_ranked_genes


# function for RNA velocity estimation

def velocity_inferring_local(adata, min_shared_counts=None, min_shared_cells=None, n_top_genes=None, mode = "stochastic", condition=None,
                       annotation = None, moments_approach=None, groupby=None, groups=None, groups_for_fit=None, perc=[5,95], min_r2=0.01):
    if min_shared_counts is None:
        scv.pp.filter_genes(adata, min_shared_cells=min_shared_cells) 
    else:
        print("\n Applying min_shared_counts filter \n")
        scv.pp.filter_genes(adata, min_shared_counts=min_shared_counts)
        
    scv.pp.normalize_per_cell(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    sc.pp.log1p(adata)

    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    if "X_diffmap" in adata.obsm:
        del adata.obsm["X_diffmap"]
    
    print("'\n Computing PCA, knn and UMAP \n ")
    # compute the PCA
    scv.pp.pca(adata, n_comps=50) # 50 is the default value
    # compute the knn
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15) # knn 15 is the default value
    # compute the UMAP
    scv.tl.umap(adata, n_components=2, min_dist=0.8) # 2 is the default value
    print("\n Computing RNA velocity \n ")
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    if mode == "stochastic":
        print("\n Applying stochastic mode \n")
        scv.tl.velocity(adata, mode=mode, groupby=groupby, groups_for_fit=groups_for_fit, perc=perc, min_r2=min_r2)
        scv.tl.velocity_graph(adata)
    elif mode == "dynamical":
        print("\n Applying dynamical mode \n")
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode=mode, groupby=groupby, groups=groups, groups_for_fit=groups_for_fit, perc=perc, min_r2=min_r2)
        scv.tl.velocity_graph(adata)
    else: 
        print("\n Applying deterministic mode \n")
        scv.tl.velocity(adata, mode="deterministic", groupby=groupby, groups_for_fit=groups_for_fit, perc=perc, min_r2=min_r2)


    #####
    """Move the diffmap coordinates to the dataset"""
    diffmap_csv = pd.read_csv("/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/data/adata_notch_control_3p_diffmap.csv", index_col=0)
    # let’s say Y is a pandas Df with diffusion components and the cell names as index.
    adata.obsm["X_diffmap"] = diffmap_csv.loc[adata.obs_names].values
    #####

    if moments_approach == "split":
     # plot some plots
    # compute and plot diffmap for notch + control samples
        sc.pl.diffmap(adata, color=annotation, components="2,3", title=f"Dmap for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments=split", ncols=1, frameon=True)
        scv.pl.velocity_embedding_stream(adata, components="1,2", basis="diffmap", color=annotation, dpi=140, legend_loc="right margin", title=f"Stream for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments split, condition={condition}")
        scv.pl.velocity_embedding_grid(adata, components="1,2", basis="diffmap", color=annotation, dpi=140, legend_loc="right margin", title=f"Stream for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments split, condition={condition}", min_mass=17, arrow_size=10, scale=.7)
    else:
     # plot some plots
    # compute and plot diffmap for notch + control samples
        sc.pl.diffmap(adata, color=annotation, components="2,3", title=f"Dmap for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments=default", ncols=1, frameon=True)
        scv.pl.velocity_embedding_stream(adata, components="1,2", basis="diffmap", color=annotation, dpi=140, legend_loc="right margin", title=f"Stream for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments default, condition={condition}")
        scv.pl.velocity_embedding_grid(adata, components="1,2", basis="diffmap", color=annotation, dpi=140, legend_loc="right margin", title=f"Stream for 3' N + C; {mode} model, \n min_s_counts={min_shared_counts}, HVGs={n_top_genes}, \n groups_for_fit={groups_for_fit}, without groupby arg \n moments default, condition={condition}", min_mass=17, arrow_size=10, scale=.7)



def velocity_inferring_wo_plot(adata, min_shared_counts=None, min_shared_cells=None, n_top_genes=None, mode = "stochastic", condition=None,
                       annotation = None, moments_approach=None, groupby=None, groups=None, groups_for_fit=None, perc=[5,95], min_r2=0.01):
    """This function without plotting the diffmap"""
    if min_shared_counts is None:
        scv.pp.filter_genes(adata, min_shared_cells=min_shared_cells) 
    else:
        print("\n Applying min_shared_counts filter \n")
        scv.pp.filter_genes(adata, min_shared_counts=min_shared_counts)
        
    scv.pp.normalize_per_cell(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    sc.pp.log1p(adata)

    print("'\n Computing PCA, knn and UMAP \n ")
    # compute the PCA
    scv.pp.pca(adata, n_comps=50) # 50 is the default value
    # compute the knn
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15) # knn 15 is the default value
    # compute the UMAP
    scv.tl.umap(adata, n_components=2, min_dist=0.8) # 2 is the default value
    print("\n Computing RNA velocity \n ")
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    if mode == "stochastic":
        print("\n Applying stochastic mode \n")
        scv.tl.velocity(adata, mode=mode, groupby=groupby, groups_for_fit=groups_for_fit, perc=perc, min_r2=min_r2)
        scv.tl.velocity_graph(adata, n_jobs=1)
    elif mode == "dynamical":
        print("\n Applying dynamical mode \n")
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode=mode, groupby=groupby, groups=groups, groups_for_fit=groups_for_fit, perc=perc, min_r2=min_r2)
        scv.tl.velocity_graph(adata, n_jobs=1)
    else: 
        print("\n Applying deterministic mode \n")
        scv.tl.velocity(adata, mode="deterministic", groupby=groupby, groups_for_fit=groups_for_fit, perc=perc, min_r2=min_r2)

    ### ranking genes 
    scv.tl.rank_velocity_genes(adata, groupby=annotation, min_corr=0.3)
    # # get the table of the top genes
    top_ranked_genes = scv.get_df(adata.uns["rank_velocity_genes"]["names"])
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    # try to get rid of the warnings
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning, module="scvelo.plotting.utils")


    return adata, top_ranked_genes


# function for RNA velocity estimation with the split moments approach
def velocity_inference_moments(adata_notch, adata_control, min_shared_counts=20, n_top_genes=2000, n_pcs=30, n_neighbors=30, mode = "stochastic"):
    """The function will be just for one estimation of the RNA velocity
    Performing moments for each then merging the datasets"""
    # for notch 
    scv.pp.normalize_per_cell(adata_notch)   
    scv.pp.pca(adata_notch)
    scv.pp.moments(adata_notch, n_pcs=n_pcs, n_neighbors=n_neighbors)

    # for control 
    scv.pp.normalize_per_cell(adata_control)
    scv.pp.pca(adata_control)
    scv.pp.moments(adata_control, n_pcs=n_pcs, n_neighbors=n_neighbors) 

    adata_concat = ad.concat([adata_control, adata_notch], join="inner")  # join = outer and inner does not work for the diffmap and rank velocity genes
                                                                        # merge = same, unique, first. It does work after deleting the ["x_diffmap"]
    # filter after merging
    scv.pp.filter_and_normalize(adata_concat, min_shared_counts=min_shared_counts, n_top_genes=n_top_genes)
    # shufling the cells on the adata object

    obs_names = adata_concat.obs_names.copy()
    shuffled = np.random.choice(obs_names, len(obs_names), replace=False) # or so, you know what I mean
    adata_concat_shuffled  = adata_concat[shuffled].copy()

    # remove the X_diffmap to calculate RNA velocity
    del adata_concat_shuffled.obsm["X_diffmap"]

    # !!!!! I should run knn here bacause it wanted it 
    sc.pp.neighbors(adata_concat_shuffled)


    if mode == "stochastic":
        scv.tl.velocity(adata_concat_shuffled, mode=mode)
        scv.tl.velocity_graph(adata_concat_shuffled)
    else:
        scv.tl.velocity(adata_concat_shuffled, mode=mode)
        scv.tl.velocity_graph(adata_concat_shuffled)

    return adata_concat_shuffled


# function for Milo DA analysis with batch correction
def differential_abundance(adata):
    adata.obs["batch"] = adata.obs["sample_name"].str[:-2]
    # correct the data with harmony
    # run the Harmony correction for the nice diffmap
    scv.pp.pca(adata) # n_comps 50 is default

    sce.pp.harmony_integrate(adata, "batch", adjusted_basis="X_pca_harmony_batch",  max_iter_harmony = 50)
    # compute the knn and UMAP for this harmony data
    sc.pp.neighbors(adata, use_rep="X_pca_harmony_batch", key_added="neighbors_harmony_batch", n_neighbors=15) 
    sc.tl.umap(adata, neighbors_key="neighbors_harmony_batch")
    """Milo analyisis"""
    # initialize the object for Milo analysis
    milo = pt.tl.Milo()
    mdata = milo.load(adata)
    """compute KNN from Harmony corrected on batch"""
    ## Build KNN graph
    sc.pp.neighbors(mdata["rna"], use_rep="X_pca_harmony_batch", n_neighbors=30)
    ## Construct neighbourhoods
    milo.make_nhoods(mdata["rna"], prop=0.1) # ? I can specify here the parameter for neighbors_key, it is specified in the pp.neighbors function
    mdata["rna"].obsm["nhoods"]

    ## visualize the distribuition of neighbourhood sizes
    nhood_size = np.array(mdata["rna"].obsm["nhoods"].sum(0)).ravel()
    plt.figure()
    plt.hist(nhood_size, bins=100)
    plt.xlabel("# cells in nhood")
    plt.ylabel("# nhoods")
    ## Count cells in neighbourhoods
    mdata = milo.count_nhoods(mdata, sample_col = "sample_name")
    
    ### Differential abundance testing with GLM
    # Reorder categories for the perturbation
    # by default, the last category is taken as the condition of interest
    mdata["rna"].obs["perturbation"] = mdata["rna"].obs["perturbation"].cat.reorder_categories(["control", "Notch gRNA"])
    milo.da_nhoods(mdata, design="~perturbation")
    plt.figure()
    old_figsize = plt.rcParams["figure.figsize"]
    plt.rcParams["figure.figsize"] = [5, 5]
    # plt.subplot(1, 2, 1)
    plt.hist(mdata["milo"].var.PValue, bins=50)
    plt.xlabel("P-Vals")
    # plt.subplot(1, 2, 2)
    plt.figure()
    plt.plot(mdata["milo"].var.logFC, -np.log10(mdata["milo"].var.SpatialFDR), ".")
    plt.xlabel("log-Fold Change")
    plt.ylabel("- log10(Spatial FDR)")
    plt.tight_layout()
    plt.rcParams["figure.figsize"] = old_figsize
    
    # assing the diffmap coordinates to new obsm for plotting
    mdata["rna"].obsm["X_diffmap_2_3"] = mdata["rna"].obsm["X_diffmap"][:, [1, 2]]
    ### Visualize results on embedding
    milo.build_nhood_graph(mdata, basis="X_diffmap_2_3")
    plt.rcParams["figure.figsize"] = [7, 5]
    milo.plot_nhood_graph(
    mdata,
    alpha=0.1,  ## SpatialFDR level (10%)
    min_size=1,  ## Size of smallest dot
    title="DA logFC with Spatial FDR 10%")
    # assign the milo results and plot them on the diffmap 
    mdata_plot = mdata["milo"].T.copy()
    mdata_plot.obs["color_logFC"] = mdata_plot.obs["logFC"]
    mdata_plot.obs["color_PValue"] = mdata_plot.obs["PValue"]
    mdata_plot.obs["color_FDR"] = mdata_plot.obs["FDR"]
    mdata_plot.obs["color_SpatialFDR"] = mdata_plot.obs["SpatialFDR"]
    mdata_plot.obs["color_Nhood_size"] = mdata_plot.obs["Nhood_size"]
    # plot all the parameters on the diffmap
    sc.pl.embedding(mdata_plot, basis="X_milo_graph", color=["color_logFC", "color_PValue", "color_FDR", "color_SpatialFDR", "color_Nhood_size"], frameon=False, ncols=1)
    ### Visualize result by celltype
    milo.annotate_nhoods(mdata, anno_col="high_res_annotation")    
    plt.figure()
    plt.hist(mdata["milo"].var["nhood_annotation_frac"], bins=30)
    plt.xlabel("celltype fraction")
    # plot the celltype log fold change
    plt.figure()
    milo.plot_da_beeswarm(mdata, alpha=0.1)
    return adata, milo, mdata

def plot_qc_metrics(adata, layer=None, seq_technique=None):
    if layer == "spliced":
        adata = adata.copy()
        adata.X = adata.layers["spliced"]
        ## calculate the number of counts in the cell to obs
        adata.obs["ncounts"] = np.ravel(np.sum(adata.X, axis=1))
        ## calculate the number of genes in the cell to obs
        adata.obs["ngenes"] = np.ravel(np.sum(adata.X > 0, axis=1))
        # perform mito and ribo counts 
        adata.var["mt"] = adata.var_names.str.startswith("mt:")
        adata.var["ribo"] = adata.var_names.str.startswith(("RpS", "RpS"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None,
                            log1p=False, inplace=True)
        plot_qc_metrics_i(adata, layer=layer, seq_technique=seq_technique)

    elif layer =="unspliced":
        adata = adata.copy()
        adata.X = adata.layers["unspliced"]
        ## calculate the number of counts in the cell to obs
        adata.obs["ncounts"] = np.ravel(np.sum(adata.X, axis=1))
        ## calculate the number of genes in the cell to obs
        adata.obs["ngenes"] = np.ravel(np.sum(adata.X > 0, axis=1))
        # perform mito and ribo counts 
        adata.var["mt"] = adata.var_names.str.startswith("mt:")
        adata.var["ribo"] = adata.var_names.str.startswith(("RpS", "RpS"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None,
                            log1p=False, inplace=True)
        plot_qc_metrics_i(adata, layer=layer, seq_technique=seq_technique)

    elif layer =="counts":
        adata = adata.copy()
        adata.X = adata.layers["counts"]
        ## calculate the number of counts in the cell to obs
        adata.obs["ncounts"] = np.ravel(np.sum(adata.X, axis=1))
        ## calculate the number of genes in the cell to obs
        adata.obs["ngenes"] = np.ravel(np.sum(adata.X > 0, axis=1))
        # perform mito and ribo counts 
        adata.var["mt"] = adata.var_names.str.startswith("mt:")
        adata.var["ribo"] = adata.var_names.str.startswith(("RpS", "RpS"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None,
                            log1p=False, inplace=True)
        plot_qc_metrics_i(adata, layer=layer, seq_technique=seq_technique)
        
        """PLOT ALL QC METRICS"""

        """Number of genes, total counts, pct count mito amd ribo for ALL samples"""
def plot_qc_metrics_i(adata, layer=None, seq_technique=None):
        plt.figure(figsize=(4, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="n_genes_by_counts", fill=False, width=0.4, gap=3)
        # xlabel
        plt.xlabel("All samples")
        plt.ylabel("Number of genes")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_ngenes.pdf", bbox_inches='tight')
        plt.show()
        # plot several QC metrics for N + control samples 3' using seaborn

        # make the violin more 
        plt.figure(figsize=(4, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="total_counts", fill=False, width=0.4, gap=3)
        # xlabel
        plt.xlabel("All samples")
        plt.ylabel("Total counts")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_counts.pdf", bbox_inches='tight')
        plt.show()
        # plot several QC metrics for N + control samples 3' using seaborn

        # make the violin more 
        plt.figure(figsize=(4, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="pct_counts_mt", fill=False, width=0.4, gap=3)
        # xlabel
        plt.xlabel("All samples")
        plt.ylabel("Percentage of miochondrial counts")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_pct_mito.pdf", bbox_inches='tight')
        plt.show()
        # plot several QC metrics for N + control samples 3' using seaborn

        # make the violin more 
        plt.figure(figsize=(4, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="pct_counts_ribo", fill=False, width=0.4, gap=3)
        # xlabel
        plt.xlabel("All samples")
        plt.ylabel("Percentage of ribosomal counts")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_pct_ribo.pdf", bbox_inches='tight')
        plt.show()
        # plot several X QC metrics for N + control samples 3' using seaborn


        """Number of genes, total counts, pct count mito and ribo for per condition"""

        plt.figure(figsize=(5, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="n_genes_by_counts", x="perturbation", fill=False, width=0.4, gap=3)
        # xlabel
        plt.xlabel("Condition")
        plt.ylabel("Number of genes")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_ngenes_condition.pdf", bbox_inches='tight')
        sns.despine()
        plt.show()

        # make the violin more 
        plt.figure(figsize=(5, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="total_counts", x="perturbation", fill=False, width=0.4, gap=3)
        # xlabel
        plt.xlabel("Condition")
        plt.ylabel("Total counts")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_counts_condition.pdf", bbox_inches='tight')
        plt.show()

        plt.figure(figsize=(5, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="pct_counts_mt", x="perturbation", fill=False, width=0.4, gap=3)
        # xlabel
        plt.xlabel("Condition")
        plt.ylabel("Percentege of miochondrial counts")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_pct_mt_condition.pdf", bbox_inches='tight')
        plt.show()

        plt.figure(figsize=(5, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="pct_counts_ribo", x="perturbation", fill=False, width=0.4, gap=3)
        # xlabel
        plt.xlabel("Condition")
        plt.ylabel("Percentege of ribosomal counts")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_pct_ribo_condition.pdf", bbox_inches='tight')
        plt.show()


        """Number of genes and total counts, pct mito and ribo per cluster"""

        fig, axs0 = plt.subplots(1, 1, figsize=(12, 5))  # Adjust figsize as needed
        sc.pl.violin(adata, ["n_genes_by_counts"], jitter=0.4, show=False, ncols = 1, groupby = "high_res_annotation", rotation=90, ax=axs0) # do not have to specify the layer as it already calculated for specific layers
        plt.suptitle(f"Quality control metrics for {layer} counts", y=1)
        axs0.spines['right'].set_visible(False)
        axs0.spines['top'].set_visible(False)
        fig.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_ngenes_cluster.pdf", bbox_inches='tight')
        plt.show()

        fig, axs0 = plt.subplots(1, 1, figsize=(12, 5))  # Adjust figsize as needed
        sc.pl.violin(adata, ["total_counts"], jitter=0.4, show=False, ncols = 1, groupby = "high_res_annotation", rotation=90, ax=axs0) # do not have to specify the layer as it already calculated for specific layers
        plt.suptitle(f"Quality control metrics for {layer} counts", y=1)
        axs0.spines['right'].set_visible(False)
        axs0.spines['top'].set_visible(False)
        fig.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_counts_cluster.pdf", bbox_inches='tight')
        plt.show()

        fig, axs0 = plt.subplots(1, 1, figsize=(12, 5))  # Adjust figsize as needed
        sc.pl.violin(adata, ["pct_counts_ribo"], jitter=0.4, show=False, ncols = 1, groupby = "high_res_annotation", rotation=90, ax=axs0) # do not have to specify the layer as it already calculated for specific layers
        plt.suptitle(f"Quality control metrics for {layer} counts", y=1)
        axs0.spines['right'].set_visible(False)
        axs0.spines['top'].set_visible(False)
        fig.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_pct_ribo_cluster.pdf", bbox_inches='tight')
        plt.show()

        fig, axs0 = plt.subplots(1, 1, figsize=(12, 5))  # Adjust figsize as needed
        sc.pl.violin(adata, ["pct_counts_mt"], jitter=0.4, show=False, ncols = 1, groupby = "high_res_annotation", rotation=90, ax=axs0) # do not have to specify the layer as it already calculated for specific layers
        plt.suptitle(f"Quality control metrics for {layer} counts", y=1)
        axs0.spines['right'].set_visible(False)
        axs0.spines['top'].set_visible(False)
        fig.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_pct_mt_cluster.pdf", bbox_inches='tight')
        plt.show()

        """Number of genes and total counts, pct mito and ribo per sample"""
        plt.figure(figsize=(8, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="n_genes_by_counts", x="sample_name", fill=False, width=0.4, gap=3)
        plt.xlabel("Sample")
        plt.ylabel("Number of genes")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_ngenes_sample.pdf", bbox_inches='tight')
        plt.show()

        plt.figure(figsize=(8, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="total_counts", x="sample_name", fill=False, width=0.4, gap=3)
        plt.xlabel("Sample")
        plt.ylabel("Total counts")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_counts_sample.pdf", bbox_inches='tight')
        plt.show()

        plt.figure(figsize=(8, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="pct_counts_ribo", x="sample_name", fill=False, width=0.4, gap=3)
        plt.xlabel("Sample")
        plt.ylabel("Percentage of ribosomal counts")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_pct_ribo_sample.pdf", bbox_inches='tight')
        plt.show()

        plt.figure(figsize=(8, 5))  # Set width to 5, height to 6
        sns.violinplot(adata.obs, y="pct_counts_mt", x="sample_name", fill=False, width=0.4, gap=3)
        plt.xlabel("Sample")
        plt.ylabel("Percentage of mito counts")
        plt.title(f"{layer}")
        sns.despine()
        plt.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_pct_mt_sample.pdf", bbox_inches='tight')
        plt.show()


        """Scatter plots pct mito vs total counts per cluster, ngenes vs total counts per cluster, ngenes vs total counts per pct mito counts"""

        ## scatter plots
        # Create a figure and axes for subplots
        fig, axs = plt.subplots(1, 3, figsize=(20, 7))  # Adjust figsize as needed
        # Plot each scatter plot in a separate subplot
        sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", color="high_res_annotation", ax=axs[0], show=False, legend_loc="best", legend_fontsize="10", frameon=False)
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="high_res_annotation", ax=axs[1], show=False, legend_loc="none", frameon=False)
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", ax=axs[2], show=False, frameon=False)
        # Add a title to the figure
        fig.suptitle(f"Quality control metrics for {layer} counts")
        # set ticks for X axis
        axs[0].set_xticks(np.arange(0, 50000, 10000))
        axs[1].set_xticks(np.arange(0, 50000, 10000)) 
        axs[2].set_xticks(np.arange(0, 50000, 10000)) 

        # remove the right and top spines
        for ax in axs:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

        # Move the legend
        fig.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust the rectangle to move the legend
        fig.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_scatter_3.pdf", bbox_inches='tight')
        plt.show()

        """ FIGURE 2 """
        # plot the number of counts vs the number of genes for perturbation and sample for X matrix
        fig, axes = plt.subplots(2,2, figsize=(20,15))
        fig.suptitle(f"QC for {layer} matrix, adata.obs")
        # plot the number of counts vs the number of genes for perturbation
        sns.scatterplot(ax=axes[0,0], data = adata.obs, x="total_counts", y="n_genes_by_counts", hue="perturbation")
        sns.despine()

        sns.scatterplot(ax=axes[1,0], data = adata.obs, x="total_counts", y="n_genes_by_counts", hue="sample_name")
        sns.despine()

        sns.scatterplot(ax=axes[0,1], data = adata.obs, x="total_counts", y="n_genes_by_counts", hue="high_res_annotation")
        sns.despine()
        axes[0,1].legend(bbox_to_anchor=(1.5, 1), borderaxespad=0)

        # add new column with the batches to obs just with 22 and 23
        # remove 2 last letters from the each row in sample name 
        adata.obs["batch"] = adata.obs["sample_name"].str[:-2]
        sns.scatterplot(ax=axes[1,1], data = adata.obs, x="total_counts", y="n_genes_by_counts", hue="batch")
        sns.despine()
        fig.savefig(f"/Users/sviatoslavkharuk/Documents/StudyBioinformatics/Bioinformatics_EMBL/Bioinformatics_project_EMBL/scRNA-seq-RNA-velocity-CRISPR-DECODE/Quarto-notebook-DECODE/figures/qc/{layer}_QC_metrics_N_control_{seq_technique}_scatter_4.pdf", bbox_inches='tight')
        plt.legend(bbox_to_anchor=(1.5, 1), borderaxespad=0)
        