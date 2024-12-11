Functions
-----------

The main functions in RidgeSpace are listed below:

.. code-block:: python

    RidgeSpace.tl_mesh(adata=None, obsm_key='spatial', input_x=None, input_y=None, max_radius=None, plot_kde=False, mesh_directory='./Ridge_tmp', save_address='', mesh_optimization=True)  

This function used for 3D mesh construction. 


.. code-block:: python

    RidgeSpace.tl_denoise(adata=None, plot_name=None, find_obs=False, use_raw=True, input_value=None, mesh_directory='./Ridge_tmp', save_address='', save_directory=None, iteration=10, L1=2/3, weight=0.9) 

This function used for spatial data denoising. 

.. code-block:: python

    RidgeSpace.tl_HE(adata=None, library_id=None, image_key='lowres', scalefactor='tissue_lowres_scalef', HE_directory='./Ridge_tmp', save_address='', filter_background=False, filter_spot=None, image_mask=None, threshold_scale=1.1, lightness=1.0, fast=True) 

This function used for HE image processing. 


.. code-block:: python

    RidgeSpace.pl_single(adata=None, find_obs=False, use_raw=True, input_z=None, plot_name='Gene', normalize_z=True, normlabels=False, out_path=None, obs_cluster=None, color_map=None, cluster=None, mesh_directory='./Ridge_tmp', save_address='', lw=0, height_ratio=1.5, elev=35, view=-70, z_height=1, bg_alpha=1, color_brighten=True, plot_HE=False, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, HE_slice=1, plot_clustering=False, clustering_z=None, clustering_size=1, clustering_alpha=0.1, plot_clustering_boundary=False, clustering_boundary_lw=1, clustering_boundary_alpha=0.4, select_c=None, xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, trunc=None, figsize_set=(9,8), return_ax=False) 

This function used for basic spatial visualization (single).

.. code-block:: python

    RidgeSpace.pl_multipleOUT(adata=None, find_obs=False, use_raw=True, input_zA=None, input_zB=None, plot_nameA='Gene', plot_nameB='Gene', normalize_z=True, normlabels=True, out_path=None, obs_cluster=None, color_map=None, cluster=None, mesh_directory='./Ridge_tmp', save_address='', lw=0, elev=35, view=-70, z_height=1, bg_alphaA=1, bg_alphaB=1, color_brighten=True, plot_HE=False, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, HE_slice=1, plot_clustering=False, clustering_z=None, clustering_size=1, clustering_alpha=0.1, select_c=None, xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, truncA=None, truncB=None, figsize_set=(8,8), return_ax=False)

This function used for multi-modal spatial comparison. 

.. code-block:: python

    RidgeSpace.pl_multipleIN(adata=None, find_obs=False, use_raw=True, input_zA=None, input_zB=None, plot_nameA='Gene', plot_nameB='Gene', height=100, normalize_z=True, normlabels=False, out_path=None, obs_cluster=None, color_map=None, cluster=None, mesh_directory='./Ridge_tmp', save_address='', lw=0, elev=35, view=-70, z_height=1, bg_alphaA=1, bg_alphaB=1, color_brighten=True, plot_HE=False, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, HE_slice=1, plot_clustering=False, clustering_z=None, clustering_size=1, clustering_alpha=0.1, plot_clustering_boundary=False, clustering_boundary_lw=1, clustering_boundary_alpha=0.4, select_c=None,xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, truncA=None, truncB=None, figsize_set=(8,8), return_ax=False) 

This function used for multi-modal spatial comparison.

.. code-block:: python

    RidgeSpace.pl_trajectory(adata=None, input_z=None, find_obs=True, use_raw=True, normalize_z=False, normlabels=True, out_path=None, obs_cluster=None, color_map=None, cluster=None, mesh_directory='./Ridge_tmp', save_address='', lw=0, height_ratio=1.5, elev=35, view=-70, z_height=1, bg_alpha=1, color_brighten=True, plot_name='Gene', plot_HE=True, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, HE_slice=1, select_c=None, xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, trunc=None, figsize_set=(8,8), density_arrow = 0.4, dense=  True, min_length = 0.15, color_arrow="k", lw_arrow=1, arrow_scale=20, return_ax=False) 

This function used for plot spatial trajectory. 
