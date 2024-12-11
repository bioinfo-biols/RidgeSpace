from .utils import *
import os
import matplotlib

def RidgeSpace_command():
    import sys
    ''' Example of taking inputs for SEVtras'''
    args = sys.argv[1:]
    if len(args) < 1:
        print("usage: convert_zarr, tl_mesh, tl_HE, tl_denoise, pl_single, pl_multipleOUT, pl_multipleIN, pl_trajectory")


def tl_mesh(adata=None, obsm_key='spatial', input_x=None, input_y=None, max_radius=None, plot_kde=False, mesh_directory='./Ridge_tmp', save_address='', mesh_optimization=True):
    import os
    import pickle

    if not os.path.exists(mesh_directory):
        os.makedirs(mesh_directory)

    if save_address != '':
        save_directory = mesh_directory + '/' + save_address + '/'
    else:
        save_directory = mesh_directory

    if adata is not None:
        x = adata.obsm[obsm_key][:,0]
        y = adata.obsm[obsm_key][:,1]
    else:
        x = input_x
        y = input_y

    triang = mtri.Triangulation(x, y)#CubicTriInterpolator(x,y)#
    triang.is_delaunay=False

    # plot only triangles with sidelength smaller some max_radius
    triangles = triang.triangles

    # Mask off unwanted triangles.
    xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
    ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
    maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)

    if max_radius is None:
        from collections import Counter
        Mode_r = Counter(maxi).most_common(1)[0][0] 
        max_radius1 = Mode_r + 0.02 * Mode_r
    else:
        max_radius1 = max_radius

    if plot_kde == True:
        try:
            import seaborn as sns
            sns.set(style="ticks", font_scale=2)
            sns.kdeplot(pd.DataFrame(maxi), gridsize=2000)
            plt.xlim(0, 100)
            plt.vlines(max_radius1, ymin=0, ymax=0.06, colors='black')
            sns.despine()
            plt.tight_layout()
            plt.show()
        except ImportError:
            raise ImportError("Plotting kde requires the `seaborn` module to be installed.")
    else:
        print('Radius   Count')
        print(pd.DataFrame(maxi).value_counts().head(5))
        print('Use default radius ' + str(max_radius1))

    triang.set_mask(maxi > max_radius1)
    import pickle
    if not os.path.exists(str(save_directory)):
        os.makedirs(str(save_directory))

    with open(str(save_directory) + "/triang.pickle", "wb") as output_file:
        pickle.dump([max_radius1, triang], output_file)

    print('Triangles mesh saved in ' + str(save_directory))

    if mesh_optimization == True:
        mesh_optimize(triang=triang, save_directory=save_directory)



def tl_HE(adata=None, library_id=None, image_key='lowres', scalefactor='tissue_lowres_scalef', HE_directory='./Ridge_tmp', save_address='', filter_background=False, filter_spot=None, image_mask=None, threshold_scale=1.1, lightness=1.0, fast=True):
    import os
    import pickle
    import numpy as np

    if not os.path.exists(HE_directory):
        os.makedirs(HE_directory)

    if save_address != '':
        save_directory = HE_directory + '/' + save_address + '/'
    else:
        save_directory = HE_directory

    if adata is None:
        print('Please provide anndata with histology images in adata.uns!')
        return 0
    
    if library_id is None:
        library_id = list(adata.uns['spatial'].keys())[0]
        print('Using ' + str(library_id))

    lowres = adata.uns["spatial"][library_id]['images'][image_key]
    rows, cols = lowres.shape[:2]

    if filter_background == True:
        if image_mask is None:
            try:
                from skimage import data,filters
                import cv2
                if np.max(lowres) > 10:
                    image_mask, contours =get_HE_mask((lowres).astype('uint8'), threshold_scale)
                else:
                    image_mask, contours =get_HE_mask((lowres*255).astype('uint8'), threshold_scale)
            except ImportError:
                raise ImportError("This method requires the `skimage` and `cv2` module to be installed. Or you can provide the binary mask of corresponding hitology image with parameter `image_mask`.")
    else:
        image_mask = np.ones((rows, cols))

    if fast == True:
        scalefactor_value = adata.uns['spatial'][library_id]['scalefactors'][scalefactor]

        sc_y, sc_x = np.mgrid[:rows, :cols] / scalefactor_value
        valid_mask = image_mask > 0
        sc_y = sc_y[valid_mask]
        sc_x = sc_x[valid_mask]
        lowres_filtered = lowres[valid_mask]

        if lightness != 1.0:
            import matplotlib.colors as mc
            # import colorsys
            # 调整亮度（矢量化操作）
            adjusted_colors = adjust_lightness_batch(lowres_filtered, lightness)

            # 批量转换 RGB 到十六进制颜色
            sc_z = np.apply_along_axis(lambda rgb: mc.to_hex(rgb, keep_alpha=False), 1, adjusted_colors)

        else:
            sc_z = lowres_filtered
    
    else:
        sc_x = []
        sc_y = []
        sc_z = []
        for i in range(lowres.shape[0]):
            for j in range(lowres.shape[1]):
                if image_mask[i][j] > 0:
                    sc_y.append(i / adata.uns['spatial'][library_id]['scalefactors'][scalefactor])
                    sc_x.append(j / adata.uns['spatial'][library_id]['scalefactors'][scalefactor])
                    if lightness != 1.0:
                        sc_z.append(adjust_lightness(matplotlib.colors.to_hex(lowres[i][j]), lightness))
                    else:
                        sc_z.append(matplotlib.colors.to_hex(lowres[i][j]))
        sc_x = np.array(sc_x)
        sc_y = np.array(sc_y)

    # filter_he = (sc_x > adata.obsm['spatial'].min(0)[0]) & (sc_x < adata.obsm['spatial'].max(0)[0]) & (sc_y > adata.obsm['spatial'].min(0)[1]) & (sc_y < adata.obsm['spatial'].max(0)[1])
    # 空间过滤条件
    spatial_min = adata.obsm['spatial'].min(0)
    spatial_max = adata.obsm['spatial'].max(0)
    filter_he = (
        (sc_x > spatial_min[0]) & (sc_x < spatial_max[0]) &
        (sc_y > spatial_min[1]) & (sc_y < spatial_max[1])
    )

    if filter_spot is None:
        sc_x = sc_x[filter_he]
        sc_y = sc_y[filter_he]
        sc_z = np.array(sc_z)[filter_he]
        save_list = [sc_x, sc_y, sc_z]

        import pickle
        if not os.path.exists(str(save_directory)):
            os.makedirs(str(save_directory))
        with open(str(save_directory) + "/HE.pickle", "wb") as output_file:
            pickle.dump(save_list, output_file)
    else:
        sc_x = sc_x[filter_he & filter_spot]
        sc_y = sc_y[filter_he & filter_spot]
        sc_z = np.array(sc_z)[filter_he & filter_spot]
        save_list = [sc_x, sc_y, sc_z]

        import pickle
        if not os.path.exists(str(save_directory)):
            os.makedirs(str(save_directory))
        with open(str(save_directory) + "/HE_filter.pickle", "wb") as output_file:
            pickle.dump(save_list, output_file)
    
    print('HE information saved in ' + str(save_directory))



def tl_denoise(adata=None, plot_name=None, find_obs=False, use_raw=True, input_value=None, mesh_directory='./Ridge_tmp', save_address='', save_directory=None, iteration=10, L1=2/3, weight=0.9):
    import os
    import pickle
    if save_address != '':
        mesh_directory = mesh_directory + '/' + save_address + '/'
    else:
        mesh_directory = mesh_directory
    
    if os.path.isfile(str(mesh_directory) + "/level1_list.pickle"):
        with open(str(mesh_directory) + "/level1_list.pickle", "rb") as input_file:
            level1_list = pickle.load(input_file)
            
        with open(str(mesh_directory) + "/level2_list.pickle", "rb") as input_file:
            level2_list = pickle.load(input_file)
    else:
        print("Please run RidgeSpace.tl_mesh() with parameter 'mesh_optimization=True' first.")
        return 0

    if save_directory is None:
        save_directory = mesh_directory
    else:
        if not os.path.exists(str(save_directory)):
            os.makedirs(str(save_directory))

        if save_address != '':
            save_directory = save_directory + '/' + save_address + '/'
        else:
            save_directory = save_directory

    if adata is None:
        original_value = input_value
        if input_value is None:
            print('Please provide molecular profile through `adata` or `input_value`')
    else:
        if find_obs == True:
            original_value = adata.obs[plot_name].values
        else:
            from scipy.sparse import csr_matrix
            if use_raw == True:
                matrix = csr_matrix(adata.raw.X)
                var_names = adata.raw.var_names
            else:
                matrix = csr_matrix(adata.X)
                var_names = adata.var_names

            original_value = matrix[:,var_names == plot_name].toarray().ravel()

    denoised_z = molecular_denoise(original_value, level1_list, level2_list, iteration=iteration, L1=L1, weight=weight)

    if plot_name is None:
        plot_name = 'Gene'

    import pickle
    if not os.path.exists(str(save_directory)):
        os.makedirs(str(save_directory))
    with open(str(save_directory) + "/" + str(plot_name) + "_denoise.pickle", "wb") as output_file:
        pickle.dump(denoised_z, output_file)

    if plot_name != 'Gene':
        print('Denoised value for ' + str(plot_name) + ' saved in ' + str(save_directory))
    else:
        print('Denoised value saved in ' + str(save_directory))
    

def pl_single(adata=None, find_obs=False, use_raw=True, input_z=None, plot_name='Gene', normalize_z=True, normlabels=False, out_path=None, obs_cluster=None, color_map=None, cluster=None, mesh_directory='./Ridge_tmp', save_address='', lw=0, height_ratio=1.5, elev=35, view=-70, z_height=1, bg_alpha=1, color_brighten=True, plot_HE=False, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, HE_slice=1, plot_clustering=False, clustering_z=None, clustering_size=1, clustering_alpha=0.1, plot_clustering_boundary=False, clustering_boundary_lw=1, clustering_boundary_alpha=0.4, select_c=None, xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, trunc=None, figsize_set=(9,8), return_ax=False):
    
    import pickle
    import os

    if save_address != '':
        mesh_directory = mesh_directory + '/' + save_address + '/'
    else:
        mesh_directory = mesh_directory

    if save_address != '':
        HE_directory = HE_directory + '/' + save_address + '/'
    else:
        HE_directory = HE_directory

    if input_z is None:
        files = os.listdir(mesh_directory)
        if str(plot_name) + "_denoise.pickle" in files:
            with open(str(mesh_directory) + "/" + str(plot_name) + "_denoise.pickle", "rb") as input_file:
                input_z = pickle.load(input_file)
        else:
            print("No file named " + str(plot_name) + " in " + str(mesh_directory))
            if adata is not None:
                print('Using noisy information in adata')
                if find_obs == True:
                    input_z = adata.obs[plot_name]
                else:
                    from scipy.sparse import csr_matrix
                    if use_raw == True:
                        matrix = csr_matrix(adata.raw.X)
                        var_names = adata.raw.var_names
                    else:
                        matrix = csr_matrix(adata.X)
                        var_names = adata.var_names

                    input_z = matrix[:,var_names == plot_name].toarray().ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_z = input_z

    import copy
    input_zP = input_z.copy()
    normalize_factor = 1
    normalize_min = 0
    if normalize_z == True:
        normalize_factor = 100 / max(input_z)
        normalize_min = np.min(input_z)
        norm_z = normalize_factor * (input_z - normalize_min) #z / np.linalg.norm(z)
    else:
        norm_z = input_z

    if trunc is not None:
        tmp = []
        for mv in norm_z:
            if mv > trunc:
                tmp.append(mv-trunc)
            else:
                tmp.append(0)
        norm_z = np.array(tmp)

    ##read triang mesh
    import pickle
    if os.path.isfile(str(mesh_directory) + "/triang.pickle"):
        with open(str(mesh_directory) + "/triang.pickle", "rb") as input_file:
            max_r, triang = pickle.load(input_file)
    else:
        print('No triangular mesh detected. Please run RidgeSpace.tl_mesh() first.')
        return 0

    if adata is not None:
        cluster = adata.obs[obs_cluster].values
        uns_cluster = obs_cluster + '_colors'
        if uns_cluster in adata.uns.keys():
            cluster_colors = adata.uns[uns_cluster]
        else:
            cluster_colors = None

    if cluster_colors is None:
        if color_map is None:
            cluster_colors = ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#e377c2', '#8c564b', '#17becf', '#b5bd61', '#aec7e8', '#aa40fc', '#ffbb78', '#98df8a','#ff9896', '#11e6df', '#e9812d', '#e5cb83', '#b2c3d8', '#bee6da', '#cccccc', '#906bac', '#f51d1f', '#c49d98']
        else:
            cluster_colors = color_map

    try:
        tmp = cluster.categories
    except AttributeError:
        cluster = pd.Categorical(cluster)

    colors, original_colors = color_alpha(input_z=norm_z, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alpha, color_brighten=color_brighten)

    fig = plt.figure( figsize=figsize_set)#
    try:
        ax = fig.add_subplot(1,1,1, projection='3d', computed_zorder=False)
    except:
        ax = fig.add_subplot(1,1,1, projection='3d')
        print('Using computed zorder!')
    

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['pdf.use14corefonts'] = True
    # plt.rcParams['font.family'] = 'sans-serif'
    # plt.rcParams['font.sans-serif'] = 'Arial'

    # matplotlib.rcParams['axes.formatter.useoffset'] = False

    if select_c is not None:
        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colors])
        colors = colors_new

    plot_Ridge(ax, triang, norm_z, lw = lw, facecolors = colors, antialiased=False, rasterized=True, zorder=1)

    if plot_HE == True:
        import pickle
        import os
        if HE_filter == False:
            if os.path.isfile(str(HE_directory) + "/HE.pickle"):
                with open(str(HE_directory) + "/HE.pickle", "rb") as input_file:
                    sc_x, sc_y, sc_z = pickle.load(input_file)
            else:
                print('Please provide hitology image and run RidgeSpace.tl_HE()')
                return 0
                
        elif HE_filter == True:
            if os.path.isfile(str(HE_directory) + "/HE_filter.pickle"):
                with open(str(HE_directory) + "/HE_filter.pickle", "rb") as input_file:
                    sc_x, sc_y, sc_z = pickle.load(input_file)
            else:
                print('Please provide hitology image and run RidgeSpace.tl_HE()')
                return 0
                
        HE_slice = int(HE_slice)
        if HE_slice < 1:
            HE_slice = 1
        if HE_z is None:
            HE_z = height_ratio*max(norm_z)
            ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, zorder=3, rasterized=True)
        else:
            if height_ratio*max(norm_z) - HE_z >  HE_z:
                ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=0)
            else:
                ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=3)

        
    if plot_clustering == True:
        if clustering_z is None:
            ax.scatter3D(xs=triang.x, ys=triang.y, zs=height_ratio*max(norm_z), s=clustering_size, facecolors=np.array(original_colors), alpha=clustering_alpha, antialiased=False, linewidth=0, zorder=4, rasterized=True)#
        else:
            ax.scatter3D(xs=triang.x, ys=triang.y, zs=clustering_z, s=clustering_size, facecolors=np.array(original_colors), alpha=clustering_alpha, antialiased=False, linewidth=0, rasterized=True, zorder=0)#

    if plot_clustering_boundary == True:

        for colori in pd.Series(original_colors).unique():
            Tx = triang.x[(original_colors == colori)]
            Ty = triang.y[(original_colors == colori)]
            T = mtri.Triangulation(Tx, Ty)#
            T.is_delaunay=False
            # plot only triangles with sidelength smaller some max_radius
            Ttriangles = T.triangles

            # Mask off unwanted triangles.
            Txtri = Tx[Ttriangles] - np.roll(Tx[Ttriangles], 1, axis=1)
            Tytri = Ty[Ttriangles] - np.roll(Ty[Ttriangles], 1, axis=1)
            Tmaxi = np.max(np.sqrt(Txtri**2 + Tytri**2), axis=1)

            T.set_mask(Tmaxi > max_r)

            # from Stack Overflow
            boundaries = []
            for i in range(len(T.triangles)):
                if T.mask[i] == False:
                    for j in range(3):
                        if T.neighbors[i,j] < 0:
                            # Triangle edge (i,j) has no neighbor so is a boundary edge.
                            boundaries.append((T.triangles[i,j],
                                            T.triangles[i,(j+1)%3]))
            boundaries = np.asarray(boundaries)

            # The following line plots the boundary edges.
            bheight = height_ratio*max(norm_z)
            for i in range(len(boundaries)):
                ax.plot([Tx[boundaries].T[0,i], Tx[boundaries].T[1,i]], [Ty[boundaries].T[0,i], Ty[boundaries].T[1,i]], \
                        [bheight, bheight], color=colori, lw=clustering_boundary_lw, alpha=clustering_boundary_alpha)#

    ax.view_init(elev=elev, azim=view)


    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_zlim(zlim_min, zlim_max)

    if normlabels == False:
        if trunc is not None:
            # input_zP norm_z normalize_min normalize_factor
            input_z_min, input_z_max = max(input_zP) * trunc / 100, max(input_zP)
            start, stop, step = get_plotlabel_trunc(0.95 * input_z_min, 1.05 * input_z_max)
            plot_labels = []
            for labels_i in np.arange(start, stop, step):
                plot_labels.append((labels_i - normalize_min) * normalize_factor - trunc)
            
        else:
            # input_zP norm_z normalize_min normalize_factor
            # plot_labels = ax.get_ztickslabel()
            input_z_min, input_z_max = min(input_zP), max(input_zP)
            start, stop, step = get_plotlabel(0.95 * input_z_min, 1.05 * input_z_max)
            plot_labels = []
            for labels_i in np.arange(start, stop, step):
                plot_labels.append((labels_i - normalize_min) * normalize_factor)
        
        ax.set_zticks(plot_labels)
        if step > 50:
            baseS = np.round(np.log10(step))
            baseSN = 10 ** (baseS - 1)
            ax.set_zticklabels([ format_number(tmpi) for tmpi in np.arange(start / baseSN, stop / baseSN, step / baseSN)] )
            ax.set_zlabel(str(plot_name) + ' expression (1e' + str(int(baseS-1)) + ')')
        else:
            ax.set_zticklabels([ format_number(tmpi) for tmpi in np.arange(start, stop, step)] )
            ax.set_zlabel(str(plot_name) + ' expression')
    
    else:
        ax.set_zlabel(str(plot_name) + ' expression')
    
    ax.set_box_aspect((1, 1 * (max(ax.get_ylim3d()) - min(ax.get_ylim3d())) / (max(ax.get_xlim3d()) - min(ax.get_xlim3d())), z_height))
    ax.set_xlabel('Spatial 1', labelpad=-6)
    ax.set_ylabel('Spatial 2', labelpad=-6)
    
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False, labelleft=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False, labelleft=False) # labels along the bottom edge are off

    import matplotlib.ticker as ticker
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

    plt.tight_layout()
    
    if out_path is None:
        pass
    else:
        plt.savefig(out_path + '/' + str(plot_name) + '_single.pdf', transparent=True, format = 'pdf', bbox_inches='tight')
    
    if return_ax == True:
        return ax
        plt.show()
    else:
        plt.show()


def pl_multipleOUT(adata=None, find_obs=False, use_raw=True, input_zA=None, input_zB=None, plot_nameA='Gene', plot_nameB='Gene', normalize_z=True, normlabels=True, out_path=None, obs_cluster=None, color_map=None, cluster=None, mesh_directory='./Ridge_tmp', save_address='', lw=0, elev=35, view=-70, z_height=1, bg_alphaA=1, bg_alphaB=1, color_brighten=True, plot_HE=False, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, HE_slice=1, plot_clustering=False, clustering_z=None, clustering_size=1, clustering_alpha=0.1, select_c=None, xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, truncA=None, truncB=None, figsize_set=(8,8), return_ax=False):
    
    import pickle
    import os

    if save_address != '':
        mesh_directory = mesh_directory + '/' + save_address + '/'
    else:
        mesh_directory = mesh_directory

    if save_address != '':
        HE_directory = HE_directory + '/' + save_address + '/'
    else:
        HE_directory = HE_directory

    if input_zA is None:
        files = os.listdir(mesh_directory)
        if str(plot_nameA) + "_denoise.pickle" in files:
            with open(str(mesh_directory) + "/" + str(plot_nameA) + "_denoise.pickle", "rb") as input_file:
                input_zA = pickle.load(input_file)
        else:
            print("No file named " + str(plot_nameA) + " in " + str(mesh_directory))
            if adata is not None:
                print('Using noisy information in adata')
                if find_obs == True:
                    input_zA = adata.obs[plot_nameA]
                else:
                    from scipy.sparse import csr_matrix
                    if use_raw == True:
                        matrix = csr_matrix(adata.raw.X)
                        var_names = adata.raw.var_names
                    else:
                        matrix = csr_matrix(adata.X)
                        var_names = adata.var_names

                    input_zA = matrix[:,var_names == plot_nameA].toarray().ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_zA = input_zA

    if input_zB is None:
        files = os.listdir(mesh_directory)
        if str(plot_nameB) + "_denoise.pickle" in files:
            with open(str(mesh_directory) + "/" + str(plot_nameB) + "_denoise.pickle", "rb") as input_file:
                input_zB = pickle.load(input_file)
        else:
            print("No file named " + str(plot_nameB) + " in " + str(mesh_directory))
            if adata is not None:
                print('Using noisy information in adata')
                if find_obs == True:
                    input_zB = adata.obs[plot_nameB]
                else:
                    from scipy.sparse import csr_matrix
                    if use_raw == True:
                        matrix = csr_matrix(adata.raw.X)
                        var_names = adata.raw.var_names
                    else:
                        matrix = csr_matrix(adata.X)
                        var_names = adata.var_names

                    input_zB = matrix[:,var_names == plot_nameB].toarray().ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_zB = input_zB

    import copy
    input_zPA = input_zA.copy()
    normalize_factorA = 1
    normalize_minA = 0
    if normalize_z == True:
        normalize_factorA = 100 / max(input_zA)
        normalize_minA = np.min(input_zA)
        norm_zA = normalize_factorA * (input_zA - normalize_minA)#z / np.linalg.norm(z)
    else:
        norm_zA = input_zA

    input_zPB = input_zB.copy()
    normalize_factorB = 1
    normalize_minB = 0
    if normalize_z == True:
        normalize_factorB = 100 / max(input_zB)
        normalize_minB = np.min(input_zB)
        norm_zB = normalize_factorB * (input_zB - normalize_minB) #z / np.linalg.norm(z)
    else:
        norm_zB = input_zB


    # if normalize_z == True:
    #     norm_zA = 100 * input_zA / max(input_zA)#z / np.linalg.norm(z)
    #     norm_zB = 100 * input_zB / max(input_zB)#z / np.linalg.norm(z)
    # else:
    #     norm_zA = input_zA
    #     norm_zB = input_zB

    ## trunc norm_z
    if truncA is not None:
        tmp = []
        for mv in norm_zA:
            if mv > truncA:
                tmp.append(mv-truncA)
            else:
                tmp.append(0)
        norm_zA = np.array(tmp)
    
    if truncB is not None:
        tmp = []
        for mv in norm_zB:
            if mv > truncB:
                tmp.append(mv-truncB)
            else:
                tmp.append(0)
        norm_zB = np.array(tmp)    


    ##read triang mesh
    import pickle
    if os.path.isfile(str(mesh_directory) + "/triang.pickle"):
        with open(str(mesh_directory) + "/triang.pickle", "rb") as input_file:
            max_r, triang = pickle.load(input_file)
    else:
        print('No triangular mesh detected. Please run RidgeSpace.tl_mesh() first.')
        return 0

    if adata is not None:
        cluster = adata.obs[obs_cluster].values
        uns_cluster = obs_cluster + '_colors'
        if uns_cluster in adata.uns.keys():
            cluster_colors = adata.uns[uns_cluster]
        else:
            cluster_colors = None

    if cluster_colors is None:
        if color_map is None:
            cluster_colors = ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#e377c2', '#8c564b', '#17becf', '#b5bd61', '#aec7e8', '#aa40fc', '#ffbb78', '#98df8a','#ff9896', '#11e6df', '#e9812d', '#e5cb83', '#b2c3d8', '#bee6da', '#cccccc', '#906bac', '#f51d1f', '#c49d98']
        else:
            cluster_colors = color_map

    try:
        tmp = cluster.categories
    except AttributeError:
        cluster = pd.Categorical(cluster)

    colorsA, original_colorsA = color_alpha(input_z=norm_zA, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alphaA, color_brighten=color_brighten)
    colorsB, original_colorsB = color_alpha(input_z=norm_zB, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alphaB, color_brighten=color_brighten)

    fig = plt.figure(figsize=figsize_set)
    try:
        ax = fig.add_subplot(1,1,1, projection='3d', computed_zorder=False)
    except:
        ax = fig.add_subplot(1,1,1, projection='3d')
        print('Using computed zorder!')

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['pdf.use14corefonts'] = True
    # plt.rcParams['font.family'] = 'sans-serif'
    # plt.rcParams['font.sans-serif'] = 'Arial'

    if select_c is not None:
        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colorsA])
        colorsA = colors_new

        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colorsB])
        colorsB = colors_new


    # plot_Ridge(ax, triang, norm_zA, lw = lw, facecolors = colorsA, antialiased=False, rasterized=True, zorder=3)
    plot_Ridge(ax, triang, -1 * norm_zB, lw = lw, facecolors = colorsB, antialiased=False, rasterized=True, zorder=0)
    
    if plot_HE == True:
        import pickle
        import os
        if HE_filter == False:
            if os.path.isfile(str(HE_directory) + "/HE.pickle"):
                with open(str(HE_directory) + "/HE.pickle", "rb") as input_file:
                    sc_x, sc_y, sc_z = pickle.load(input_file)
            else:
                print('Please provide hitology image and run RidgeSpace.tl_HE()')
                return 0
                
        elif HE_filter == True:
            if os.path.isfile(str(HE_directory) + "/HE_filter.pickle"):
                with open(str(HE_directory) + "/HE_filter.pickle", "rb") as input_file:
                    sc_x, sc_y, sc_z = pickle.load(input_file)
            else:
                print('Please provide hitology image and run RidgeSpace.tl_HE()')
                return 0

        HE_slice = int(HE_slice)
        if HE_slice < 1:
            HE_slice = 1
        if HE_z is None:
            HE_z = 0
            ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=1)
        else:
            ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=1)

        
    if plot_clustering == True:
        if clustering_z is None:
            clustering_z = 0
            ax.scatter3D(xs=triang.x, ys=triang.y, zs=clustering_z, s=clustering_size, facecolors=np.array(original_colorsA), alpha=clustering_alpha, antialiased=False, linewidth=0, rasterized=True, zorder=1)#
        else:
            ax.scatter3D(xs=triang.x, ys=triang.y, zs=clustering_z, s=clustering_size, facecolors=np.array(original_colorsA), alpha=clustering_alpha, antialiased=False, linewidth=0, rasterized=True, zorder=1)#

    plot_Ridge(ax, triang, norm_zA, lw = lw, facecolors = colorsA, antialiased=False, rasterized=True, zorder=2)
    ax.view_init(elev=elev, azim=view)
    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_zlim(zlim_min, zlim_max)
    
    X1 = [min(ax.get_xlim3d()), min(ax.get_xlim3d()), max(ax.get_xlim3d()), max(ax.get_xlim3d()), min(ax.get_xlim3d())]
    Y1 = [min(ax.get_ylim3d()), max(ax.get_ylim3d()), max(ax.get_ylim3d()), min(ax.get_ylim3d()), min(ax.get_ylim3d())]
    for i in range(4):
        j = i + 1
        if j == 4:
            j = 0
        ax.plot([X1[i], X1[j]], [Y1[i], Y1[j]], [-0.5, -0.5], color='gray', lw=1, alpha=1, zorder = 1)#

    if normlabels == False:
        if truncA is not None:
            # input_zPA norm_zA normalize_minA normalize_factorA
            input_z_minA, input_z_maxA = max(input_zPA) * truncA / 100, max(input_zPA)
            startA, stopA, stepA = get_plotlabel_trunc(0.95 * input_z_minA, 1.05 * input_z_maxA)
            plot_labelsA = []
            for labels_iA in np.arange(startA, stopA, stepA):
                plot_labelsA.append((labels_iA - normalize_minA) * normalize_factorA - truncA)
            
        else:
            # input_zPA norm_zA normalize_minA normalize_factorA
            input_z_minA, input_z_maxA = min(input_zPA), max(input_zPA)
            startA, stopA, stepA = get_plotlabel(0.95 * input_z_minA, 1.05 * input_z_maxA)
            plot_labelsA = []
            for labels_iA in np.arange(startA, stopA, stepA):
                plot_labelsA.append((labels_iA - normalize_minA) * normalize_factorA)

        if truncB is not None:
            # input_zPB norm_zB normalize_minB normalize_factorB
            input_z_minB, input_z_maxB = max(input_zPB) * truncB / 100, max(input_zPB)
            startB, stopB, stepB = get_plotlabel_trunc(0.95 * input_z_minB, 1.05 * input_z_maxB)
            plot_labelsB = []
            for labels_iB in np.arange(startB, stopB, stepB):
                plot_labelsB.append((labels_iB - normalize_minB) * normalize_factorB - truncB)
            
        else:
            # input_zPB norm_zB normalize_minB normalize_factorB
            input_z_minB, input_z_maxB = min(input_zPB), max(input_zPB)
            startB, stopB, stepB = get_plotlabel(0.95 * input_z_minB, 1.05 * input_z_maxB)
            plot_labelsB = []
            for labels_iB in np.arange(startB, stopB, stepB):
                plot_labelsB.append((labels_iB - normalize_minB) * normalize_factorB)
        
        plot_labels = []
        for plot_labelsBi in plot_labelsB:
            if plot_labelsBi >= 0:
                plot_labels.append(-1*plot_labelsBi)
            else:
                plot_labels.append(plot_labelsBi)

        for plot_labelsAi in plot_labelsA:
            if plot_labelsAi >= 0:
                plot_labels.append(plot_labelsAi)
            else:
                plot_labels.append(-1 * plot_labelsAi)

        ax.set_zticks(plot_labels)
        if stepA > 50:
            baseSA = np.round(np.log10(stepA))
            baseSNA = 10 ** (baseSA - 1)
            set_zticklabelsA = np.arange(startA / baseSNA, stopA / baseSNA, stepA / baseSNA)
            set_zlabelA = str(plot_nameA) + ' (1e' + str(int(baseSA - 1)) + ')'
        else:
            set_zticklabelsA = np.arange(startA, stopA, stepA)
            set_zlabelA = str(plot_nameA)

        if stepB > 50:
            baseSB = np.round(np.log10(stepB))
            baseSNB = 10 ** (baseSB - 1)
            set_zticklabelsB = np.arange(startB / baseSNB, stopB / baseSNB, stepB / baseSNB)
            set_zlabelB = str(plot_nameB) + ' (1e' + str(int(baseSB-1)) + ')'
        else:
            set_zticklabelsB = np.arange(startB, stopB, stepB)
            set_zlabelB = str(plot_nameB)

        ax.set_zticklabels([ format_number(tmpi) for tmpi in np.concatenate([set_zticklabelsB, set_zticklabelsA]) ])
        ax.set_zlabel(set_zlabelB + ' & ' + set_zlabelA + ' expression')
    
    else:
        # ticks = ax.get_zticklabels()
        # ax.set_zticklabels( [format_number(np.abs(float(tmpi))) for tmpi in ticks])
        ax.set_zlabel(str(plot_nameB) + ' & '+ str(plot_nameA) + ' expression')

    ax.set_box_aspect((1, 1 * (max(ax.get_ylim3d()) - min(ax.get_ylim3d())) / (max(ax.get_xlim3d()) - min(ax.get_xlim3d())), z_height))
    ax.set_xlabel('Spatial 1', labelpad=-6)
    ax.set_ylabel('Spatial 2', labelpad=-6)
    

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False, labelleft=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False, labelleft=False) # labels along the bottom edge are off

    import matplotlib.ticker as ticker
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

    plt.tight_layout()
    if out_path is None:
        pass
    else:
        plt.savefig(out_path + '/' + str(plot_nameA) + '_' + str(plot_nameB) + '_multipleOUT.pdf', transparent=True, format = 'pdf', bbox_inches='tight')
    
    if return_ax == True:
        return ax
        plt.show()
    else:
        plt.show()


def pl_multipleIN(adata=None, find_obs=False, use_raw=True, input_zA=None, input_zB=None, plot_nameA='Gene', plot_nameB='Gene', height=100, normalize_z=True, normlabels=False, out_path=None, obs_cluster=None, color_map=None, cluster=None, mesh_directory='./Ridge_tmp', save_address='', lw=0, elev=35, view=-70, z_height=1, bg_alphaA=1, bg_alphaB=1, color_brighten=True, plot_HE=False, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, HE_slice=1, plot_clustering=False, clustering_z=None, clustering_size=1, clustering_alpha=0.1, plot_clustering_boundary=False, clustering_boundary_lw=1, clustering_boundary_alpha=0.4, select_c=None,xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, truncA=None, truncB=None, figsize_set=(8,8), return_ax=False):

    import pickle
    import os
    if save_address != '':
        mesh_directory = mesh_directory + '/' + save_address + '/'
    else:
        mesh_directory = mesh_directory

    if save_address != '':
        HE_directory = HE_directory + '/' + save_address + '/'
    else:
        HE_directory = HE_directory

    if input_zA is None:
        files = os.listdir(mesh_directory)
        if str(plot_nameA) + "_denoise.pickle" in files:
            with open(str(mesh_directory) + "/" + str(plot_nameA) + "_denoise.pickle", "rb") as input_file:
                input_zA = pickle.load(input_file)
        else:
            print("No file named " + str(plot_nameA) + " in " + str(mesh_directory))
            if adata is not None:
                print('Using noisy information in adata')
                if find_obs == True:
                    input_zA = adata.obs[plot_nameA]
                else:
                    from scipy.sparse import csr_matrix
                    if use_raw == True:
                        matrix = csr_matrix(adata.raw.X)
                        var_names = adata.raw.var_names
                    else:
                        matrix = csr_matrix(adata.X)
                        var_names = adata.var_names

                    input_zA = matrix[:,var_names == plot_nameA].toarray().ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_zA = input_zA

    if input_zB is None:
        files = os.listdir(mesh_directory)
        if str(plot_nameB) + "_denoise.pickle" in files:
            with open(str(mesh_directory) + "/" + str(plot_nameB) + "_denoise.pickle", "rb") as input_file:
                input_zB = pickle.load(input_file)
        else:
            print("No file named " + str(plot_nameB) + " in " + str(mesh_directory))
            if adata is not None:
                print('Using noisy information in adata')
                if find_obs == True:
                    input_zB = adata.obs[plot_nameB]
                else:
                    from scipy.sparse import csr_matrix
                    if use_raw == True:
                        matrix = csr_matrix(adata.raw.X)
                        var_names = adata.raw.var_names
                    else:
                        matrix = csr_matrix(adata.X)
                        var_names = adata.var_names

                    input_zB = matrix[:,var_names == plot_nameB].toarray().ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_zB = input_zB

    import copy
    input_zPA = input_zA.copy()
    normalize_factorA = 1
    normalize_minA = 0
    if normalize_z == True:
        normalize_factorA = 100 / max(input_zA)
        normalize_minA = np.min(input_zA)
        norm_zA = normalize_factorA * (input_zA - normalize_minA)#z / np.linalg.norm(z)
    else:
        norm_zA = input_zA

    input_zPB = input_zB.copy()
    normalize_factorB = 1
    normalize_minB = 0
    if normalize_z == True:
        normalize_factorB = 100 / max(input_zB)
        normalize_minB = np.min(input_zB)
        norm_zB = normalize_factorB * (input_zB - normalize_minB) #z / np.linalg.norm(z)
    else:
        norm_zB = input_zB

    # if normalize_z == True:
    #     norm_zA = 100 * input_zA / max(input_zA)#z / np.linalg.norm(z)
    #     norm_zB = 100 * input_zB / max(input_zB)#z / np.linalg.norm(z)
    # else:
    #     norm_zA = input_zA
    #     norm_zB = input_zB

    ## trunc norm_z
    if truncA is not None:
        tmp = []
        for mv in norm_zA:
            if mv > truncA:
                tmp.append(mv-truncA)
            else:
                tmp.append(0)
        norm_zA = np.array(tmp)
    
    if truncB is not None:
        tmp = []
        for mv in norm_zB:
            if mv > truncB:
                tmp.append(mv-truncB)
            else:
                tmp.append(0)
        norm_zB = np.array(tmp)    


    ##read triang mesh
    import pickle
    if os.path.isfile(str(mesh_directory) + "/triang.pickle"):
        with open(str(mesh_directory) + "/triang.pickle", "rb") as input_file:
            max_r, triang = pickle.load(input_file)
    else:
        print('No triangular mesh detected. Please run RidgeSpace.tl_mesh() first.')
        return 0

    if adata is not None:
        cluster = adata.obs[obs_cluster].values
        uns_cluster = obs_cluster + '_colors'
        if uns_cluster in adata.uns.keys():
            cluster_colors = adata.uns[uns_cluster]
        else:
            cluster_colors = None

    if cluster_colors is None:
        if color_map is None:
            cluster_colors = ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#e377c2', '#8c564b', '#17becf', '#b5bd61', '#aec7e8', '#aa40fc', '#ffbb78', '#98df8a','#ff9896', '#11e6df', '#e9812d', '#e5cb83', '#b2c3d8', '#bee6da', '#cccccc', '#906bac', '#f51d1f', '#c49d98']
        else:
            cluster_colors = color_map

    try:
        tmp = cluster.categories
    except AttributeError:
        cluster = pd.Categorical(cluster)

    colorsA, original_colorsA = color_alpha(input_z=norm_zA, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alphaA, color_brighten=color_brighten)
    colorsB, original_colorsB = color_alpha(input_z=norm_zB, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alphaB, color_brighten=color_brighten)

    fig = plt.figure(figsize=figsize_set)
    try:
        ax = fig.add_subplot(1,1,1, projection='3d', computed_zorder=False)
    except:
        ax = fig.add_subplot(1,1,1, projection='3d')
        print('Using computed zorder!')

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['pdf.use14corefonts'] = True
    # plt.rcParams['font.family'] = 'sans-serif'
    # plt.rcParams['font.sans-serif'] = 'Arial'

    if select_c is not None:
        print(np.unique([i[:-2] for i in colorsA]))
        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colorsA])
        colorsA = colors_new
        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colorsB])
        colorsB = colors_new

    plot_Ridge(ax, triang, height - norm_zA, lw = lw, facecolors = colorsA, antialiased=False, rasterized=True, zorder=2)
    
    
    if plot_HE == True:
        import pickle
        import os
        if HE_filter == False:
            if os.path.isfile(str(HE_directory) + "/HE.pickle"):
                with open(str(HE_directory) + "/HE.pickle", "rb") as input_file:
                    sc_x, sc_y, sc_z = pickle.load(input_file)
            else:
                print('Please provide hitology image and run RidgeSpace.tl_HE()')
                return 0
                
        elif HE_filter == True:
            if os.path.isfile(str(HE_directory) + "/HE_filter.pickle"):
                with open(str(HE_directory) + "/HE_filter.pickle", "rb") as input_file:
                    sc_x, sc_y, sc_z = pickle.load(input_file)
            else:
                print('Please provide hitology image and run RidgeSpace.tl_HE()')
                return 0

        HE_slice = int(HE_slice)
        if HE_slice < 1:
            HE_slice = 1                
        if HE_z is None:
            HE_z = height
            ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=3)
            HE_z = 0
            ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=0)
        elif type(HE_z) == list:
            ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=min(HE_z), facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=0)
            ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=max(HE_z), facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=3)
        else:
            if height - clustering_z >  clustering_z:
                ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=0)
            else:
                ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=3)

    if plot_clustering == True:
        if clustering_z is None:
            clustering_z = height
            ax.scatter3D(xs=triang.x, ys=triang.y, zs=clustering_z, s=clustering_size, facecolors=np.array(original_colorsA), alpha=clustering_alpha, antialiased=False, linewidth=0, rasterized=True, zorder=3)#
        elif type(clustering_z) == list:
            ax.scatter3D(xs=triang.x, ys=triang.y, zs=min(clustering_z), s=clustering_size, facecolors=np.array(original_colorsA), alpha=clustering_alpha, antialiased=False, linewidth=0, rasterized=True, zorder=0)
            ax.scatter3D(xs=triang.x, ys=triang.y, zs=max(clustering_z), s=clustering_size, facecolors=np.array(original_colorsA), alpha=clustering_alpha, antialiased=False, linewidth=0, rasterized=True, zorder=3)
        else:
            if height - clustering_z >  clustering_z:
                ax.scatter3D(xs=triang.x, ys=triang.y, zs=clustering_z, s=clustering_size, facecolors=np.array(original_colorsA), alpha=clustering_alpha, antialiased=False, linewidth=0, rasterized=True, zorder=0)#
            else:
                ax.scatter3D(xs=triang.x, ys=triang.y, zs=clustering_z, s=clustering_size, facecolors=np.array(original_colorsA), alpha=clustering_alpha, antialiased=False, linewidth=0, rasterized=True, zorder=3)#

    if plot_clustering_boundary == True:

        for colori in pd.Series(original_colorsA).unique():
            Tx = triang.x[(original_colorsA == colori)]
            Ty = triang.y[(original_colorsA == colori)]
            T = mtri.Triangulation(Tx, Ty)#
            T.is_delaunay=False
            # plot only triangles with sidelength smaller some max_radius
            Ttriangles = T.triangles

            # Mask off unwanted triangles.
            Txtri = Tx[Ttriangles] - np.roll(Tx[Ttriangles], 1, axis=1)
            Tytri = Ty[Ttriangles] - np.roll(Ty[Ttriangles], 1, axis=1)
            Tmaxi = np.max(np.sqrt(Txtri**2 + Tytri**2), axis=1)

            T.set_mask(Tmaxi > max_r)


            # from Stack Overflow
            boundaries = []
            for i in range(len(T.triangles)):
                if T.mask[i] == False:
                    for j in range(3):
                        if T.neighbors[i,j] < 0:
                            # Triangle edge (i,j) has no neighbor so is a boundary edge.
                            boundaries.append((T.triangles[i,j],
                                            T.triangles[i,(j+1)%3]))
            boundaries = np.asarray(boundaries)

            # The following line plots the boundary edges.
            bheight = height
            for i in range(len(boundaries)):
                ax.plot([Tx[boundaries].T[0,i], Tx[boundaries].T[1,i]], [Ty[boundaries].T[0,i], Ty[boundaries].T[1,i]], \
                        [bheight, bheight], color=colori, lw=clustering_boundary_lw, alpha=clustering_boundary_alpha, zorder=3)#
                        
    plot_Ridge(ax, triang, norm_zB, lw = lw, facecolors = colorsB, antialiased=False, rasterized=True, zorder=1)
    ax.view_init(elev=elev, azim=view)
    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_zlim(zlim_min, zlim_max)
    
    if normlabels == False:
        if truncA is not None:
            # input_zPA norm_zA normalize_minA normalize_factorA
            input_z_minA, input_z_maxA = max(input_zPA) * truncA / 100, max(input_zPA)
            startA, stopA, stepA = get_plotlabel_trunc(0.95 * input_z_minA, 1.05 * input_z_maxA)
            plot_labelsA = []
            for labels_iA in np.arange(startA, stopA, stepA):
                plot_labelsA.append((labels_iA - normalize_minA) * normalize_factorA - truncA)
            
        else:
            # input_zPA norm_zA normalize_minA normalize_factorA
            input_z_minA, input_z_maxA = min(input_zPA), max(input_zPA)
            startA, stopA, stepA = get_plotlabel(0.95 * input_z_minA, 1.05 * input_z_maxA)
            plot_labelsA = []
            for labels_iA in np.arange(startA, stopA, stepA):
                plot_labelsA.append((labels_iA - normalize_minA) * normalize_factorA)

        if truncB is not None:
            # input_zPB norm_zB normalize_minB normalize_factorB
            input_z_minB, input_z_maxB = max(input_zPB) * truncB / 100, max(input_zPB)
            startB, stopB, stepB = get_plotlabel_trunc(0.95 * input_z_minB, 1.05 * input_z_maxB)
            plot_labelsB = []
            for labels_iB in np.arange(startB, stopB, stepB):
                plot_labelsB.append((labels_iB - normalize_minB) * normalize_factorB - truncB)
            
        else:
            # input_zPB norm_zB normalize_minB normalize_factorB
            input_z_minB, input_z_maxB = min(input_zPB), max(input_zPB)
            startB, stopB, stepB = get_plotlabel(0.95 * input_z_minB, 1.05 * input_z_maxB)
            plot_labelsB = []
            for labels_iB in np.arange(startB, stopB, stepB):
                plot_labelsB.append((labels_iB - normalize_minB) * normalize_factorB)
        
        plot_labelsBF = []
        for plot_labelsBi in plot_labelsB:
            if plot_labelsBi >= 0:
                plot_labelsBF.append(height - plot_labelsBi)
            else:
                plot_labelsBF.append(height + plot_labelsBi)

        plot_labelsAF = []
        for plot_labelsAi in plot_labelsA:
            if plot_labelsAi >= 0:
                plot_labelsAF.append(plot_labelsAi)
            else:
                plot_labelsAF.append(-1 * plot_labelsAi)

        if stepA > 50:
            baseSA = np.round(np.log10(stepA))
            baseSNA = 10 ** (baseSA - 1)
            set_zticklabelsA = np.arange(startA / baseSNA, stopA / baseSNA, stepA / baseSNA)
            set_zlabelA = str(plot_nameA) + ' (1e' + str(int(baseSA - 1)) + ')'
        else:
            set_zticklabelsA = np.arange(startA, stopA, stepA)
            set_zlabelA = str(plot_nameA)

        if stepB > 50:
            baseSB = np.round(np.log10(stepB))
            baseSNB = 10 ** (baseSB - 1)
            set_zticklabelsB = np.arange(startB / baseSNB, stopB / baseSNB, stepB / baseSNB)
            set_zlabelB = str(plot_nameB) + ' (1e' + str(int(baseSB-1)) + ')'
        else:
            set_zticklabelsB = np.arange(startB, stopB, stepB)
            set_zlabelB = str(plot_nameB)

        ax1 = fig.add_axes(MyAxes3D(ax, 'lr'))
        ax1.set_custom_ticks('l', ticks=plot_labelsBF, ticklabels=[format_number(tmpi) for tmpi in set_zticklabelsB])
        ax1.set_zaxis_label('l', set_zlabelB + ' expression')
        ax1.set_custom_ticks('r', ticks=plot_labelsAF, ticklabels=[format_number(tmpi) for tmpi in set_zticklabelsA])
        ax1.set_zaxis_label('r', set_zlabelA + ' expression')
        # ax.set_zticklabels([ format_number(tmpi) for tmpi in np.concatenate([set_zticklabelsB, set_zticklabelsA]) ])
        # ax.set_zlabel(set_zlabelB + ' & ' + set_zlabelA + ' expression')
    
    else:
        
        if truncA is not None:
            # input_zPA norm_zA normalize_minA normalize_factorA
            input_z_minA, input_z_maxA = max(norm_zA) * truncA / 100, max(norm_zA)
            startA, stopA, stepA = get_plotlabel_trunc(0.95 * input_z_minA, 1.05 * input_z_maxA)
            plot_labelsA = np.arange(startA, stopA, stepA)
            
        else:
            # input_zPA norm_zA normalize_minA normalize_factorA
            input_z_minA, input_z_maxA = min(norm_zA), max(norm_zA)
            startA, stopA, stepA = get_plotlabel(0.95 * input_z_minA, 1.05 * input_z_maxA)
            plot_labelsA = np.arange(startA, stopA, stepA)

        if truncB is not None:
            # input_zPB norm_zB normalize_minB normalize_factorB
            input_z_minB, input_z_maxB = max(norm_zB) * truncB / 100, max(norm_zB)
            startB, stopB, stepB = get_plotlabel_trunc(0.95 * input_z_minB, 1.05 * input_z_maxB)
            plot_labelsB = np.arange(startB, stopB, stepB)
            
        else:
            # input_zPB norm_zB normalize_minB normalize_factorB
            input_z_minB, input_z_maxB = min(norm_zB), max(norm_zB)
            startB, stopB, stepB = get_plotlabel(0.95 * input_z_minB, 1.05 * input_z_maxB)
            plot_labelsB = np.arange(startB, stopB, stepB)

        plot_labelsBF = []
        for plot_labelsBi in plot_labelsB:
            if plot_labelsBi >= 0:
                plot_labelsBF.append(height - plot_labelsBi)
            else:
                plot_labelsBF.append(height + plot_labelsBi)

        plot_labelsAF = []
        for plot_labelsAi in plot_labelsA:
            if plot_labelsAi >= 0:
                plot_labelsAF.append(plot_labelsAi)
            else:
                plot_labelsAF.append(-1 * plot_labelsAi)

        if stepA > 50:
            baseSA = np.round(np.log10(stepA))
            baseSNA = 10 ** (baseSA - 1)
            set_zticklabelsA = np.arange(startA / baseSNA, stopA / baseSNA, stepA / baseSNA)
            set_zlabelA = str(plot_nameA) + ' (1e' + str(int(baseSA - 1)) + ')'
        else:
            set_zticklabelsA = np.arange(startA, stopA, stepA)
            set_zlabelA = str(plot_nameA)

        if stepB > 50:
            baseSB = np.round(np.log10(stepB))
            baseSNB = 10 ** (baseSB - 1)
            set_zticklabelsB = np.arange(startB / baseSNB, stopB / baseSNB, stepB / baseSNB)
            set_zlabelB = str(plot_nameB) + ' (1e' + str(int(baseSB-1)) + ')'
        else:
            set_zticklabelsB = np.arange(startB, stopB, stepB)
            set_zlabelB = str(plot_nameB)

        ax1 = fig.add_axes(MyAxes3D(ax, 'lr'))
        ax1.set_custom_ticks('l', ticks=plot_labelsBF, ticklabels=[format_number(tmpi) for tmpi in set_zticklabelsB])
        ax1.set_zaxis_label('l', set_zlabelB + ' expression')
        ax1.set_custom_ticks('r', ticks=plot_labelsAF, ticklabels=[format_number(tmpi) for tmpi in set_zticklabelsA])
        ax1.set_zaxis_label('r', set_zlabelA + ' expression')

        # ticks = ax.get_zticklabels()
        # ax.set_zticklabels( [format_number(np.abs(float(tmpi))) for tmpi in ticks])
        # ax.set_zlabel(str(plot_nameB) + ' & '+ str(plot_nameA) + ' expression')

    ax.set_box_aspect((1, 1 * (max(ax.get_ylim3d()) - min(ax.get_ylim3d())) / (max(ax.get_xlim3d()) - min(ax.get_xlim3d())), z_height))
    ax.set_xlabel('Spatial 1', labelpad=-6)
    ax.set_ylabel('Spatial 2', labelpad=-6)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False, labelleft=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False, labelleft=False) # labels along the bottom edge are off

    import matplotlib.ticker as ticker
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

    plt.tight_layout()
    if out_path is None:
        pass
    else:
        plt.savefig(out_path + '/' + str(plot_nameA) + '_' + str(plot_nameB) + '_multipleOUT.pdf', transparent=True, format = 'pdf', bbox_inches='tight')
    
    if return_ax == True:
        return ax
        plt.show()
    else:
        plt.show()


def pl_trajectory(adata=None, input_z=None, find_obs=True, use_raw=True, normalize_z=False, normlabels=True, out_path=None, obs_cluster=None, color_map=None, cluster=None, mesh_directory='./Ridge_tmp', save_address='', lw=0, height_ratio=1.5, elev=35, view=-70, z_height=1, bg_alpha=1, color_brighten=True, plot_name='Gene', plot_HE=True, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, HE_slice=1, select_c=None, xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, trunc=None, figsize_set=(8,8), density_arrow = 0.4, dense=  True, min_length = 0.15, color_arrow="k", lw_arrow=1, arrow_scale=20, return_ax=False):

    import pickle
    import os

    if save_address != '':
        mesh_directory = mesh_directory + '/' + save_address + '/'
    else:
        mesh_directory = mesh_directory

    if save_address != '':
        HE_directory = HE_directory + '/' + save_address + '/'
    else:
        HE_directory = HE_directory

    if input_z is None:
        files = os.listdir(mesh_directory)
        if str(plot_name) + "_denoise.pickle" in files:
            with open(str(mesh_directory) + "/" + str(plot_name) + "_denoise.pickle", "rb") as input_file:
                input_z = pickle.load(input_file)
        else:
            print("No file named " + str(plot_name) + " in " + str(mesh_directory))
            if adata is not None:
                print('Using noisy information in adata')
                if find_obs == True:
                    input_z = adata.obs[plot_name]
                else:
                    from scipy.sparse import csr_matrix
                    if use_raw == True:
                        matrix = csr_matrix(adata.raw.X)
                        var_names = adata.raw.var_names
                    else:
                        matrix = csr_matrix(adata.X)
                        var_names = adata.var_names

                    input_z = matrix[:,var_names == plot_name].toarray().ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_z = input_z

    import copy
    input_zP = input_z.copy()
    normalize_factor = 1
    normalize_min = 0
    if normalize_z == True:
        normalize_factor = 100 / max(input_z)
        normalize_min = np.min(input_z)
        norm_z = normalize_factor * (input_z - normalize_min) #z / np.linalg.norm(z)
    else:
        norm_z = input_z


    if trunc is not None:
        tmp = []
        for mv in norm_z:
            if mv > trunc:
                tmp.append(mv-trunc)
            else:
                tmp.append(0)
        norm_z = np.array(tmp)

    ##read triang mesh
    import pickle
    if os.path.isfile(str(mesh_directory) + "/triang.pickle"):
        with open(str(mesh_directory) + "/triang.pickle", "rb") as input_file:
            max_r, triang = pickle.load(input_file)
    else:
        print('No triangular mesh detected. Please run RidgeSpace.tl_mesh() first.')
        return 0

    if adata is not None:
        cluster = adata.obs[obs_cluster].values
        uns_cluster = obs_cluster + '_colors'
        if uns_cluster in adata.uns.keys():
            cluster_colors = adata.uns[uns_cluster]
        else:
            cluster_colors = None

    if cluster_colors is None:
        if color_map is None:
            cluster_colors = ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#e377c2', '#8c564b', '#17becf', '#b5bd61', '#aec7e8', '#aa40fc', '#ffbb78', '#98df8a','#ff9896', '#11e6df', '#e9812d', '#e5cb83', '#b2c3d8', '#bee6da', '#cccccc', '#906bac', '#f51d1f', '#c49d98']
        else:
            cluster_colors = color_map

    try:
        tmp = cluster.categories
    except AttributeError:
        cluster = pd.Categorical(cluster)

    colors, original_colors = color_alpha(input_z=norm_z, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alpha, color_brighten=color_brighten)

    fig = plt.figure(figsize=figsize_set)
    try:
        ax = fig.add_subplot(1,1,1, projection='3d', computed_zorder=False)
    except:
        ax = fig.add_subplot(1,1,1, projection='3d')
        print('Using computed zorder!')

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['pdf.use14corefonts'] = True
    # plt.rcParams['font.family'] = 'sans-serif'
    # plt.rcParams['font.sans-serif'] = 'Arial'

    if select_c is not None:
        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colors])
        colors = colors_new

    plot_Ridge(ax, triang, norm_z, lw = lw, facecolors = colors, antialiased=False, rasterized=True, zorder=1)

    if plot_HE == True:
        import pickle
        import os
        if HE_filter == False:
            if os.path.isfile(str(HE_directory) + "/HE.pickle"):
                with open(str(HE_directory) + "/HE.pickle", "rb") as input_file:
                    sc_x, sc_y, sc_z = pickle.load(input_file)
            else:
                print('Please provide hitology image and run RidgeSpace.tl_HE()')
                return 0
                
        elif HE_filter == True:
            if os.path.isfile(str(HE_directory) + "/HE_filter.pickle"):
                with open(str(HE_directory) + "/HE_filter.pickle", "rb") as input_file:
                    sc_x, sc_y, sc_z = pickle.load(input_file)
            else:
                print('Please provide hitology image and run RidgeSpace.tl_HE()')
                return 0

        HE_slice = int(HE_slice)
        if HE_slice < 1:
            HE_slice = 1                
        if HE_z is None:
            HE_z = 0
            ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=0)
        else:
            ax.scatter3D(sc_x[::HE_slice], sc_y[::HE_slice], zs=HE_z, facecolors=sc_z[::HE_slice], s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=0)

    plot_p, plot_lc = infer_trajectory(triang, norm_z, mesh_directory=mesh_directory, density = density_arrow, dense = dense, min_length = min_length, arrow_scale = arrow_scale, lw = lw_arrow, color = color_arrow)
    for pi in plot_p:
        ax.add_artist(pi)
    
    ax.add_artist(plot_lc)#add_collection

    ax.view_init(elev=elev, azim=view)

    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_zlim(zlim_min, zlim_max)

    if normlabels == False:
        if trunc is not None:
            # input_zP norm_z normalize_min normalize_factor
            input_z_min, input_z_max = max(input_zP) * trunc / 100, max(input_zP)
            start, stop, step = get_plotlabel_trunc(0.95 * input_z_min, 1.05 * input_z_max)
            plot_labels = []
            for labels_i in np.arange(start, stop, step):
                plot_labels.append((labels_i - normalize_min) * normalize_factor - trunc)
            
        else:
            # input_zP norm_z normalize_min normalize_factor
            # plot_labels = ax.get_ztickslabel()
            input_z_min, input_z_max = min(input_zP), max(input_zP)
            start, stop, step = get_plotlabel(0.95 * input_z_min, 1.05 * input_z_max)
            plot_labels = []
            for labels_i in np.arange(start, stop, step):
                plot_labels.append((labels_i - normalize_min) * normalize_factor)
        
        ax.set_zticks(plot_labels)
        if step > 50:
            baseS = np.round(np.log10(step))
            baseSN = 10 ** (baseS - 1)
            ax.set_zticklabels([ format_number(tmpi) for tmpi in np.arange(start / baseSN, stop / baseSN, step / baseSN)] )
            ax.set_zlabel(str(plot_name) + ' trajectory (1e' + str(int(baseS-1)) + ')')
        else:
            ax.set_zticklabels([ format_number(tmpi) for tmpi in np.arange(start, stop, step)] )
            ax.set_zlabel(str(plot_name) + ' trajectory')
    
    else:
        ax.set_zlabel(str(plot_name) + ' trajectory')

    ax.set_box_aspect((1, 1 * (max(ax.get_ylim3d()) - min(ax.get_ylim3d())) / (max(ax.get_xlim3d()) - min(ax.get_xlim3d())), z_height))#
    ax.set_xlabel('Spatial 1', labelpad=-6)
    ax.set_ylabel('Spatial 2', labelpad=-6)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False, labelleft=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False, labelleft=False) # labels along the bottom edge are off

    import matplotlib.ticker as ticker
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

    plt.tight_layout()
    if out_path is None:
        pass
    else:
        plt.savefig(out_path + '/' + str(plot_name) + '_trajectory.pdf', transparent=True, format = 'pdf', bbox_inches='tight')
    
    if return_ax == True:
        return ax
        plt.show()
    else:
        plt.show()


def convert_zarr(file_path):
    import numpy as np

    zarr_data = read_zarr(file_path)
    table_keys = [i for i in zarr_data.keys() if 'table' in i]
    if len(table_keys) >= 1:
        table_key = table_keys[0]
        if 'table' in zarr_data[table_key].keys():
            if 'X' in zarr_data[table_key]['table'].keys():
                Xe = zarr_data[table_key]['table']['X']
                from scipy.sparse import csr_matrix
                try:
                    expression_matrix = csr_matrix((Xe["data"], Xe["indices"], Xe["indptr"]))
                except Exception as e:
                    expression_matrix = csr_matrix(Xe)
                    
                obs = zarr_data[table_key]['table'].get('obs', {})
                var = zarr_data[table_key]['table'].get('var', {})
                uns = zarr_data[table_key]['table']['uns']
                spatial_coords = zarr_data[table_key]['table']['obsm'].get('spatial', None)

            else:
                raise KeyError("'X' (expression matrix) is missing in the zarr data.")
        else:
            if 'X' in zarr_data[table_key].keys():
                Xe = zarr_data[table_key]['X']
                from scipy.sparse import csr_matrix
                try:
                    expression_matrix = csr_matrix((Xe["data"], Xe["indices"], Xe["indptr"]))
                except Exception as e:
                    raise KeyError('Please check zarr expression matrix index.')

                obs = zarr_data[table_key].get('obs', {})
                var = zarr_data[table_key].get('var', {})
                uns = zarr_data[table_key]['uns']
                spatial_coords = zarr_data[table_key]['obsm'].get('spatial', None)       
            else:
                raise KeyError("'X' (expression matrix) is missing in the zarr data.")
    else:
        raise KeyError("Table is missing in the zarr data.")
        
    images = zarr_data.get('images', None)
    try:
        lowres = images[[i for i in images.keys() if 'lowres' in i][0]]
    except Exception as e:
        lowres = images[[i for i in zarr_data.get('images', None).keys()][0]]

    try:
        import anndata

    except Exception as e:
        raise ImportError('Please pip install anndata.')

    adata = anndata.AnnData(X=expression_matrix, obs=obs, var=var, uns=uns)

    if images is not None:
        if 'spatial' not in adata.uns.keys():
            adata.uns['spatial'] = dict()
        adata.uns['spatial']['library_id'] = dict()
        adata.uns['spatial']['library_id']['images'] = dict()
        adata.uns['spatial']['library_id']['images']['lowres'] = lowres['0']

    if spatial_coords is not None:
        adata.obsm['spatial'] = spatial_coords

    adata.var.index = adata.var['_index']
    adata.obs.index = adata.obs['_index']
    return adata