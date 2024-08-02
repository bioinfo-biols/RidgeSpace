from .utils import *
import os
import sys

def RidgeSpace_command():
    ''' Example of taking inputs for SEVtras'''
    args = sys.argv[1:]
    if len(args) < 1:
        print("usage: tl_mesh, tl_HE, tl_denoise, pl_single, pl_multipleOUT, pl_multipleIN, pl_trajectory")


def tl_mesh(adata=None, input_x=None, input_y=None, max_radius=None, plot_kde=False, save_directory='./Ridge_tmp', mesh_optimization=False, degree=2):

    if adata is not None:
        x = adata.obsm['spatial'][:,0]
        y = adata.obsm['spatial'][:,1]
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
    import os
    if not os.path.exists(str(save_directory)):
        os.makedirs(str(save_directory))

    with open(str(save_directory) + "/triang.pickle", "wb") as output_file:
        pickle.dump([max_radius1, triang], output_file)

    print('Triangles mesh saved in ' + str(save_directory))

    if mesh_optimization == True:
        mesh_optimize(triang=triang, save_directory=save_directory, degree = degree)



def tl_HE(adata=None, library_id=None, image_mask=None, save_directory='./Ridge_tmp', image_key='lowres', scalefactor='tissue_lowres_scalef', vscale=255, lightness=1.2, filter_spot=None):
    if adata is None:
        print('Please provide anndata with histology images in adata.uns!')
        return 0
    
    if library_id is None:
        library_id = list(adata.uns['spatial'].keys())[0]
        print('Using ' + str(library_id))

    lowres = adata.uns["spatial"][library_id]['images'][image_key]

    if image_mask is None:
        try:
            from skimage import data,filters
            import cv2
            image_mask, contours =get_HE_mask((lowres*vscale).astype('uint8'))
        except ImportError:
            raise ImportError("This method requires the `skimage` and `cv2` module to be installed. Or you can provide the binary mask of corresponding hitology image with parameter `image_mask`.")
        
    sc_x = []
    sc_y = []
    sc_z = []
    for i in range(lowres.shape[0]):
        for j in range(lowres.shape[1]):
            if image_mask[i][j] > 0:
                sc_y.append(i / adata.uns['spatial'][library_id]['scalefactors'][scalefactor])
                sc_x.append(j / adata.uns['spatial'][library_id]['scalefactors'][scalefactor])
                sc_z.append(adjust_lightness(matplotlib.colors.to_hex(lowres[i][j]), lightness))

    filter_he = (sc_x > adata.obsm['spatial'].min(0)[0]) & (sc_x < adata.obsm['spatial'].max(0)[0]) & (sc_y > adata.obsm['spatial'].min(0)[1]) & (sc_y < adata.obsm['spatial'].max(0)[1])

    if filter_spot is None:
        sc_x = np.array(sc_x)[filter_he]
        sc_y = np.array(sc_y)[filter_he]
        sc_z = np.array(sc_z)[filter_he]
        save_list = [sc_x, sc_y, sc_z]

        import pickle
        if not os.path.exists(str(save_directory)):
            os.makedirs(str(save_directory))
        with open(str(save_directory) + "/HE.pickle", "wb") as output_file:
            pickle.dump(save_list, output_file)
    else:
        sc_x = np.array(sc_x)[filter_he & filter_spot]
        sc_y = np.array(sc_y)[filter_he & filter_spot]
        sc_z = np.array(sc_z)[filter_he & filter_spot]
        save_list = [sc_x, sc_y, sc_z]

        import pickle
        if not os.path.exists(str(save_directory)):
            os.makedirs(str(save_directory))
        with open(str(save_directory) + "/HE_filter.pickle", "wb") as output_file:
            pickle.dump(save_list, output_file)
    
    print('HE information saved in ' + str(save_directory))



def tl_denoise(adata=None, plot_name=None, find_obs=False, use_raw=True, input_value=None, Mesh_directory='./Ridge_tmp', save_directory=None, iteration=10, L1=2/3, weight=0.9):
    import pickle
    if os.path.isfile(str(Mesh_directory) + "/level1_list.pickle"):
        with open(str(Mesh_directory) + "/level1_list.pickle", "rb") as input_file:
            level1_list = pickle.load(input_file)
            
        with open(str(Mesh_directory) + "/level2_list.pickle", "rb") as input_file:
            level2_list = pickle.load(input_file)
    else:
        print("Please run RidgeSpace.tl_mesh() with parameter 'mesh_optimization=True' first.")
        return 0

    if save_directory is None:
        save_directory = Mesh_directory

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

            original_value = matrix[:,var_names == plot_name].A.ravel()

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
    

def pl_single(adata=None, input_z=None, plot_name='Gene', normalize_z=True, out_path=None, obs_cluster=None, color_map=None, cluster=None, Mesh_directory='./Ridge_tmp', lw=1, height_ratio=1.5, elev=35, azim=-70, z_height=1, bg_alpha=1, plot_HE=False, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_bg=False, HE_z=None, plot_clustering=False, clustering_size=1, clustering_alpha=0.1, plot_clustering_boundary=False, clustering_boundary_lw=1, clustering_boundary_alpha=0.4, select_c=None, xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, trunc=None, figsize_set=(8,8)):
    import os
    import pickle
    import sys
    if input_z is None:
        files = os.listdir(Mesh_directory)
        if str(plot_name) + "_denoise.pickle" in files:
            with open(str(Mesh_directory) + "/" + str(plot_name) + "_denoise.pickle", "rb") as input_file:
                input_z = pickle.load(input_file)
        else:
            print("No file named " + str(plot_name) + " in " + str(Mesh_directory))
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

                    input_z = matrix[:,var_names == plot_name].A.ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_z = input_z

    if normalize_z == True:
        norm_z = 100 * input_z / max(input_z)#z / np.linalg.norm(z)
    else:
        norm_z = input_z

    if trunc is not None:
        tmp = []
        for mv in norm_z:
            if mv > trunc:
                tmp.append(mv)
            else:
                tmp.append(0)
        norm_z = np.array(tmp)

    ##read triang mesh
    import pickle
    if os.path.isfile(str(Mesh_directory) + "/triang.pickle"):
        with open(str(Mesh_directory) + "/triang.pickle", "rb") as input_file:
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

    colors, original_colors = color_alpha(input_z=norm_z, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alpha)

    fig = plt.figure(figsize=figsize_set)
    ax = fig.add_subplot(1,1,1, projection='3d', computed_zorder=False)

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'

    if select_c is not None:
        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colors])
        colors = colors_new

    plot_Ridge(ax, triang, norm_z, lw = lw, facecolors = colors, antialiased=False, rasterized=True)

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
                

        ax.scatter3D(sc_x, sc_y, zs=height_ratio*max(norm_z), facecolors=sc_z, s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True)

        if HE_bg == True:
            if HE_z is None:
                HE_z = 0
            ax.scatter3D(sc_x, sc_y, zs=HE_z, facecolors=sc_z, s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=0)

        
    if plot_clustering == True:
        ax.scatter3D(xs=triang.x, ys=triang.y, zs=height_ratio*max(norm_z), s=clustering_size, facecolors=np.array(original_colors), alpha=clustering_alpha, antialiased=False, linewidth=0, rasterized=True)#

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

    ax.view_init(elev=elev, azim=azim)

    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_zlim(zlim_min, zlim_max)

    ax.set_box_aspect((1, 1 * (max(ax.get_ylim3d()) - min(ax.get_ylim3d())) / (max(ax.get_xlim3d()) - min(ax.get_xlim3d())), z_height))
    ax.set_xlabel('Spatial 1', labelpad=-6)
    ax.set_ylabel('Spatial 2', labelpad=-6)
    ax.set_zlabel(str(plot_name) + 'expression')

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
    plt.show()

def pl_multipleOUT(adata=None, input_zA=None, input_zB=None, plot_nameA='Gene', plot_nameB='Gene', normalize_z=True, out_path=None, obs_cluster=None, color_map=None, cluster=None, Mesh_directory='./Ridge_tmp', lw=1, elev=35, azim=-70, z_height=1, bg_alphaA=1, bg_alphaB=1, plot_HE_bg=False, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, plot_clustering=False, clustering_size=1, clustering_alpha=0.1, select_c=None, xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, truncA=None, truncB=None, figsize_set=(8,8)):
    import os
    import pickle
    import sys
    if input_zA is None:
        files = os.listdir(Mesh_directory)
        if str(plot_nameA) + "_denoise.pickle" in files:
            with open(str(Mesh_directory) + "/" + str(plot_nameA) + "_denoise.pickle", "rb") as input_file:
                input_zA = pickle.load(input_file)
        else:
            print("No file named " + str(plot_nameA) + " in " + str(Mesh_directory))
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

                    input_zA = matrix[:,var_names == plot_nameA].A.ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_zA = input_zA

    if input_zB is None:
        files = os.listdir(Mesh_directory)
        if str(plot_nameB) + "_denoise.pickle" in files:
            with open(str(Mesh_directory) + "/" + str(plot_nameB) + "_denoise.pickle", "rb") as input_file:
                input_zB = pickle.load(input_file)
        else:
            print("No file named " + str(plot_nameB) + " in " + str(Mesh_directory))
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

                    input_zB = matrix[:,var_names == plot_nameB].A.ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_zB = input_zB

    if normalize_z == True:
        norm_zA = 100 * input_zA / max(input_zA)#z / np.linalg.norm(z)
        norm_zB = 100 * input_zB / max(input_zB)#z / np.linalg.norm(z)
    else:
        norm_zA = input_zA
        norm_zB = input_zB

    ## trunc norm_z
    if truncA is not None:
        tmp = []
        for mv in norm_zA:
            if mv > truncA:
                tmp.append(mv)
            else:
                tmp.append(0)
        norm_zA = np.array(tmp)
    
    if truncB is not None:
        tmp = []
        for mv in norm_zB:
            if mv > truncB:
                tmp.append(mv)
            else:
                tmp.append(0)
        norm_zB = np.array(tmp)    


    ##read triang mesh
    import pickle
    if os.path.isfile(str(Mesh_directory) + "/triang.pickle"):
        with open(str(Mesh_directory) + "/triang.pickle", "rb") as input_file:
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

    colorsA, original_colorsA = color_alpha(input_z=norm_zA, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alphaA)
    colorsB, original_colorsB = color_alpha(input_z=norm_zB, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alphaB)

    fig = plt.figure(figsize=figsize_set)
    ax = fig.add_subplot(1,1,1, projection='3d', computed_zorder=False)

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'

    if select_c is not None:
        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colorsA])
        colorsA = colors_new

        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colorsB])
        colorsB = colors_new

    tmp = -1 * norm_zB
    norm_zB = tmp

    plot_Ridge(ax, triang, norm_zA, lw = lw, facecolors = colorsA, antialiased=False, rasterized=True, zorder=1)
    plot_Ridge(ax, triang, norm_zB, lw = lw, facecolors = colorsB, antialiased=False, rasterized=True, zorder=0)
    
    if plot_HE_bg == True:
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
                
        if HE_z is None:
            HE_z = 0
        ax.scatter3D(sc_x, sc_y, zs=HE_z, facecolors=sc_z, s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=0)

    ax.view_init(elev=elev, azim=azim)
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

    ax.set_box_aspect((1, 1 * (max(ax.get_ylim3d()) - min(ax.get_ylim3d())) / (max(ax.get_xlim3d()) - min(ax.get_xlim3d())), z_height))
    ax.set_xlabel('Spatial 1', labelpad=-6)
    ax.set_ylabel('Spatial 2', labelpad=-6)
    ax.set_zlabel(str(plot_nameA) + ' & '+ str(plot_nameB) + ' expression')

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
    plt.show()


def pl_multipleIN(adata=None, input_zA=None, input_zB=None, plot_nameA='Gene', plot_nameB='Gene', height=100, normalize_z=True, out_path=None, obs_cluster=None, color_map=None, cluster=None, Mesh_directory='./Ridge_tmp', lw=1, elev=35, azim=-70, z_height=1, bg_alphaA=1, bg_alphaB=1,  plot_HE=False, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, plot_clustering=False, clustering_size=1, clustering_alpha=0.1, plot_clustering_boundary=False, clustering_boundary_lw=1, clustering_boundary_alpha=0.4, select_c=None,xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, truncA=None, truncB=None, HE_bg=False, figsize_set=(8,8)):
    import os
    import pickle
    import sys
    if input_zA is None:
        files = os.listdir(Mesh_directory)
        if str(plot_nameA) + "_denoise.pickle" in files:
            with open(str(Mesh_directory) + "/" + str(plot_nameA) + "_denoise.pickle", "rb") as input_file:
                input_zA = pickle.load(input_file)
        else:
            print("No file named " + str(plot_nameA) + " in " + str(Mesh_directory))
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

                    input_zA = matrix[:,var_names == plot_nameA].A.ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_zA = input_zA

    if input_zB is None:
        files = os.listdir(Mesh_directory)
        if str(plot_nameB) + "_denoise.pickle" in files:
            with open(str(Mesh_directory) + "/" + str(plot_nameB) + "_denoise.pickle", "rb") as input_file:
                input_zB = pickle.load(input_file)
        else:
            print("No file named " + str(plot_nameB) + " in " + str(Mesh_directory))
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

                    input_zB = matrix[:,var_names == plot_nameB].A.ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_zB = input_zB

    if normalize_z == True:
        norm_zA = 100 * input_zA / max(input_zA)#z / np.linalg.norm(z)
        norm_zB = 100 * input_zB / max(input_zB)#z / np.linalg.norm(z)
    else:
        norm_zA = input_zA
        norm_zB = input_zB

    ## trunc norm_z
    if truncA is not None:
        tmp = []
        for mv in norm_zA:
            if mv > truncA:
                tmp.append(mv)
            else:
                tmp.append(0)
        norm_zA = np.array(tmp)
    
    if truncB is not None:
        tmp = []
        for mv in norm_zB:
            if mv > truncB:
                tmp.append(mv)
            else:
                tmp.append(0)
        norm_zB = np.array(tmp)    


    ##read triang mesh
    import pickle
    if os.path.isfile(str(Mesh_directory) + "/triang.pickle"):
        with open(str(Mesh_directory) + "/triang.pickle", "rb") as input_file:
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

    colorsA, original_colorsA = color_alpha(input_z=norm_zA, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alphaA)
    colorsB, original_colorsB = color_alpha(input_z=norm_zB, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alphaB)

    fig = plt.figure(figsize=figsize_set)
    ax = fig.add_subplot(1,1,1, projection='3d', computed_zorder=False)

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'

    if select_c is not None:
        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colorsA])
        colorsA = colors_new
        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colorsB])
        colorsB = colors_new


    tmp = height - norm_zA
    norm_zA = tmp

    plot_Ridge(ax, triang, norm_zA, lw = lw, facecolors = colorsA, antialiased=False, rasterized=True, zorder=2)
    plot_Ridge(ax, triang, norm_zB, lw = lw, facecolors = colorsB, antialiased=False, rasterized=True, zorder=1)
    
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
                

        ax.scatter3D(sc_x, sc_y, zs=height, facecolors=sc_z, s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=3)
        if HE_bg == True:
            if HE_z is None:
                HE_z = 0
            ax.scatter3D(sc_x, sc_y, zs=HE_z, facecolors=sc_z, s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=0)

    if plot_clustering == True:
        ax.scatter3D(xs=triang.x, ys=triang.y, zs=height, s=clustering_size, facecolors=np.array(original_colorsA), alpha=clustering_alpha, antialiased=False, linewidth=0, rasterized=True, zorder=3)#

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

    ax.view_init(elev=elev, azim=azim)
    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_zlim(zlim_min, zlim_max)
    
    ax.set_box_aspect((1, 1 * (max(ax.get_ylim3d()) - min(ax.get_ylim3d())) / (max(ax.get_xlim3d()) - min(ax.get_xlim3d())), z_height))
    ax.set_xlabel('Spatial 1', labelpad=-6)
    ax.set_ylabel('Spatial 2', labelpad=-6)
    ax.set_zlabel(str(plot_nameA) + ' & '+ str(plot_nameB) + ' expression')

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
    plt.show()

def pl_trajectory(adata=None, input_z=None, find_obs=True, use_raw=True, normalize_z=False, out_path=None, obs_cluster=None, color_map=None, cluster=None, Mesh_directory='./Ridge_tmp', lw=1, height_ratio=1.5, elev=35, azim=-70, z_height=1, bg_alpha=1, plot_name='Gene', plot_HE_bg=True, HE_directory='./Ridge_tmp', HE_filter=False, HE_size=1, HE_alpha=1, HE_z=None, select_c=None, xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None, zlim_min=None, zlim_max=None, trunc=None, figsize_set=(8,8), color_arrow="black", lw_arrow=1, mutation_scale=15, density_arrow = .08, smooth_arrow = 1, neighbors_arrow = 50, splice_arrow=50):
    import os
    import pickle
    import sys
    if input_z is None:
        files = os.listdir(Mesh_directory)
        if str(plot_name) + "_denoise.pickle" in files:
            with open(str(Mesh_directory) + "/" + str(plot_name) + "_denoise.pickle", "rb") as input_file:
                input_z = pickle.load(input_file)
        else:
            print("No file named " + str(plot_name) + " in " + str(Mesh_directory))
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

                    input_z = matrix[:,var_names == plot_name].A.ravel()
            else:
                print("Please check input value for plotting or use parameter 'input_z'")
                return 0
    else:
        input_z = input_z

    if normalize_z == True:
        norm_z = 100 * input_z / max(input_z)#z / np.linalg.norm(z)
    else:
        norm_z = input_z

    if trunc is not None:
        tmp = []
        for mv in norm_z:
            if mv > trunc:
                tmp.append(mv)
            else:
                tmp.append(0)
        norm_z = np.array(tmp)

    ##read triang mesh
    import pickle
    if os.path.isfile(str(Mesh_directory) + "/triang.pickle"):
        with open(str(Mesh_directory) + "/triang.pickle", "rb") as input_file:
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

    colors, original_colors = color_alpha(input_z=norm_z, cluster=cluster, cluster_colors=cluster_colors, bg_alpha=bg_alpha)

    fig = plt.figure(figsize=figsize_set)
    ax = fig.add_subplot(1,1,1, projection='3d', computed_zorder=False)

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'

    if select_c is not None:
        colors_new = np.array([i if i[:-2] in select_c else '#cccccc' + i[-2:] for i in colors])
        colors = colors_new

    plot_Ridge(ax, triang, norm_z, lw = lw, facecolors = colors, antialiased=False, rasterized=True, zorder=1)

    X_i, Y_i, Z_i, U_i, V_i, W_i = infer_trajectory(triang, norm_z, density = density_arrow, smooth = smooth_arrow, neighbors = neighbors_arrow, splice = splice_arrow)
    for i in range(len(X_i)):
        arw = Arrow3D([X_i[i], U_i[i]],[Y_i[i], V_i[i]],[Z_i[i], W_i[i]], arrowstyle="-|>", color=color_arrow, lw = lw_arrow, mutation_scale=mutation_scale, zorder=1)
        ax.add_artist(arw)

    if plot_HE_bg == True:
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
                
        if HE_z is None:
            HE_z = 0
        ax.scatter3D(sc_x, sc_y, zs=HE_z, facecolors=sc_z, s=HE_size, alpha=HE_alpha, linewidth=0, antialiased=False, rasterized=True, zorder=0)

    ax.view_init(elev=elev, azim=azim)

    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_zlim(zlim_min, zlim_max)

    ax.set_box_aspect((1, 1 * (max(ax.get_ylim3d()) - min(ax.get_ylim3d())) / (max(ax.get_xlim3d()) - min(ax.get_xlim3d())), z_height))
    ax.set_xlabel('Spatial 1', labelpad=-6)
    ax.set_ylabel('Spatial 2', labelpad=-6)
    ax.set_zlabel(str(plot_name) + 'trajectory')

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
    plt.show()
    