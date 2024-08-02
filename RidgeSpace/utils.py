import pandas as  pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import matplotlib.tri as mtri
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.axes3d import *

import sys
import os

def plot_Ridge(self, *args, norm=None, vmin=None, vmax=None, lightsource=None, **kwargs):
        """
        Create a 3D RidgeSpace plot.

        Parameters
        ----------
        X, Y, Z : 2D arrays
            Data values.

        rcount, ccount : int
            Maximum number of samples used in each direction.  If the input
            data is larger, it will be downsampled (by slicing) to these
            numbers of points.  Defaults to 50.

        cmap : Colormap
            Colormap of the surface patches.

        facecolors : array-like of colors.
            Colors of each individual patch.

        norm : Normalize
            Normalization for the colormap.

        vmin, vmax : float
            Bounds for the normalization.

        shade : bool, default: True
            Whether to shade the facecolors.  Shading is always disabled when
            *cmap* is specified.

        lightsource : `~matplotlib.colors.LightSource`
            The lightsource to use when *shade* is True.

        """

        had_data = self.has_data()
        
        ##Triangle
        tri, args, kwargs = \
            Triangulation.get_from_args_and_kwargs(*args, **kwargs)
        try:
            z = kwargs.pop('Z')
            if z.ndim != 1:
                raise ValueError("Argument z must be 1-dimensional.")
        except KeyError:
            # We do this so Z doesn't get passed as an arg to PolyCollection
            z, *args = args
        z = np.asarray(z)

        triangles = tri.get_masked_triangles()
        xt = tri.x[triangles]
        yt = tri.y[triangles]
        zt = z[triangles]
        verts = np.stack((xt, yt, zt), axis=-1)
        ##--
        
        
        fcolors = kwargs.pop('facecolors', None)

        cmap = kwargs.get('cmap', None)
        shade = kwargs.pop('shade', cmap is None)
        if shade is None:
            raise ValueError("shade cannot be None.")
        
        colset = []# the sampled facecolor

        if fcolors is not None:
            
            # get triangles index
            masked_triangles = tri.get_masked_triangles().astype(int)

            # from fcolors
            fcolors_selected = np.array(fcolors)[masked_triangles]

            # zip for each colors
            fcolors_tuples = [tuple(row) for row in fcolors_selected]

            # Counter for mode
            from collections import Counter
            mode_result = [Counter(row).most_common(1)[0][0] for row in fcolors_tuples]

            fcolors1 = mode_result

            colset=list(fcolors1)#
            polyc = art3d.Poly3DCollection(verts, edgecolors=colset, facecolor=colset, *args, **kwargs)#
            if isinstance(verts, np.ndarray):
                avg_z = verts[:, :, 2].mean(axis=1)

            if vmin is not None or vmax is not None:
                polyc.set_clim(vmin, vmax)
            if norm is not None:
                polyc.set_norm(norm)
                
        elif cmap:
            ## from tri_surface
            polyc = art3d.Poly3DCollection(verts, *args, **kwargs)
            # can't always vectorize, because polys might be jagged
            if isinstance(verts, np.ndarray):
                avg_z = verts[..., 2].mean(axis=-1)

            polyc.set_array(avg_z)
            
            if vmin is not None or vmax is not None:
                polyc.set_clim(vmin, vmax)
            if norm is not None:
                polyc.set_norm(norm)

        self.add_collection(polyc)
        self.auto_scale_xyz(tri.x, tri.y, z, had_data)

        return polyc
    
def color_variant(hex_color, brightness_offset=2):
    """ takes a color like #87c95f and produces a lighter or darker variant """
    if len(hex_color) != 7:
        raise Exception("Passed %s into color_variant(), needs to be in #87c95f format." % hex_color)
    rgb_hex = [hex_color[x:x+2] for x in [1, 3, 5]]
    new_rgb_int = [int(hex_value, 16) + brightness_offset for hex_value in rgb_hex]
    new_rgb_int = [min([255, max([0, i])]) for i in new_rgb_int] # between 0 and 255

    out = '#'
    for i in new_rgb_int:
        hex_c = hex(i)[2:]
        if len(hex_c) == 2:
            out += hex_c
        else:
            out += '0' + hex_c
    return out

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return mc.to_hex(colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2]), keep_alpha=False)

def rgb2gray(rgb):
    import numpy as np
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140]).astype('uint8') 

def get_HE_mask(img_HE):
    import numpy as np
    ## change
    from skimage import data,filters
    import cv2
    image=rgb2gray(img_HE)
    thresh = filters.threshold_otsu(image)   #返回一个阈值

    _, image_b = cv2.threshold(image, thresh, 250, cv2.THRESH_BINARY_INV)
    mask=np.zeros([image.shape[0],image.shape[1]])
    # find all contours
    contours, _ = cv2.findContours(image_b, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    area = []
    # find largest
    for k in range(len(contours)):
        area.append(cv2.contourArea(contours[k]))
    max_idx = np.argmax(np.array(area))
    # get largest
    mask = cv2.drawContours(mask, contours, max_idx, 1, cv2.FILLED)
    return mask, contours[max_idx]

def mesh_optimize(triang, save_directory, degree = 2):
    import numpy as np
    from scipy.sparse import lil_matrix

    num_vertices = triang.x.shape[0]
    level1_list = []
    level2_list = []

    # Create a sparse adjacency matrix for triangles
    adj_matrix = lil_matrix((num_vertices, num_vertices), dtype=bool)
    for idx, triangle in enumerate(triang.triangles):
        for v in triangle:
            adj_matrix[v, triangle] = True

    # Compute level1_list and level2_list
    for i in range(num_vertices):
        # Get adjacent triangles of vertex i
        adj_triangles_i = adj_matrix[i].rows[0]
        
        # Level 1: adjacent triangles excluding triangle i itself
        level1_set = set(adj_triangles_i) - {i}
        level1_list.append(level1_set)
        
        # Level 2: neighbors of neighbors excluding level 1
        level2_set = set()
        for j in level1_set:
            adj_triangles_j = adj_matrix[j].rows[0]
            level2_set.update(adj_triangles_j)
        level2_set -= {i}  # Remove triangle i itself
        level2_set -= level1_set  # Remove level 1 triangles
        level2_list.append(level2_set)

    # Convert lists of sets to lists of lists (if needed)
    level1_list = [list(s) for s in level1_list]
    level2_list = [list(s) for s in level2_list]

    import pickle
    if not os.path.exists(str(save_directory)):
        os.makedirs(str(save_directory))
    with open(str(save_directory) + "/level1_list.pickle", "wb") as output_file:
        pickle.dump(level1_list, output_file)
    with open(str(save_directory) + "/level2_list.pickle", "wb") as output_file:
        pickle.dump(level2_list, output_file)

    print('Mesh optimization done')

def molecular_denoise(original_value, level1_list, level2_list, iteration=10, L1=2/3, weight=0.9):

    z = original_value
    tmp_z = z.copy()
    for i in range(iteration):
        # print(i)
        dz1 = []
        dz2 = []
        for m in range(z.shape[0]):
            dz1i = tmp_z[list(level1_list[m])].sum() / len(level1_list[m])
            dz1.append(dz1i)

            dz2i = tmp_z[list(level2_list[m])].sum() / len(level2_list[m])
            dz2.append(dz2i)   

        dz1_a = np.array(dz1)
        dz2_a = np.array(dz2)

        dz = L1 * dz1_a + (1-L1) * dz2_a
        tmp_z = weight*tmp_z + (1-weight)*dz

    return(tmp_z)


def color_alpha(input_z, cluster, cluster_colors, bg_alpha=1):
    # change color alpha based on z value

    color_list = [color_variant(i[:7], 2) for i in cluster_colors]
    colors_pd = pd.DataFrame(np.array(cluster)).applymap(lambda x: dict(zip(cluster.categories, color_list))[x])#
    colors_c = colors_pd[0].values.copy()
    add_alpha=bg_alpha

    norm1 = input_z
    face_alphas = 1*norm1 / max(norm1)
    
    expected = { 100: 'FF', 99: 'FC', 98: 'FA', 97: 'F7',
    96: 'F5',95: 'F2', 94: 'F0',93: 'ED', 92: 'EB',91: 'E8', 90: 'E6',89: 'E3',  
    88: 'E0',87: 'DE',   86: 'DB',85: 'D9',
    84: 'D6',83: 'D4',    82: 'D1',81: 'CF',
    80: 'CC',79: 'C9',    78: 'C7',77: 'C4',
    76: 'C2',75: 'BF',    74: 'BD',73: 'BA',
    72: 'B8',71: 'B5',
    70: 'B3',69: 'B0',
    68: 'AD',67: 'AB',
    66: 'A8', 65: 'A6',
    64: 'A3',63: 'A1',
    62: '9E', 61: '9C',
    60: '99', 59: '96',
    58: '94', 57: '91',
    56: '8F', 55: '8C',
    54: '8A', 53: '87', 52: '85', 51: '82',
    50: '80', 49: '7D',   48: '7A', 47: '78',
    46: '75', 45: '73',  44: '70', 43: '6E',
    42: '6B', 41: '69', 40: '66', 39: '63', 
    38: '61', 37: '5E', 36: '5C', 35: '59',
    34: '57', 33: '54',    32: '52', 31: '4F',
    30: '4D', 29: '4A', 28: '47', 27: '45',
    26: '42', 25: '40',
    24: '3D', 23: '3B',
    22: '38', 21: '36',
    20: '33', 19: '30', 18: '2E',  17: '2B',
    16: '29', 15: '26',
    14: '24', 13: '21',
    12: '1F', 11: '1C',
    10: '1A',  9: '17',
    8: '14', 7: '12', 6: '0F', 5: '0D',
    4: '0A', 3: '08', 2: '05', 1: '03', 0: '00'}

    for i in range(colors_c.shape[0]):
        if colors_c[i] != 'None':
            colors_c[i] = colors_c[i] + expected[min(100, add_alpha + int(100 * (np.abs(face_alphas[i]))) )]

    return([colors_c, colors_pd[0].values])

def sort_xvgrid(i, neighbors, triang, z):
    tx, ty = i
    x = triang.x - tx
    y = triang.y - ty
    xy = np.abs(x) + np.abs(y)
    zi = np.argmin(xy)#[0]
    
    min_neighbor = np.partition(xy, neighbors)[:neighbors]
    neighbor_index = np.where(xy <= min_neighbor.max())[0]
    
    Temp_zi = np.argmax((z[zi] - z[neighbor_index]))#np.abs
    Vzi = neighbor_index[Temp_zi]
    return((zi, Vzi))

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
#         FancyArrowPatch.draw(self, renderer)
        return np.min(zs)

def infer_trajectory(triang, norm_z, density = .08, smooth = 1, neighbors = 50, splice=50):
    z = norm_z

    grs = []
    m, M = np.min(triang.x), np.max(triang.x)
    m = m - 0.01 * np.abs(M - m)
    M = M + 0.01 * np.abs(M - m)
    gr = np.linspace(m, M, int(50 * density))
    grs.append(gr)

    m, M = np.min(triang.y), np.max(triang.y)
    m = m - 0.01 * np.abs(M - m)
    M = M + 0.01 * np.abs(M - m)
    gr = np.linspace(m, M, int(50 * density))
    grs.append(gr)
    meshes_tuple = np.meshgrid(*grs)
    X_grid = np.vstack([i.flat for i in meshes_tuple]).T

    scale = np.mean([(g[1] - g[0]) for g in grs]) * smooth

    min_neighborM = np.partition(z, 5)[-5:]
    neighbor_indexM = np.where(z >= min_neighborM.min())[0][::splice]


    OX_i = []
    OY_i = []
    OZ_i = []
    OU_i = []
    OV_i = []
    OW_i = []

    for neighbor_indexi in neighbor_indexM:
        X_i = []
        Y_i = []
        Z_i = []
        U_i = []
        V_i = []
        W_i = []
        i = [triang.x[neighbor_indexi], triang.y[neighbor_indexi]]
        count = 0
        while 1:
            count += 1
            zi, Vzi = sort_xvgrid(i, neighbors, triang, z)

            if z[zi] == z[Vzi]:
                zi, Vzi = sort_xvgrid(i, neighbors+1, triang, z)
                break

            if z[zi] > z[Vzi]:
                X_i.append(triang.x[zi])
                Y_i.append(triang.y[zi])
                Z_i.append(z[zi])
                U_i.append(triang.x[Vzi])#-triang.x[zi])
                V_i.append(triang.y[Vzi])#-triang.y[zi])
                W_i.append(z[Vzi])#-z[zi])
            elif z[zi] < z[Vzi]:
                X_i.append(triang.x[Vzi])
                Y_i.append(triang.y[Vzi])
                Z_i.append(z[Vzi])
                U_i.append(triang.x[zi])#-triang.x[Vzi])
                V_i.append(triang.y[zi])#-triang.y[Vzi])
                W_i.append(z[zi])#-z[Vzi])

            if W_i[-1] == z.min():
                count = 0
                break
            else:
                i = [U_i[-1], V_i[-1]]
            if count > 30:
                count = 0
                break
        
        OX_i += X_i
        OY_i += Y_i
        OZ_i += Z_i
        OU_i += U_i
        OV_i += V_i
        OW_i += W_i

    for neighbor_indexi in X_grid:#neighbor_indexM:
        X_i = []
        Y_i = []
        Z_i = []
        U_i = []
        V_i = []
        W_i = []
        i = neighbor_indexi#[triang.x[neighbor_indexi], triang.y[neighbor_indexi]]
        count = 0
        while 1:
            count += 1
            zi, Vzi = sort_xvgrid(i, neighbors, triang, z)

            if z[zi] == z[Vzi]:
                zi, Vzi = sort_xvgrid(i, neighbors+1, triang, z)
                break

            if z[zi] > z[Vzi]:
                X_i.append(triang.x[zi])
                Y_i.append(triang.y[zi])
                Z_i.append(z[zi])
                U_i.append(triang.x[Vzi])#-triang.x[zi])
                V_i.append(triang.y[Vzi])#-triang.y[zi])
                W_i.append(z[Vzi])#-z[zi])
            elif z[zi] < z[Vzi]:
                X_i.append(triang.x[Vzi])
                Y_i.append(triang.y[Vzi])
                Z_i.append(z[Vzi])
                U_i.append(triang.x[zi])#-triang.x[Vzi])
                V_i.append(triang.y[zi])#-triang.y[Vzi])
                W_i.append(z[zi])#-z[Vzi])

            if W_i[-1] == z.min():
                count = 0
                break
            else:
                i = [U_i[-1], V_i[-1]]
            if count > 50:
                count = 0
                break
        
        OX_i += X_i
        OY_i += Y_i
        OZ_i += Z_i
        OU_i += U_i
        OV_i += V_i
        OW_i += W_i

    X_i = OX_i
    Y_i = OY_i
    Z_i = OZ_i
    U_i = OU_i
    V_i = OV_i
    W_i = OW_i

    return([X_i, Y_i, Z_i, U_i, V_i, W_i])