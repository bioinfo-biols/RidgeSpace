import pandas as  pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import axes3d, art3d
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection 
# from mpl_toolkits.mplot3d.axes3d import *  Axes3D,

import sys
import os

def plot_Ridge(ax, *args, norm=None, vmin=None, vmax=None, lightsource=None, **kwargs):
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

        had_data = ax.has_data()

        ##Triangle
        tri, args, kwargs = \
            axes3d.Triangulation.get_from_args_and_kwargs(*args, **kwargs)
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

        ax.add_collection(polyc)
        ax.auto_scale_xyz(tri.x, tri.y, z, had_data)

        return polyc


class MyAxes3D(object):
    ## modified based on https://stackoverflow.com/questions/15042129/changing-position-of-vertical-z-axis-of-3d-plot
    def __init__(self, baseObject, sides_to_draw):
        self.__class__ = type(baseObject.__class__.__name__,
                              (self.__class__, baseObject.__class__),
                              {})
        self.__dict__ = baseObject.__dict__
        self.sides_to_draw = list(sides_to_draw)
        self.mouse_init()
        self.left_ticks = []
        self.left_ticklabels = []
        self.right_ticks = []
        self.right_ticklabels = []
        self.left_label = ""
        self.right_label = ""

    def set_custom_ticks(self, side, ticks=None, ticklabels=None):
        """
        Set custom ticks and ticklabels for the specified side ('l' or 'r').
        """
        if side == 'l':
            self.left_ticks = ticks
            self.left_ticklabels = ticklabels

        elif side == 'r':
            self.right_ticks = ticks
            self.right_ticklabels = ticklabels

    def set_zaxis_label(self, side, label):
        """
        Set a custom label for the specified z-axis side ('l' or 'r').
        """
        if side == 'l':
            self.left_label = label
        elif side == 'r':
            self.right_label = label

    def draw(self, renderer):
        # Hide default z-axis features
        self.zaxis.set_ticks([])  # remove original ticks
        self.zaxis.set_ticklabels([])  # remove original ticklabels  
        self.zaxis.line.set_visible(False)
        self.zaxis.label.set_visible(False)
#         self.w_zaxis.pane.set_visible(False)
        # Draw axes
        super(MyAxes3D, self).draw(renderer)

        zaxis = self.zaxis
        draw_grid_old = zaxis.axes._draw_grid
        zaxis.axes._draw_grid = False
        tmp_planes = zaxis._PLANES
        
        self.zaxis.label.set_visible(True)
        self.zaxis.line.set_visible(True)

        zaxis._PLANES = (tmp_planes[2], tmp_planes[3],
                         tmp_planes[0], tmp_planes[1],
                         tmp_planes[4], tmp_planes[5])
        if self.left_ticks:
            zaxis.set_ticks(self.left_ticks)
            if self.left_ticklabels and len(self.left_ticklabels) == len(self.left_ticks):
                zaxis.set_ticklabels(self.left_ticklabels)
        zaxis.label.set_text(self.left_label)  # Set left z-axis label       
        zaxis.draw(renderer)

        zaxis._PLANES = (tmp_planes[3], tmp_planes[2],
                         tmp_planes[1], tmp_planes[0],
                         tmp_planes[4], tmp_planes[5])

        if self.right_ticks:
            zaxis.set_ticks(self.right_ticks)
            if self.right_ticklabels and len(self.right_ticklabels) == len(self.right_ticks):
                zaxis.set_ticklabels(self.right_ticklabels)
        zaxis.label.set_text(self.right_label)  # Set right z-axis label
        zaxis.draw(renderer)
        
        # Reset planes and grid
        zaxis._PLANES = tmp_planes
        zaxis.axes._draw_grid = draw_grid_old


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

def adjust_lightness_batch(rgb_array, amount=0.5):
    rgb_array = np.clip(rgb_array, 0, 1)
    # 按比例调整亮度
    return np.clip(rgb_array * amount, 0, 1)

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

def get_HE_mask(img_HE, scale=1):
    import numpy as np
    ## change
    from skimage import data,filters
    import cv2
    image=rgb2gray(img_HE)
    thresh = filters.threshold_otsu(image)   #返回一个阈值

    _, image_b = cv2.threshold(image, thresh * scale, 250, cv2.THRESH_BINARY_INV)
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


def find_neighbor0(triang):
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
    return([level1_list, level2_list])

def find_neighbor1(triang):
    neighbor1 = []
    neighbor2 = []
    triange1 = triang.edges[:,0]
    triange2 = triang.edges[:,1]
    for ti in range(len(triang.x)):
        tmp = np.unique(triang.edges[np.in1d(triange1, ti) + np.in1d(triange2, ti),:].flatten())
        neighbor1.append([ i for i in tmp if i not in [ti, -1]])

        if len(neighbor1[-1]) > 1:
            tmp = np.unique(triang.edges[np.in1d(triange1, neighbor1[-1]) + np.in1d(triange2, neighbor1[-1]),:].flatten())
            neighbor2.append([ i for i in tmp if i not in neighbor1[-1] + [-1]])

    return([neighbor1, neighbor2])

import numpy as np
def find_neighbor(triang):
    n = len(triang.x)
    
    # 提取三角形的边信息
    triange1 = triang.edges[:, 0]
    triange2 = triang.edges[:, 1]
    
    # 通过字典记录每个节点的邻居
    neighbor1 = {i: set() for i in range(n)}
    neighbor2 = {i: set() for i in range(n)}
    
    # 记录每条边的两个端点，构建邻接关系
    for e1, e2 in zip(triange1, triange2):
        neighbor1[e1].add(e2)
        neighbor1[e2].add(e1)
    
    # 获取一级邻居
    for ti in range(n):
        neighbor1[ti].discard(ti)  # 移除自身
        # neighbor1[ti].discard(-1)

    # 获取二级邻居
    for ti in range(n):
        # 通过一级邻居的邻居来找二级邻居
        second_degree_neighbors = set()
        for n1 in neighbor1[ti]:
            second_degree_neighbors.update(neighbor1[n1])
        
        # second_degree_neighbors.discard(-1)
        second_degree_neighbors.discard(ti)  # 移除自身
        second_degree_neighbors -= neighbor1[ti]  # 移除一级邻居
        neighbor2[ti] = second_degree_neighbors
    
    # 将结果转换为列表形式
    return [list(neighbor1[i]) for i in range(n)], [list(neighbor2[i]) for i in range(n)]

def mesh_optimize(triang, save_directory):
    
    neighbor1, neighbor2 = find_neighbor(triang)

    import pickle
    if not os.path.exists(str(save_directory)):
        os.makedirs(str(save_directory))
    with open(str(save_directory) + "/level1_list.pickle", "wb") as output_file:
        pickle.dump(neighbor1, output_file)
    with open(str(save_directory) + "/level2_list.pickle", "wb") as output_file:
        pickle.dump(neighbor2, output_file)

    print('Mesh optimization done')

def molecular_denoise(original_value, level1_list, level2_list, iteration=10, L1=2/3, weight=0.9):

    z = original_value
    tmp_z = z.copy()
    for i in range(iteration):
        # print(i)
        dz1 = []
        dz2 = []
        for m in range(z.shape[0]):
            if len(level1_list[m]) > 0:
                dz1i = tmp_z[list(level1_list[m])].sum() / len(level1_list[m])
                dz1.append(dz1i)
            else:
                dz1.append(tmp_z[m])

            if len(level2_list[m]) > 0:
                dz2i = tmp_z[list(level2_list[m])].sum() / len(level2_list[m])
                dz2.append(dz2i)
            else:
                dz2.append(tmp_z[m])  
                

        dz1_a = np.array(dz1)
        dz2_a = np.array(dz2)

        dz = L1 * dz1_a + (1-L1) * dz2_a
        tmp_z = weight*tmp_z + (1-weight)*dz

    return(tmp_z)



def color_alpha(input_z, cluster, cluster_colors, bg_alpha=1, color_brighten = True):
    # change color alpha based on z value
    if color_brighten == True:
        color_list = [color_variant(i[:7], 2) for i in cluster_colors]
    else:
        color_list = cluster_colors
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
            colors_c[i] = colors_c[i] + expected[min(100, max(0, add_alpha + int(100 * (np.abs(face_alphas[i])))))]

    return([colors_c, colors_pd[0].values])


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

def sort_xvgrid(i, indezi, Neighbors, Mask_mesh, z):
    tx, ty = i
    zi = indezi#np.argmin(xy)#[0]
    neighbor_tmp = np.array(Neighbors[indezi])
    if len(neighbor_tmp) > 0:
        neighbor_index = neighbor_tmp[Mask_mesh[neighbor_tmp, 2] == 0]
        if len(neighbor_index) > 0:
            Temp_zi = np.argmax((z[zi] - z[neighbor_index]))#np.abs
            Vzi = neighbor_index[Temp_zi]
            return((zi, Vzi))
        else:
            return((zi, zi))
    else:
        return((zi, zi))

def get_plotlabel(vmin, vmax):

    # 计算数据范围的大小
    range = vmax - vmin
    # 选择一个合适的间隔
    base = 10 ** np.floor(np.log10(range))# 选择基于数据范围的10的幂次间隔
    step = base# 初步选择间隔为基数

    # 保证间隔的“整洁性”
    if range / step > 10:
        step *= 4
    elif range / step > 5:
        step *= 2
    
    # 生成刻度的位置
    start = np.floor(vmin / step) * step
    stop = np.ceil(vmax / step) * step
    return [start, stop, step]#np.arange(start, stop, step)

def format_number(value):
    return f"{int(value)}" if value.is_integer() else f"{value}"

def get_plotlabel_trunc(vmin, vmax):
    import numpy as np
    # 计算数据范围的大小
    range = vmax - vmin
    # 选择一个合适的间隔
    base = 10 ** np.floor(np.log10(range))# 选择基于数据范围的10的幂次间隔
    step = base# 初步选择间隔为基数

    # 保证间隔的“整洁性”
    if range / step > 10:
        step *= 2
    elif range / step > 5:
        step *= 1.5

    # 生成刻度的位置
    start = np.floor(vmin / step) * step + np.floor(vmin / (step / 10)) * int((step / 10))
    stop = np.ceil(vmax / step) * step
    return [start, stop, step]#np.arange(start, stop, step)

def infer_trajectory(triang, norm_z, mesh_directory, density = 0.5, dense = True, min_length = 0.12, arrow_scale = 22, lw = 1, color = 'k'):
    z = norm_z
    spotsN = len(triang.x)
    Mask_mesh = np.vstack([triang.x, triang.y, [0]*triang.x.shape[0], range(spotsN)]).T
    Z_mesh = np.vstack([triang.x, triang.y, z, range(spotsN)]).T
    #Z_mesh_order = Z_mesh[:, 2].argsort()[::-1]

    split_k = int(int(np.sqrt(spotsN) / 2) * density)
    minlength = min_length * np.abs((z.max() - z.min()))
    neighbors = 1

    # 
    import pickle
    if os.path.isfile(str(mesh_directory) + "/level1_list.pickle"):
        with open(str(mesh_directory) + "/level1_list.pickle", "rb") as input_file:
            neighbor1 = pickle.load(input_file)
            
        with open(str(mesh_directory) + "/level2_list.pickle", "rb") as input_file:
            neighbor2 = pickle.load(input_file)
    else:
        neighbor1, neighbor2 = find_neighbor(triang)

    Neighbor_mesh = np.vstack([range(len(triang.x)), np.array(neighbor1, dtype=object), np.array(neighbor2, dtype=object)]).T
    
    
    grs = []
    m, M = np.min(triang.x), np.max(triang.x)
    gr1 = np.linspace(m, M, split_k)
    m, M = np.min(triang.y), np.max(triang.y)
    gr2 = np.linspace(m, M, split_k)

    Start = []
    for i in range(len(gr1)-1):
        for j in range(len(gr2)-1):
            Start_tmp = Z_mesh[(Z_mesh[:,0] >= gr1[i]) & (Z_mesh[:,0] < gr1[i+1]) & (Z_mesh[:,1] >= gr2[j]) & (Z_mesh[:,1] < gr2[j+1]),2:]
            if len(Start_tmp) > 0:
                Start += Start_tmp[Start_tmp[:,0].argsort()[::-1][0:1], 1].astype(int).tolist()
                
    traj = []
    for indexFi in Start:
        if Mask_mesh[indexFi,2] == 1:
            continue
        tmpx, tmpy, tmpz, tmpo = Z_mesh[indexFi,:]
        indexi = [tmpx, tmpy]
        route = []
        routeI = []
        route.append([tmpx, tmpy, tmpz])
        routeI.append(indexFi)
        count = 0
        while 1:
            count += 1
            zi, Vzi = sort_xvgrid(indexi, routeI[-1], Neighbor_mesh[:,neighbors], Mask_mesh, z)

            if z[zi] == z[Vzi]:
                if z[zi] < np.min(z) + 0.01 * (np.max(z) - np.min(z)):
                    break

                zi, Vzi = sort_xvgrid(indexi, routeI[-1], Neighbor_mesh[:,neighbors+1], Mask_mesh, z)
                if count > 20:
                    count = 0
                    break
                if z[zi] == z[Vzi]:
                    count = 0
                    break

            if z[zi] > z[Vzi]:
                route.append([triang.x[Vzi], triang.y[Vzi], z[Vzi]])#triang.x[zi], triang.y[zi], z[zi]])#, 
                routeI.append(Vzi)

            elif z[zi] < z[Vzi]:
                count = 0
                break

            if route[-1][2] == z.min():
                count = 0
                break
            else:
                indexi = [route[-1][0], route[-1][1]]

            if len(route) > 20:
                count = 0
                break

        if len(route) > 1:
            tx, ty, tz = np.array(route).T
            if np.abs(route[-1][2] - route[0][2]) > minlength: 
                #np.hypot(np.abs(route[-1][0] - route[0][0]), np.abs(route[-1][1] - route[0][1]))
                traj.append(np.array(route))
                for ri in routeI:
                    Mask_mesh[ri, 2] = 1
                    for rj in Neighbor_mesh[ri, 1]:
                        Mask_mesh[rj, 2] = 1
                        if dense == True:
                            for rk in Neighbor_mesh[rj, 1]:
                                Mask_mesh[rk, 2] = 1
                                
    streamlines = []
    plot_p = []
    for t in traj:
        tx, ty, tz = t.T
        # Create multiple tiny segments if varying width or color is given
        points = np.transpose([tx, ty, tz])
        streamlines.append(points)

        # Add arrows halfway along each trajectory.
        s = np.cumsum(np.abs(np.diff(tz) / (tz.max() - tz.min() + 1e-5)))#np.diff(tx) / (tx.max() - tx.min() + 1e-5), np.diff(ty) / (ty.max() - ty.min() + 1e-5), np.hyplot
        n = np.searchsorted(s, s[-1] / 2.)
        # if n > 0:
        #     arrow_tail = (tx[n-1], ty[n-1], tz[n-1])
        #     arrow_head = (tx[n], ty[n], tz[n])
        # else:
        #     arrow_tail = (tx[0], ty[0], tz[0])
        #     arrow_head = (tx[1], ty[1], tz[1])
        if len(tx) > n+1:
            arrow_tail = (tx[n], ty[n], tz[n])
            arrow_head = (tx[n+1], ty[n+1], tz[n+1])
        else:
            arrow_tail = (tx[n-1], ty[n-1], tz[n-1])
            arrow_head = (tx[n], ty[n], tz[n])


        # if len(tx) > n+1:
        #     arrow_head = (np.mean(tx[n+1:]), np.mean(ty[n+1:]), np.mean(tz[n+1:]))#n + 2
        # else:
        #     arrow_head = (np.mean(tx[n:]), np.mean(ty[n:]), np.mean(tz[n:]))

        # min_flag = 0
        # for ni in range(n, len(tx)-1):
        #     A = np.array([tx[ni], ty[ni], tz[ni]])
        #     B = np.array(arrow_head)
        #     C = np.array([tx[ni+1], ty[ni+1], tz[ni+1]])

        #     AB = B - A
        #     AC = C - A
        #     crossABC = np.cross(AB, AC)

        #     if np.allclose(crossABC, [0, 0, 0]):
        #         arrow_tail = (tx[ni], ty[ni], tz[ni])
        #         min_flag = 0
        #         break

        # if min_flag != 0:
        #     arrow_head = (tx[n+1], ty[n+1], tz[n+1])

        p = Arrow3D([arrow_tail[0], arrow_head[0]], [arrow_tail[1], arrow_head[1]], [arrow_tail[2], arrow_head[2]], arrowstyle="-|>", mutation_scale=arrow_scale, lw=0, color=color)#
        plot_p.append(p)
        
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    plot_lc = Line3DCollection(streamlines, colors=color, lw=lw)
    
    
    return([plot_p, plot_lc])


def read_zarr(file_path):
    import os
    import numpy as np

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file at path '{file_path}' does not exist.")

    # 检查文件是否为 .zip 格式
    if not file_path.endswith(".zip"):
        raise ValueError("Input file must be a .zip file.")
    
    try:
        import zarr
    except Exception as e:
        ImportError('Please pip install zarr.')
        
    try:
        import zipfile, tempfile
    except Exception as e:
        ImportError('Please install zipfile and tempfile.')   
    # 解压 .zarr.zip 文件到临时目录
    temp_dir = tempfile.TemporaryDirectory()
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(temp_dir.name)

    # 假设解压后的目录中含有 .zarr 文件夹
    zarr_folder = [f for f in os.listdir(temp_dir.name) if f.endswith(".zarr")]
    if not zarr_folder:
        raise RuntimeError("No .zarr folder found in the zip archive.")

    zarr_path = os.path.join(temp_dir.name, zarr_folder[0])

    # 打开 zarr 文件夹
    try:
        zarr_root = zarr.open(zarr_path, mode='r')
    except Exception as e:
        raise RuntimeError(f"Failed to read the zarr file: {e}")

    # 遍历 zarr 数据集，读取内容
    data = {}
    for key in zarr_root.keys():
        dataset = zarr_root[key]
        if isinstance(dataset, zarr.core.Array):
            data[key] = np.array(dataset)
        elif isinstance(dataset, zarr.hierarchy.Group):
            # 递归读取子组
            data[key] = read_zarr_group(dataset)
        else:
            raise TypeError(f"Unsupported dataset type for key '{key}'")

    # 清理临时目录
    temp_dir.cleanup()

    return data


def read_zarr_group(group):
    import zarr
    import numpy as np

    data = {}
    for key in group.keys():
        item = group[key]
        if isinstance(item, zarr.core.Array):
            data[key] = np.array(item)
        elif isinstance(item, zarr.hierarchy.Group):
            data[key] = read_zarr_group(item)
    return data


def sort_xvgrid_old(i, neighbors, triang, z):
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

def infer_trajectory_old(triang, norm_z, density = .08, smooth = 1, neighbors = 50, splice=50):
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
            zi, Vzi = sort_xvgrid_old(i, neighbors, triang, z)

            if z[zi] == z[Vzi]:
                zi, Vzi = sort_xvgrid_old(i, neighbors+1, triang, z)
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
            zi, Vzi = sort_xvgrid_old(i, neighbors, triang, z)

            if z[zi] == z[Vzi]:
                zi, Vzi = sort_xvgrid_old(i, neighbors+1, triang, z)
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