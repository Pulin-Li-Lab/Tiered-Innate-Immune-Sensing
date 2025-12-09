"""
3D Segmentation and Read Assignment
Functions to perform 2.5D segmentation of DAPI stains and DAPI/amplicon channel overlays using StarDist
"""

# Import packages
from stardist.models import StarDist2D
from csbdeep.utils import normalize
from glob import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from tifffile import imread, imwrite
from skimage import filters, morphology, segmentation
from skimage.measure import regionprops, regionprops_table
from scipy import ndimage as ndi
from scipy.io import loadmat
import argparse
import timeit


#================================= LOAD DATA =================================

# 0. Load channels
def getChannels(input_path):
    """
    Load data 
        - load dapi channel 
        - load each of the 4 channels
        - take maximum of 4 channels and DAPI = overlay
    Return:
        - dapi channel (3D arrays)
        - overlay (3D array maximum of 5 channels)
    """
    # Get channel paths
    ch_paths = glob(os.path.join(input_path, 'registered/ch[0-3].tif'))
    dapi_path = glob(os.path.join(input_path, '*ch04.tif'))
    # Import and store images as 3D arrays
    chs = [imread(path) for path in ch_paths]
    dapi = imread(dapi_path)

    # Take maximum of all 4 channels and DAPI = overlay
    overlay = dapi
    for ch in chs:
        overlay = np.maximum(overlay, ch)

    return dapi, overlay

#================================= HELPER FNS =================================

# Log function
def log(msg, output_path, first=False):
    """
    Given a message, print message to screen for visual tracking and appends to log.txt
    @param msg: string to print and log
    @param output_path: path to tile folder
    @param first: boolean for whether to write (first log of the run, overwrites existing log) or append
    """
    print(msg)
    if first:
        fid = open(os.path.join(output_path, 'log.txt'), "w")
    else: 
        fid = open(os.path.join(output_path, 'log.txt'), "a")
    fid.write(msg)
    fid.close()

# Directory checker
def checkdir(output_path):
    """
    Check if an output directory exists, if not make it
    Used for output/segassign and output/segassign/tile_X directories
    """
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    

#================================= DAPI SEGMENTATION =================================

# 1. Perform 2D segmentation using StarDist2D pre-trained model
def segmentDapi2D(dapi, preprocess=False, save_path=None):
    """
    2D segmentation of nuclei
    returns a 2D label array
    """
    # Get StarDist2D pretrained versatile fluorescent nuclei model
    model = StarDist2D.from_pretrained('2D_versatile_fluo')

    # Take 2D sum projection of DAPI and segment with model
    dapi_sum = dapi.sum(axis=0)
    if preprocess:
        dapi_sum = filters.gaussian(dapi_sum, sigma=4)
    labels2D, _ = model.predict_instances(normalize(dapi_sum))

    if save_path:
        plt.figure(figsize=(20,10))
        plt.subplot(121)
        plt.title('DAPI')
        plt.imshow(dapi.max(axis=0))
        plt.subplot(122)
        plt.title('StarDist labels')
        plt.imshow(labels2D)
        plt.savefig(os.path.join(save_path, 'dapi_segmentation.png'))

    return labels2D

# 2. Process and binarize DAPI
def im_process(image, strel, show=False):
    """
    2D or 3D Image processing function 
    Performs sequential gaussian blur, median filter, thresholding, and morphological closing
    @param image is the image being processed
    @param strel is a structuring elem (2D or 3D depending on needs, ie ball(7) or disk(5))
    """
    
    blur = filters.gaussian(image, sigma=3)
    denoised = ndi.median_filter(blur, size=2)
    thresholded = denoised > filters.threshold_otsu(denoised)
    closed = morphology.closing(thresholded, selem=strel)
    
    if show == True:
        plt.figure(figsize=(10,10))
        plt.imshow(closed)
        
    return closed
    
# 3. 2.5D segment by multiplying each layer of binary threshold DAPI z-stack with StarDist 2D labels
def segmentDapi2_5D(dapi, labels2D, strel):
    """
    Perform 2.5D segmentation on 3D DAPI channel:
    - Extract each layer of DAPI z-stack
    - Process and binarize layer
    - Multiply binarized layer with 2D labels array to transfer labels to DAPI
    Return z-stack (3D array) of labeled DAPI channel
    """
    labels3D = np.zeros(dapi.shape)
    for z in range(0, dapi.shape[0]):
        layer = dapi[z] # extract
        processed = im_process(layer, strel, show=False)
        labels3D[z] = processed * labels2D

    return labels3D

# 4. Get centroid labels from labeled image (DAPI or cell) and use as watershed seeds for later use in overlay segmentation
def getCentroids(labels3D, save_df_path=None):
    """
    Use skimage regionprops to locate centroid of each labeled nuclei
    Returns a nx3 array, where n is the total number of nuclei labeled, and columns are z,y,x coords respectively
    Returns a df with centroid coordinates, area/volume, and corresponding region label
    """
    centroids = []
    labels3D = labels3D.astype(int)
    for i, region in enumerate(regionprops(labels3D)):
        centroids.append(region.centroid)
    centroids = np.array(centroids).astype(int)

    centroids_df = pd.DataFrame(
        regionprops_table(labels3D.astype(int),
                          properties = ('label', 'centroid', 'area'))
    )
    centroids_df.columns = ['cell_barcode', 'z', 'y', 'x', 'volume']
    centroids_df = centroids_df.astype({
        'z':int, 
        'y':int, 
        'x':int
    })
    if save_df_path:
        centroids_df.to_csv(save_df_path)

    return centroids



#================================= CELL SEGMENTATION =================================

# 5. Assign DAPI centroid as marker in shape of image
def getMarkerArray(overlay, centroids):
    """
    Returns an array of markers in the shape of the overlay
    Each corresponding labeled nuclei in labels3D has a single label at its centroid location, or marker
    """
    numCells = centroids.shape[0]
    markers = np.zeros(overlay.shape, dtype=np.uint8)
    for i in range(numCells):
        z,y,x = centroids[i,:]
        if z < overlay.shape[0] and y < overlay.shape[1] and x < overlay.shape[2]:
            markers[z,y,x] = i+1

    return markers, numCells

# 6. Process overlay image and watershed segment
def segmentCells3D(overlay, markers, strel3D):
    """
    Perform seed-based watershed segmentation of DAPI/amplicon overlay using markers
    @param strel is a 3D structuring element defining dilation parameter (i.e. ball)
    Returns a dilated 3D labeled image of cells
    """
    overlay_processed = im_process(overlay, strel3D)
    cellLabels = segmentation.watershed(overlay_processed, markers, mask=overlay_processed)
    cellLabels = morphology.dilation(cellLabels, selem=strel3D)
    
    return cellLabels

# 6.1 If no overlay image or DAPI signals are extremely dense, just dilate dapi labels to cell labels
def dilateCells3D(labels3D, strel3D):
    return morphology.dilation(labels3D, strel3D)

# 7. Filter by cell size
def filter_cellSize(labels3D, size_thres=0):
    if size_thres == 0:
        return
    regions = regionprops(labels3D.astype(int))
    for i in regions:
        if i.area < size_thres:
            labels3D[labels3D==i.label] = 0
    return labels3D


#================================= READ ASSIGNMENT =================================

# 8. Load reads.csv (sitk output)
def loadReads(base_path, output_path, output_tile_path, tile_num):
    """
    Loads and returns a file containing amplicon read coordinates and genes
    """
    reads = pd.read_csv(os.path.join(base_path, f'output/tile_{tile_num}.csv'), index_col=0)
    #log(f"\tNumber of reads: {len(reads)}\n", output_tile_path)

    return reads

# 8.1 Load goodPoints_max3D.mat (matlab output)
def loadReadsMat(base_path, tile_num, output_path, assign_dir):
    dots = loadmat(base_path, f'output/max/tile_{tile_num}')
    bases = [i[0] for i in dots["goodReads"]]
    bases = np.array(bases)
    temp = dots["goodSpots"]
    temp = temp[:,:3]
    points = np.zeros(temp.shape)
    points[:,0] = temp[:,0]-1
    points[:,1] = temp[:,1]-1
    points[:,2] = temp[:,2]-1

    # Log
    assign_output_path = os.path.join(output_path, assign_dir)
    if not os.path.exists(assign_output_path):
        os.makedirs(assign_output_path)
    #log(f"\tNumber of reads: {len(bases)}\n", output_tile_path)
    
    return bases, points

# 8.2 Load goodPoints_max3d.csv (output of 3-part registration pipeline)
def loadReadsNew(base_path, tile_num):
    # Load reads (currently not 0-based, so adjust)
    reads = pd.read_csv(os.path.join(base_path, f'02.processed_data/Position{tile_num}/goodPoints_max3d.csv'))
    for c in ['x', 'y', 'z']:
        reads[c] -= 1

    return reads


# 9. Load genes.csv (with the raw data)
def loadGenes(base_path):
    """
    Loads genes.csv located in root data directory
    """
    gene_path = os.path.join(base_path, '01.data/genes.csv') 
    genes = pd.read_csv(gene_path, header=None, names=["Gene Name", "Barcode"])["Gene Name"]
    gene2idx = {gene:idx for idx, gene in enumerate(genes)}
    return genes, gene2idx

# 10. Assign reads and write to file
def assignReads(numCells, reads, genes, cellLabels, gene2idx, output_tile_path):
    matrix = np.zeros((numCells, len(genes)))
    for i in range(len(reads)):
        x,y,z,gene = reads.iloc[i].values
        cell = cellLabels[z,y,x]
        if cell != 0:
            matrix[cell.astype(int)-1, gene2idx[gene]] += 1

    # Write matrix to file
    np.savetxt(os.path.join(output_tile_path, 'readMatrix.csv'), matrix, fmt='%d', delimiter=',')

    return matrix # just for logging purposes

# 10.1 Assign reads and write to file (for MATLAB input)
def assignReadsMat(numCells, bases, points, cellLabels, gene2idx, output_path, assign_dir, tile_num):
    matrix = np.zeros((numCells, len(gene2idx.keys())))
    genecounts_mat = {}
    for i in range(len(points)):
        x,y,z = points[i,:].astype(int)
        gene = ''.join(bases[i])
        cell = cellLabels[z,y,x]
        if cell != 0:
            matrix[cell.astype(int)-1, gene2idx[gene]] += 1
            genecounts_mat[gene] += 1
    matrix = matrix.astype(int)

    # Write matrix to file
    np.savetxt(os.path.join(output_path, assign_dir, f'tile{tile_num}_readmatrix_mat.csv'))

    return matrix

# 10.2 Assign reads and create a remain_reads.csv file 
def assignRemainReads(reads, cellLabels, output_tile_path):
    cell_barcode = []
    for i in range(len(reads)):
        x,y,z,gene = reads.iloc[i].values
        cell = cellLabels[z,y,x] # extract whether read overlaps with a cell
        if cell != 0:
            cell_barcode.append(cell) # append corresponding cell label
        else:
            cell_barcode.append(0)
    reads['cell_barcode'] = cell_barcode

    # Save assigned reads
    remain_reads = reads[reads['cell_barcode']!=0]
    remain_reads.to_csv(os.path.join(output_tile_path, 'remain_reads.csv'))

    log(f"Total reads: {len(cell_barcode)}\n", output_tile_path)
    log(f"Assigned reads: {cell_barcode.count(0)}\n", output_tile_path)
    pct_assigned = round((len(cell_barcode) - cell_barcode.count(0)) / len(cell_barcode) * 100, 2)
    log(f"Assignment percentage: {pct_assigned}%\n", output_tile_path)
    
    return remain_reads, pct_assigned

# 11. Plot results
def plotSegAssignResults(dapi, cell_centroids, cellLabels, reads, remain_reads, output_tile_path):
    fig = plt.figure(figsize=(30,30))
    fig.subplots_adjust(hspace=0.05, wspace=0.05)

    # Plot dapi sum projection
    plt.subplot(221)
    plt.title('DAPI sum projection', fontdict={'fontsize':30})
    plt.axis('off')
    plt.imshow(dapi.sum(axis=0))

    # Plot cell segmentation 2D max projection
    plt.subplot(222)
    plt.title(f'Cell segmentation ({len(cell_centroids)} cells)', fontdict={'fontsize':30})
    plt.axis('off')
    plt.imshow(cellLabels.max(axis=0))

    # Plot dapi with cell centers
    plt.subplot(223)
    plt.title('DAPI with cell centers', fontdict={'fontsize':30})
    plt.axis('off')
    plt.imshow(dapi.sum(axis=0))
    plt.scatter(cell_centroids[:,2], cell_centroids[:,1], s=25, c='red')

    # Plot dapi with remain reads and all reads
    plt.subplot(224)
    pct_assigned = round((len(reads) - len(remain_reads) / len(reads) * 100), 2)
    plt.title(f'DAPI with all/assigned reads ({pct_assigned}%)', fontdict={'fontsize':30})
    plt.axis('off')
    plt.imshow(dapi.sum(axis=0))
    plt.scatter(reads['x'], reads['y'], alpha=0.1, c='red', s=10)
    plt.scatter(remain_reads['x'],remain_reads['y'],s = 10,alpha = 0.8,c=pd.Categorical(np.array(remain_reads['cell_barcode'])).codes, cmap= matplotlib.colors.ListedColormap ( np.random.rand ( 256,3)))

    plt.savefig(os.path.join(output_tile_path, 'segmentation_assignment_results.png'))


#================================= MASTER FUNCTIONS =================================

# Segmentation 
def segment(input_path, output_tile_path, strel2D, strel3D):
    """
    Wrapping together all segmentation functions
    @param input_path: base_path/round1/tile_num
    @param output_tile_path: base_path/output/segassign/tile_X
    @param strel2D: structuring element used in morphological closing for 2.5D DAPI segmentation
    @param strel3D: structuring element used in morphological closing and dilation for 3D overlay segmentation
    #@param assign: boolean for whether assignment will be performed in this job
    """
    start = timeit.default_timer()

    # DAPI segmentation
    log("Performing DAPI Segmentation...", output_tile_path)
    dapi, overlay = getChannels(input_path)
    labels2D = segmentDapi2D(dapi)
    labels3D = segmentDapi2_5D(dapi, labels2D, strel2D)
    centroids = getDapiCentroids(labels3D)

    # Cell segmentation
    log("Performing cell segmentation...", output_tile_path)
    markers, numCells = getMarkerArray(overlay, centroids)
    cellLabels = segmentCells3D(overlay, markers, strel3D)
    saveCellLabels(cellLabels, output_tile_path)
    
    stop = timeit.default_timer()

    # Log
    msg = f"Segmentation results:\n\tNumber of cells: {numCells}\n\tTime: {round(stop-start,3)} seconds\n"
    log(msg, output_tile_path)


# Read assignment
def readAssignment(base_path, output_path, output_tile_path, tile_num, mat):
    """
    Wrapping together all read assignment functions 
    @param output_tile_path: base_path/output/segassign/tile_X
    @param tile_num: number of tile
    @param mat: boolean for whether reads locations are output from MATLAB
    @param numCells: number of cells labeled in segmentation mask
    @param seg: boolean for whether or not seg was performed in this job
    @param seg_dir: base_path/output/segassign/seg
    """
    start = timeit.default_timer()

    log("Performing read assignment...", output_tile_path)
    
    # Load segmentation mask
    if os.path.exists(output_tile_path): 
        cellLabels = imread(os.path.join(output_tile_path, 'cellLabels.tiff'))
        numCells = np.unique(cellLabels).max().astype(int)
    else:
        raise ValueError("Cannot perform read assignment without segmentation mask")
        
    if not mat:
        reads = loadReads(base_path, output_path, output_tile_path, tile_num)
        genes, gene2idx = loadGenes(base_path)
        matrix = assignReads(numCells, reads, genes, cellLabels, gene2idx, output_tile_path)
    else:
        bases, points = loadReadsMat(base_path, tile_num, output_path)
        _, gene2idx = loadGenes(base_path)
        matrix = assignReadsMat(numCells, bases, points, cellLabels, gene2idx, output_path, tile_num)

    stop = timeit.default_timer()

    # Log
    msg = f"Read assignment result:\n\t{round(matrix.sum()/len(reads) * 100, 4)}% reads assigned ({matrix.sum().astype(int)} out of {len(reads)} reads)\n\tTime: {round(stop-start, 3)} seconds\n"
    log(msg, output_tile_path)

# Segmentation and assignment v2 (using results from 3-part pipeline):
def segAssign3D(base_path, tile_num, input_path, output_tile_path, strel2D, strel3D):
    start = timeit.default_timer()

    # DAPI segmentation
    log("Performing DAPI Segmentation...\n", output_tile_path)
    dapi, _ = getChannels(input_path)
    labels2D = segmentDapi2D(dapi, preprocess=True, save_path=os.path.join(output_tile_path))
    labels3D = segmentDapi2_5D(dapi, labels2D, strel2D)
    dapi_centroids = getCentroids(labels3D, save_df_path=os.path.join(output_tile_path, 'dapi_centroids.csv'))

    # Cell segmentation
    log("Performing cell segmentation...\n", output_tile_path)
    cellLabels = dilateCells3D(labels3D, strel3D)
    cellLabels = filter_cellSize(cellLabels, size_thres=500)
    cell_centroids = getCentroids(cellLabels, save_df_path=os.path.join(output_tile_path, 'cell_center.csv'))

    # Save cell labels and dapi_segmentation
    np.save(os.path.join(output_tile_path, 'cellLabels.npy'), cellLabels)
    np.save(os.path.join(output_tile_path, 'dapiLabels.npy'), labels3D)  

    stop = timeit.default_timer()

    # Log
    msg = f"Segmentation results:\n\tNumber of cells: {cell_centroids.shape[0]}\n\tTime: {round(stop-start,3)} seconds\n"
    log(msg, output_tile_path)

    # Reads assignment
    start = timeit.default_timer()

    log("Performing read assignment...\n", output_tile_path)
    
    reads = loadReadsNew(base_path, tile_num)
    remain_reads, pct_assigned = assignRemainReads(reads, cellLabels, output_tile_path)
    plotSegAssignResults(dapi, cell_centroids, cellLabels, reads, remain_reads, output_tile_path)

    stop = timeit.default_timer()

    # Log
    msg = f"Read assignment result:\n\t{pct_assigned}% reads assigned ({len(remain_reads)} out of {len(reads)} reads)\n\tTime: {round(stop-start, 3)} seconds\n"
    log(msg, output_tile_path)

    print("Segmentation and reads assignment done!")

    
    

#================================= SHELL EXECUTION =================================

def main(args):

    # Parse arguments
    base_path = args.base_path
    tile_num = str(args.tile_num).zfill(3)
    #seg = args.seg
    #assign = args.assign
    mat = args.mat
    strel2D = args.strel2D
    strel3D = args.strel3D
    input_path = os.path.join(base_path, f'01.data/round1/Position{tile_num}')
    output_path = os.path.join(base_path, f'03.segmentation/watershed_sum')
    os.makedirs(output_path, exist_ok=True)
    output_tile_path = os.path.join(output_path, f'Position{tile_num}')
    os.makedirs(output_tile_path, exist_ok=True)

    # Check if tile exists
    if not os.path.exists(input_path):
        raise NameError('Tile data does not exist or is not accessible')

    log(f"=================================== TILE {tile_num} ===================================\n", output_tile_path, first=True)

    # Segment
    #segment(input_path, output_tile_path, strel2D, strel3D)

    # Assign Reads
    # readAssignment(base_path, output_path, output_tile_path, tile_num, mat)

    segAssign3D(base_path, tile_num, input_path, output_tile_path, strel2D, strel3D)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Segment and Assign Reads from registered STARmap images')
    parser.add_argument('base_path', type=str, help='path to data directory')
    parser.add_argument('tile_num', type=int, help='number of tile being processed')
    #parser.add_argument('-seg', action='store_true', help='flag for performing cell segmentation and creating a mask')
    #parser.add_argument('-assign', action='store_true', help='flag for performing read assignment')
    parser.add_argument('-mat', action='store_true', help='flag for if reads locations are output by MATLAB')
    parser.add_argument('--strel2D', type=np.ndarray, default=morphology.disk(7), help='structuring element for 2.5D morphological operations')
    parser.add_argument('--strel3D', type=np.ndarray, default=morphology.ball(7), help='structuring element for 3D morphological operations for overlay processing and label dilation')

    args = parser.parse_args()

    main(args)















