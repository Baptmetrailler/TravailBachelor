#from osgeo import gdal
import rasterio as rio
import numpy as np
import sys
import os
import random
import math
import time
import quadtreed3 as qtree


USAGE = 'sim.py raster_emplois raster_log nsimulations outdir'




def main():
    if len(sys.argv) < 5:
        print('ERROR. Not enough arguments.')
        print('Usage:', USAGE)
        sys.exit()

    f_empl = sys.argv[1]
    f_log = sys.argv[2]
    n = int(sys.argv[3])
    d_out = sys.argv[4]

    if not os.path.exists(f_empl):
        print('ERROR. File "%s" not found.' % f_empl)
        sys.exit()

    if not os.path.exists(f_log):
        print('ERROR. File "%s" not found.' % f_log)
        sys.exit()

    if n < 0 or n > 1000000:
        print('ERROR. nsim must be in range 1 to 1000000')
        sys.exit()

    os.makedirs(d_out, exist_ok=True)

    empl = read_raster(f_empl)
    log = read_raster(f_log)

    ds = rio.open(f_log)
    transform = ds.transform
    crs = ds.crs
    t0 = time.time()
    for i in range(n):
        run(empl, log, d_out, i+1, transform, crs)
        print('t', time.time() - t0)
        t0 = time.time()


def read_raster(rast):
    #ds = gdal.Open(rast)
    #band = ds.GetRasterBand(1)
    #rast_arr = np.array(band.ReadAsArray(), dtype=np.int32)
    #ds = None
    ds = rio.open(rast)
    rast_arr = np.array(ds.read()[0], dtype=np.int32)

    return rast_arr


def run(empl, log, d_out, cnt, transform, crs):
    if log.shape[0] != empl.shape[0] or log.shape[1] != empl.shape[1]:
        print('ERROR. Size of emplois and logement not identical.')
        sys.exit()



    m = log.shape[0]
    n = log.shape[1]
    empl_liste = rast_to_list(empl)
    distances_moyennes = np.zeros((m, n), dtype=np.float32)
    n_empl = len(empl_liste)
    print('Number of emplois found:', n_empl)

    # Random sort empl_liste
    np.random.shuffle(empl_liste)
    log_liste = rast_to_list(log)
    n_log = len(log_liste)
    print('Number of logements found:', n_log)

    log_to_add = []
    for p in range(max(0, n_empl - n_log)):
        log_to_add.append(random.choice(log_liste))
    log_liste += log_to_add
    print('Number of logements after correction:', len(log_liste))

    #update the log array
    updated_log = np.zeros((m, n), dtype=np.int32)

    for logement in log_liste:
        updated_log[logement[0], logement[1]] += 1

    print(sum(sum(updated_log)))
    #updated_log_orig = [k for k in updated_log]
    log = np.array(updated_log)
    #dist = compute_dist(empl_liste, log_liste)
    #k = 1

    sum_distances = np.zeros((m, n), dtype=np.float64)

    for emploi in empl_liste:
        occupied = True
        r = 0
        potential_log = []
        i,j = emploi
        while occupied:
            r_imin, r_imax = max(0, i - r), i + r + 1
            r_jmin, r_jmax = max(0, j - r), j + r + 1
            if np.sum(log[r_imin:r_imax, r_jmin:r_jmax]) > 0:
                potential_cells = np.nonzero(log[r_imin:r_imax, r_jmin:r_jmax])
                c = random.randint(0, len(potential_cells[0]) - 1)
                cell = (potential_cells[0][c], potential_cells[1][c])
                p, q = (r_imin + cell[0], r_jmin + cell[1]) # raster coordinates du logement
                potential_log.append((p, q))
                log[p, q] -= 1
                d = math.sqrt((i - p)**2 + (j - q)**2)
                sum_distances[p, q] += d
                occupied = False
            else:
                r += 1

    # Divide sum_distances by nemplois pour obtenir la moyenne
    distances_moyennes = np.where(updated_log == 0, 0, sum_distances / np.array(updated_log))
    
    name = d_out + "/dist_moy_"+str(cnt).zfill(4)+".tif"

    new_dataset = rio.open(name, 'w', driver='GTiff',
                                height=distances_moyennes.shape[0], width=distances_moyennes.shape[1],
                                count=1, dtype=str(distances_moyennes.dtype),
                                crs=crs,
                                transform=transform)

    new_dataset.write(distances_moyennes, 1)
    new_dataset.close()
    print(cnt)









    # Créer un array avec les distances, et un autre avec l'écart-type
    # Faire une boucle à travers dist pour calculer la distance moyenne et l'écart-type
    # Ecrire arrays dans dossier de sortie, e.g. dist_moy_0001.tif, dist_std_0001.tif

def dist (a,b):
    return math.sqrt((a[0] - b[0])** 2 + (a[1] - b[1])** 2)



def compute_dist(emplois, log):
    dist = {}
    for emp in emplois:
        # Chercher le logement le plus proche
        idx = min_dist(emp[0], emp[1], log)
        i, j = log.pop(idx)
        d = dist.get('%i/%i' % emp, [])
        d.append(math.sqrt((emp[0] - i)**2 + (emp[1] - j)**2))
        dist['%i/%i' % emp] = d
    return dist


def min_dist(i, j, liste):
    # Chercher le voisin le plus proche...
    min_d, min_k = np.Infinity, None
    for k in range(len(liste)):
        x,y = liste[k]
        d = (x - i)**2 + (y - j)**2
        if d < min_d:
            min_d = d
            min_k = k
    return min_k

def rast_to_list(rast):
    h, w = rast.shape
    liste = []
    for i in range(h):
        for j in range(w):
            n = rast[i, j]
            for k in range(n):
                liste.append((i, j))
    return liste


if __name__ == '__main__':
    main()
