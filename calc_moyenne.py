import rasterio as rio
import numpy as np
import sys
import os

USAGE = 'calc_moyenne.py simdir'

def main():
    if len(sys.argv) < 2:
        print('ERROR. Not enough arguments.')
        print('Usage:', USAGE)
        sys.exit()

    rast_sum = None
    rast_sum2 = None

    simdir = sys.argv[1]
    files = os.listdir(simdir)
    cnt = 0
    for f in files:
        if f.endswith('.tif'):
            cnt += 1
            rast = read_raster(os.path.join(simdir, f))
            if rast_sum is None:
                rast_sum = rast
                rast_sum2 = np.power(np.array(rast), 2)
            else:
                rast_sum += rast
                rast_sum2 += np.power(rast, 2)

    moy = rast_sum / cnt
    stdev = np.sqrt((rast_sum2 - (np.power(rast_sum, 2) / cnt)) / cnt)

    new_dataset = rio.open(
        os.path.join(simdir, '_moyenne.tif'),
        'w',
        driver='GTiff',
        height=moy.shape[0],
        width=moy.shape[1],
        count=1, dtype=str(moy.dtype),
        crs=crs,
        transform=transform)
    new_dataset.write(moy, 1)
    new_dataset.close()

    new_dataset = rio.open(
        os.path.join(simdir, '_stdev.tif'),
        'w',
        driver='GTiff',
        height=moy.shape[0],
        width=moy.shape[1],
        count=1, dtype=str(moy.dtype),
        crs=crs,
        transform=transform)
    new_dataset.write(stdev, 1)
    new_dataset.close()

    print('Moyenne est dans fichier _moyenne.tif')
    print('Ecart-type est dans fichier _stdev.tif')

transform = None
crs = None

def read_raster(rast):
    global transform, crs
    ds = rio.open(rast)
    transform = ds.transform
    crs = ds.crs
    rast_arr = np.array(ds.read()[0], dtype=np.float64)
    return rast_arr


if __name__ == '__main__':
    main()
