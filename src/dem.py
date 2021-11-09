import numpy as np
import math


class DigitalElevationModel:
    def __init__(self, filepath, headerpath):
        self.width = 0
        self.height = 0
        self.Xll = 0
        self.Yll = 0
        self.CellSize = 0
        self.NoData = 0
        self.filename = filepath
        self.header_path = headerpath
        self.s = None

    def read_dem(self):
        """Reads a raster DEM, returns it as numpy array"""
        try:
            file = open(self.header_path, "r")
            self.width = int(file.readline().rstrip('\n'))
            self.height = int(file.readline().rstrip('\n'))
            self.Xll = float(file.readline().rstrip('\n'))
            self.Yll = float(file.readline().rstrip('\n'))
            self.CellSize = float(file.readline().rstrip('\n'))
            self.NoData = float(file.readline().rstrip('\n'))
            file.close()
        except FileNotFoundError:
            print("Reading the hdr failed")

    def read_dem_to_array(self):
        data = np.fromfile(self.filename, dtype='h') * 0.01
        self.s = data.reshape(self.width, self.height)

    def get_height(self, x_kkj, y_kkj):
        x_utm = -2471441.562 + 0.9987798071 * x_kkj + 0.04612734592 * y_kkj
        y_utm = 124518.3273 - 0.04613846192 * x_kkj + 0.9987750048 * y_kkj
        dem_x = int(math.floor((x_utm - 355133.5) / 1))
        dem_y = int(math.floor(((y_utm - 6855143.5) / 1)))
        dem_z = self.s[int(dem_x)][int(dem_y)] - 18.67 + 0.32
        return dem_z
