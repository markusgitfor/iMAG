import numpy as np


class Point:
    def __init__(self):
        self.x = None
        self.y = None
        self.z = None


def vector_angle(vec1: Point, vec2: Point):
    return np.arccos((vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z) / (vector_length(vec1) *
                                                                              vector_length(vec2)))


def vector_length(vector: Point):
    return np.power((np.power(vector.x, 2) + np.power(vector.y, 2) + np.power(vector.z, 2)), 0.5)


class ImageHeader:
    def __init__(self, filename: str):
        self.filename = filename
        self.image_type = None
        self.image_code = None
        self.sub_width = None
        self.sub_height = None
        self.color = None
        self.c_col = None
        self.c_row = None
        self.o_col = None
        self.o_row = None
        self.width = None
        self.height = None
        self.constant = None
        self.x_ps = None
        self.y_ps = None
        self.lambda_ = None
        self.alpha = None
        self.mean_x = None
        self.mean_y = None
        self.X_mean = None
        self.Y_mean = None
        self.a_ = None
        self.b_ = None
        self.c_ = None
        self.d_ = None
        self.e_ = None
        self.f_ = None
        self.omega = None
        self.phi = None
        self.kappa = None
        self.x0 = None
        self.y0 = None
        self.z0 = None
        self.sun_azimuth = None
        self.sun_elevation = None
        self.start_of = None
        self.r = None

    def read_header(self):
        file = open(self.filename, mode='r')
        self.image_type = file.readline().rstrip('\n')
        self.image_code = int(file.readline().rstrip('\n'))
        self.sub_width = int(file.readline().rstrip('\n'))
        self.sub_height = int(file.readline().rstrip('\n'))
        self.filename = file.readline().rstrip('\n')
        self.color = int(file.readline().rstrip('\n'))
        self.c_col = int(file.readline().rstrip('\n'))
        self.c_row = int(file.readline().rstrip('\n'))
        self.o_col = int(file.readline().rstrip('\n'))
        self.o_row = int(file.readline().rstrip('\n'))
        self.width = int(file.readline().rstrip('\n'))
        self.height = int(file.readline().rstrip('\n'))
        self.constant = float(file.readline().rstrip('\n'))
        self.x_ps = float(file.readline().rstrip('\n'))
        self.y_ps = float(file.readline().rstrip('\n'))
        self.lambda_ = float(file.readline().rstrip('\n'))
        self.alpha = float(file.readline().rstrip('\n'))
        self.mean_x = float(file.readline().rstrip('\n'))
        self.mean_y = float(file.readline().rstrip('\n'))
        self.X_mean = float(file.readline().rstrip('\n'))
        self.Y_mean = float(file.readline().rstrip('\n'))
        self.a_ = float(file.readline().rstrip('\n'))
        self.b_ = float(file.readline().rstrip('\n'))
        self.c_ = float(file.readline().rstrip('\n'))
        self.d_ = float(file.readline().rstrip('\n'))
        self.e_ = float(file.readline().rstrip('\n'))
        self.f_ = float(file.readline().rstrip('\n'))
        self.omega = float(file.readline().rstrip('\n'))
        self.phi = float(file.readline().rstrip('\n'))
        self.kappa = float(file.readline().rstrip('\n'))
        self.x0 = float(file.readline().rstrip('\n'))
        self.y0 = float(file.readline().rstrip('\n'))
        self.z0 = float(file.readline().rstrip('\n'))
        self.sun_azimuth = float(file.readline().rstrip('\n'))
        self.sun_elevation = float(file.readline().rstrip('\n'))
        self.start_of = file.readline().rstrip('\n')
        self.r_transform_matrix()

    def get_filename(self):
        return self.filename

    def get_width(self):
        return self.width

    def get_height(self):
        return self.height

    def r_transform_matrix(self):
        """
        Returns camera rotation matrix given omega, phi, kappa rotations

        https://engineering.purdue.edu/~bethel/rot2.pdf
        """
        self.r = np.empty((3, 3))
        self.r[0][0] = np.cos(self.phi) * np.cos(self.kappa)
        self.r[0][1] = - np.cos(self.phi) * np.sin(self.kappa)
        self.r[0][2] = np.sin(self.phi)
        self.r[1][0] = np.cos(self.omega) * np.sin(self.kappa) + np.sin(self.omega) * np.sin(self.phi) *\
            np.cos(self.kappa)
        self.r[1][1] = np.cos(self.omega) * np.cos(self.kappa) - np.sin(self.omega) * np.sin(self.phi) *\
            np.sin(self.kappa)
        self.r[1][2] = -np.sin(self.omega) * np.cos(self.phi)
        self.r[2][0] = np.sin(self.omega) * np.sin(self.kappa) - np.cos(self.omega) * np.sin(self.phi) *\
            np.cos(self.kappa)
        self.r[2][1] = np.sin(self.omega) * np.cos(self.kappa) + np.cos(self.omega) * np.sin(self.phi) *\
            np.sin(self.kappa)
        self.r[2][2] = np.cos(self.omega) * np.cos(self.phi)

    def r_transform_ground_to_pixel(self, x, y, z):
        """
        Transform ground coordinates to image-coordinates

        :param i:
        :param x:
        :param y:
        :param z:
        :return:
        """
        camera_x, camera_y = self.r_transform_3d(x, y, z)
        direction = 0
        p_x, p_y = self.a_transform_affine(direction, camera_x, camera_y)
        return p_x, self.height - p_y

    def r_transform_3d(self, x, y, z):
        """
        This function calculates 3D -> 2D transformation

        :param i:
        :param x:
        :param y:
        :param z:
        :return:
        """
        camera_x = 0.0
        camera_y = 0.0
        k = (self.r[0][2] * (x - self.x0) + self.r[1][2] * (y - self.y0) + self.r[2][2] * (z - self.z0))
        if k != 0:
            camera_x = -self.constant * (self.r[0][0] * (x - self.x0) + self.r[1][0] * (y - self.y0) + self.r[2][0] *
                                         (z - self.z0)) / k
            camera_y = -self.constant * (self.r[0][1] * (x - self.x0) + self.r[1][1] * (y - self.y0) + self.r[2][1] *
                                         (z - self.z0)) / k
        return camera_x, camera_y

    def a_transform_affine(self, direction, x, y):
        """
        Return camera coordinates

        When direction is 0, transform is from  camera coordinates to image coordinates (pixels)
        When direction is 1, transform is from image coordinates (pixels) to camera coordinates (mm)

        :param direction:
        :param x:
        :param y:
        :return:
        """
        p_x = None
        p_y = None
        if direction == 0:
            p_x = self.a_ * x + self.b_ * y + self.c_
            p_y = self.d_ * x + self.e_ * y + self.f_
        if direction == 1:
            p_x = (-self.e_ * self.c_ + self.e_ * x + self.f_ * self.b_ - y * self.b_) / (self.a_ * self.e_ - self.b_ *
                                                                                          self.d_)
            p_y = -(self.a_ * self.f_ - self.a_ * y - self.c_ * self.d_ + x * self.d_) / (self.a_ * self.e_ - self.b_ *
                                                                                          self.d_)
        return p_x, p_y

    def view_illumination(self, xsol, ysol, zsol):
        """

        :return:
        """
        fii = np.radians(90 - np.degrees(self.sun_azimuth))
        theta = np.radians(90 - np.degrees(self.sun_elevation))

        sun_vector = Point()
        sun_vector.x = 1 * np.sin(theta) * np.cos(fii)
        sun_vector.y = 1 * np.sin(theta) * np.sin(fii)
        sun_vector.z = 1 * np.cos(theta)

        camera_vector = Point()
        camera_vector.x = self.x0 - xsol
        camera_vector.y = self.y0 - ysol
        camera_vector.z = self.z0 - zsol

        plumb_vector = Point()
        plumb_vector.x = 0
        plumb_vector.y = 0
        plumb_vector.z = 1

        t1 = Point()
        t1.x = sun_vector.x
        t1.y = sun_vector.y
        t1.z = 0

        t2 = Point()
        t2.x = camera_vector.x
        t2.y = camera_vector.y
        t2.z = 0

        sun_azimuth = np.degrees(fii)
        sun_elevation = 90 - np.degrees(theta)
        view_zenith = np.degrees(vector_angle(plumb_vector, camera_vector))
        azimuth_difference = np.degrees(vector_angle(t1, t2))
        view_azimuth = azimuth_difference - sun_azimuth
        return sun_azimuth, sun_elevation, view_zenith, azimuth_difference, view_azimuth
















