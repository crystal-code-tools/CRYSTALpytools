#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The module for 3D / 2D elastic properties.
"""
from CRYSTALpytools import units
import numpy as np

__all__ = ['Tensor3D', 'Tensor2D',
           'young', 'comp', 'shear', 'poisson',
           'tensor_from_file', 'voigt2cart', 'cart2voigt']


class Tensor3D():
    """
    3D elastic tensor and related properties.

    Args:
        matrix (numpy.ndarray): 6\*6 compliance or stiffness matrix. Unit: GPa :math:`^{-1}` or GPa.
        lattice (numpy.ndarray): lattice matrix.
        is_compliance (bool): Compliance or stiffness matrix used as input.
    """
    def __init__(self, matrix, lattice=None, is_compliance=True):
        import warnings
        from pymatgen.core.lattice import Lattice

        if np.all(lattice==None):
            warnings.warn("Lattice information not available, Miller indices and fractional coordinates disabled.",
                          stacklevel=2)
            self.lattice = None
        else:
            self.lattice = Lattice(np.array(lattice, dtype=float))

        if is_compliance == True:
            self._S = np.array(matrix, dtype=float)
            self._C = np.linalg.inv(matrix)
        else:
            self._C = np.array(matrix, dtype=float)
            self._S = np.linalg.inv(matrix)

        # convert matrices from Voigt notation to Cartesian notation
        self._Ccart = voigt2cart(self._C)
        self._Scart = voigt2cart(self._S)

    @property
    def stiffness(self):
        """
        6\*6 stiffness matrix in Voigt annotation. Unit: GPa.
        """
        return self._C

    @property
    def compliance(self):
        """
        6\*6 compliance matrix in Voigt annotation. Unit: GPa :math:`^{-1}`.
        """
        return self._S

    @classmethod
    def from_file(cls, output, conventional_lattice=True):
        """
        Read elastic tensor from CRYSTAL output file and generate ``Tensor3D``
        object. Calls the ``crystal_io.Crystal_output.get_elatensor()`` method.
        Lattice information is obtained.

        Args:
            output (str): CRYSTAL output file.
            conventional_lattice (bool): Get conventional lattice.

        Returns:
            cls (Tensor3D)
        """
        from CRYSTALpytools.crystal_io import Crystal_output
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        out = Crystal_output(output)
        if out.get_dimensionality() != 3:
            raise Exception("Dimensionality error. Input system is not a 3D system.")

        if conventional_lattice == True:
            struc = SpacegroupAnalyzer(out.get_geometry(initial=False)).get_conventional_standard_structure()
            latt = struc.lattice.matrix
        else:
            latt = out.get_lattice(initial=False)
        return cls(matrix=out.get_elatensor(), lattice=latt, is_compliance=False)

    def get_1D(self, property, u, nchi=180, use_cartesian=True):
        """
        Compute properties along the given vector. Available properties are
        defined by the ``property`` variable. To get shear modulus / Poisson
        ratio of 2 specified vectors, use ``get_pair``.

        Options:

        * "young": Young's modulus in GPa.
        * "comp": Linear compressibility.
        * "shear": Shear modulus on the 'loop' normal to vector, in GPa.
        * "shear avg": Averaged shear modulus on the 'loop', in GPa.
        * "shear min": Minimum shear modulus on the 'loop', in GPa.
        * "shear max": Maximum shear modulus on the 'loop', in GPa.
        * "poisson": Poisson ratio on the 'loop' normal to vector.
        * "poisson avg": Averaged Poisson ratio on the 'loop'.
        * "poisson min": Minimum Poisson ratio on the 'loop'.
        * "poisson max": Maximum Poisson ratio on the 'loop'.

        Args:
            property (str): Property to compute. See options above.
            u (numpy.ndarray): 3\*1 or nu\*3 array of vectors.
            nchi (int): Resolution of auxiliary angle :math:`\\chi`, in radian
                :math:`[0, 2\\pi)`. Shear modulus and Poisson ratio only.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)

        Returns:
            e (float | numpy.ndarray): float or nu\*1, some sort of property.
        """
        property_list = ['young', 'comp', 'shear', 'shear avg', 'shear min',
                         'shear max', 'poisson', 'poisson avg', 'poisson min',
                         'poisson max']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        if use_cartesian != True and np.all(self.lattice==None):
            raise Exception('Lattice matrix not available.')

        u = np.array(u, dtype=float, ndmin=2)

        if use_cartesian != True:
            for i in range(len(u)):
                u[i, :] = self.lattice.get_cartesian_coords(u[i, :])
                u[i, :] = u[i, :] / np.linalg.norm(u[i, :])
        else:
            for i in range(len(u)):
                u[i, :] = u[i, :] / np.linalg.norm(u[i, :])
        # define the loop around u
        if 'shear' in property or 'poisson' in property:
            v = np.zeros([len(u), nchi+1, 3], dtype=float)
            theta = np.arccos(u[:,2])
            phi = np.arctan2(u[:,1], u[:,0])
            chi = np.linspace(0, 2*np.pi, nchi+1)
            for i in range(len(u)):
                v[i, :, :] = np.vstack([
                    np.cos(phi[i]) * np.cos(theta[i]) * np.cos(chi) - np.sin(phi[i]) * np.sin(chi),
                    np.sin(phi[i]) * np.cos(theta[i]) * np.cos(chi) + np.cos(phi[i]) * np.sin(chi),
                    -np.sin(theta[i]) * np.cos(chi)
                ]).transpose()

        if property == 'young':
            e = young(self._Scart, u, 3)
        elif property == 'comp':
            e = comp(self._Scart, u, 3)
        elif property == 'shear':
            e = shear(self._Scart, u, v, 3)
        elif property == 'shear avg':
            e = np.average(shear(self._Scart, u, v, 3), axis=1)
        elif property == 'shear min':
            e = np.min(shear(self._Scart, u, v, 3), axis=1)
        elif property == 'shear max':
            e = np.max(shear(self._Scart, u, v, 3), axis=1)
        elif property == 'poisson':
            e = poisson(self._Scart, u, v, 3)
        elif property == 'poisson avg':
            e = np.average(poisson(self._Scart, u, v, 3), axis=1)
        elif property == 'poisson min':
            e = np.min(poisson(self._Scart, u, v, 3), axis=1)
        elif property == 'poisson max':
            e = np.max(poisson(self._Scart, u, v, 3), axis=1)

        if len(u) == 1:
            e = e[0]
        return e

    def get_pair(self, property, u, v, use_cartesian=True):
        """
        Compute properties between given pairs of vectors. Used for shear
        modulus and Poisson ratio. For properties along the normal vector, use
        ``get_1D``.

        Options:

        * "shear": Shear modulus between ``u`` and ``v``, in GPa.
        * "poisson": Poisson ratio between ``u`` and ``v``.

        Args:
            property (str): Property to compute. See options above.
            u (numpy.ndarray): 3\*1 or nu\*3 array of vector 1.
            v (numpy.ndarray): 3\*1 or nu\*3 array of vector 2.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)

        Returns:
            e (float | numpy.ndarray): float or nu\*1, some sort of property.
        """
        property_list = ['shear', 'poisson']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        if use_cartesian != True and np.all(self.lattice==None):
            raise Exception('Lattice matrix not available.')

        u = np.array(u, dtype=float, ndmin=2)
        v = np.array(v, dtype=float, ndmin=2)
        if u.shape[0] != v.shape[0]:
            raise ValueError("Inconsistent length of u and v.")

        if use_cartesian != True:
            for i in range(np.shape(u)[0]):
                u[i, :] = self.lattice.get_cartesian_coords(u[i, :])
                v[i, :] = self.lattice.get_cartesian_coords(v[i, :])
                u[i, :] = u[i, :] / np.linalg.norm(u[i, :])
                v[i, :] = v[i, :] / np.linalg.norm(v[i, :])
        else:
            for i in range(len(u)):
                u[i, :] = u[i, :] / np.linalg.norm(u[i, :])
                v[i, :] = v[i, :] / np.linalg.norm(v[i, :])
        # orthogonality
        for i in range(len(u)):
            if np.abs(np.dot(u[i, :], v[i, :])) > 1e-6:
                raise ValueError("Vector pair {:d} (from 0) is non-orthogonal.".format(i))

        # add 1 extra dimension to v (nu*1*3)
        v = np.expand_dims(v, axis=1)
        if property == 'shear':
            e = shear(self._Scart, u, v, 3)
        elif property == 'poisson':
            e = poisson(self._Scart, u, v, 3)

        if len(u) == 1:
            e = e[0, 0]
        else:
            e = e[:, 0]
        return e

    def plot_2D(
        self, property, plane, ntheta=90, nchi=180, plane_definition='miller',
        u=None, utext=None, use_cartesian=True, plot_lattice=True,
        loop_color=None, loop_linestyle=None, loop_linewidth=None,
        figsize=[6.4, 4.8], get_data=False, **kwargs):
        """
        Plot 2D crystal elastic properties across the given plane. The plane is
        defined by any of the following methods:

        * 'miller': Miller indices.  
        * 'cartesian': The plane normal vector in Cartesian coordinates.  
        * 'fractional': The plane normal vector in fractional coordinates.

        Options:

        * "young": Young's modulus in GPa.  
        * "comp": Linear compressibility.  
        * "shear": Shear modulus in the plane, in GPa.  
        * "shear avg": Averaged shear modulus of vectors in the plane, in GPa.  
        * "shear min": Minimum shear modulus of vectors in the plane, in GPa.  
        * "shear max": Maximum shear modulus of vectors in the plane, in GPa.  
        * "poisson": Poisson ratio in the plane.  
        * "poisson avg": Averaged Poisson ratio of vectors in the plane.  
        * "poisson min": Minimum Poisson ratio of vectors in the plane.  
        * "poisson max": Maximum Poisson ratio of vectors in the plane.

        This method can be used for multi-plane, multi-property plots. The same
        scale is enforced for subplots of the same property. Subplot titles are
        added automatically. The layout of subplots is, by default,
        nProp\*nPlane or 1\*nProp, which is not editable.

        Args:
            property (str|list): Property(ies) to compute. See options above.
            plane (numpy.ndarray): 3\*1 or nplane\*3 array of planes.
            ntheta (int): Resolution of azimuth angle :math:`\\theta` on plane,
                in radian :math:`[0, 2\\pi)`.
            nchi (int): Resolution of auxiliary angle  :math:`\\chi`, in radian
                :math:`[0, 2\\pi)`. Shear modulus and Poisson ratio only.
            plane_definition (str): The method used to define the plane. See
                options above.
            u (numpy.ndarray): 3\*1 or nu\*3 array of vectors to be added into
                the figure. A line segment with doubled radius is plotted to
                indicate the vector. Or 'max' / 'min' / 'bothends', plot
                vectors corresponding to the max and min of elastic properties.
                For 'bothends', min first.
            utext (list[str] | str): A string or a list of string. Used to mark
                The vector. If ``None`` or 'value', the value is annotated.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)
            plot_lattice (bool): Draw lattice base vectors to indicate
                orientation of the plane.
            loop_color (str|list): Color of 2D loops. If only one string is
                given, apply it to all loops. 'None' for default values
                (matplotlib Tableau Palette).
            loop_linestyle (str|list): Linestyle of 2D loops. If only one
                string is given, apply it to all plots. 'None' for default
                values ('-').
            loop_linewidth (str|list): Linewidth of band structure. If only one
                number is given, apply it to all plots. For the 'multi' mode,
                1\*nSystem list of linewidth numbers. 'None' for default values
                (1.0).
            figsize (list): The figure size specified as \[width, height\].
            get_data (bool): *Developer only* Return data or figure.
            \*\*kwargs: Parameters passed to ``PolarAxes.plot`` to customize
                the loops. Applied to all the loops.

        Returns:
            fig (Figure): Matplotlib ``Figure`` object.
        """
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        import warnings

        # sanity check
        property_list = ['young', 'comp', 'shear', 'shear avg', 'shear min',
                         'shear max', 'poisson', 'poisson avg', 'poisson min',
                         'poisson max']
        self.property = np.array(property, ndmin=1)
        nprop = len(self.property)
        for i in range(nprop):
            if self.property[i].lower() not in property_list:
                raise ValueError("Unknown property input: '{}'".format(self.property[i]))
            self.property[i] = self.property[i].lower()

        planedef_list = ['miller', 'cartesian', 'fractional']
        plane_definition = plane_definition.lower()
        if plane_definition not in planedef_list:
            raise ValueError("Unknown plane definition: '{}'".format(plane_definition))

        if np.all(self.lattice==None):
            if use_cartesian != True or plot_lattice == True:
                raise Exception('Lattice matrix not available.')

        # plane definition: use normalised cartesian for plot
        plane = np.array(plane, ndmin=2); self.norm = np.array(plane, dtype=float)
        if plane_definition == 'miller':
            for i in range(len(self.norm)):
                self.norm[i, :] = self.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(self.norm[i, :])
                self.norm[i, :] = self.norm[i, :] / np.linalg.norm(self.norm[i, :])
        elif plane_definition == 'fractional':
            for i in range(len(self.norm)):
                self.norm[i, :] = self.lattice.get_cartesian_coords(self.norm[i, :])
                self.norm[i, :] = self.norm[i, :] / np.linalg.norm(self.norm[i, :])
        else:
            for i in range(len(self.norm)):
                self.norm[i, :] = self.norm[i, :] / np.linalg.norm(self.norm[i, :])

        nplane = len(self.norm)

        # Compute each 2D plane separately
        self.r = np.zeros([nprop, nplane, ntheta+1], dtype=float)
        chi = np.linspace(0, 2*np.pi, ntheta+1) # pole angle is actually chi.
        self.chi = chi

        for iplot, v0 in enumerate(self.norm):
            # get pole angles and 3D vectors. use (v0) v1. Data is always associated with v1
            theta = np.arccos(v0[2])
            phi = np.arctan2(v0[1], v0[0])
            v1 = np.vstack([
                np.cos(phi) * np.cos(theta) * np.cos(chi) - np.sin(phi) * np.sin(chi),
                np.sin(phi) * np.cos(theta) * np.cos(chi) + np.cos(phi) * np.sin(chi),
                -np.sin(theta) * np.cos(chi)
            ]).transpose()
            for ip, p in enumerate(self.property):
                # get pole angles and 3D vectors. use v1, v2. Data is always associated with v1
                if 'avg' in p or 'min' in p or 'max' in p:
                    v2 = np.zeros([ntheta+1, nchi+1, 3], dtype=float) # another chi.
                    chi2 = np.linspace(0, 2*np.pi, nchi+1)
                    for i in range(ntheta+1):
                        theta2 = np.arccos(v1[i, 2])
                        phi2 = np.arctan2(v1[i, 1], v1[i, 0])
                        v2[i, :, :] = np.vstack([
                            np.cos(phi2) * np.cos(theta2) * np.cos(chi2) - np.sin(phi2) * np.sin(chi2),
                            np.sin(phi2) * np.cos(theta2) * np.cos(chi2) + np.cos(phi2) * np.sin(chi2),
                            -np.sin(theta2) * np.cos(chi2)
                        ]).transpose()
                # compute
                if p == 'young':
                    self.r[ip, iplot, :] = young(self._Scart, v1, 3)
                elif p == 'comp':
                    self.r[ip, iplot, :] = comp(self._Scart, v1, 3)
                elif p == 'shear':
                    self.r[ip, iplot, :] = shear(self._Scart, [v0], [v1], 3)
                elif p == 'shear avg':
                    self.r[ip, iplot, :] = np.average(shear(self._Scart, v1, v2, 3), axis=1)
                elif p == 'shear min':
                    self.r[ip, iplot, :] = np.min(shear(self._Scart, v1, v2, 3), axis=1)
                elif p == 'shear max':
                    self.r[ip, iplot, :] = np.max(shear(self._Scart, v1, v2, 3), axis=1)
                elif p == 'poisson':
                    self.r[ip, iplot, :] = poisson(self._Scart, [v0], [v1], 3)
                elif p == 'poisson avg':
                    self.r[ip, iplot, :] = np.average(poisson(self._Scart, v1, v2, 3), axis=1)
                elif p == 'poisson min':
                    self.r[ip, iplot, :] = np.min(poisson(self._Scart, v1, v2, 3), axis=1)
                elif p == 'poisson max':
                    self.r[ip, iplot, :] = np.max(poisson(self._Scart, v1, v2, 3), axis=1)

        # Plot
        if plot_lattice == False:
            self.lattice_plot = None
        else:
            self.lattice_plot = self.lattice.matrix

        if get_data == True:
            return self
        else:
            fig = _plot2D(
                self, None, use_cartesian, u, utext, None, loop_color,
                loop_linestyle, loop_linewidth, None, figsize, None, None,
                **kwargs
            )
            # set title
            if nplane == 1 and nprop > 1:
                for ip, p in enumerate(self.property):
                    if 'shear' in p or 'young' in p: ptitle = '{} (GPa)'.format(p)
                    else: ptitle = p
                    fig.axes[ip].set_title(ptitle, y=1.05)

                if plane_definition == 'miller':
                    title = '({:<3d}{:^3d}{:>3d})'.format(plane[0,0], plane[0,1], plane[0,2])
                else:
                    title = '({:<4.1f}{:^4.1f}{:>4.1f})'.format(plane[0,0], plane[0,1], plane[0,2])

                fig.axes[0].text(
                    -0.1, 0.5, title, rotation='vertical', horizontalalignment='right',
                    verticalalignment='center', transform=fig.axes[0].transAxes,
                    fontsize=fig.axes[0].title._fontproperties._size
                )
            else:
                for inorm, norm in enumerate(plane):
                    if plane_definition == 'miller':
                        title = '({:<3d}{:^3d}{:>3d})'.format(norm[0], norm[1], norm[2])
                    else:
                        title = '({:<4.1f}{:^4.1f}{:>4.1f})'.format(norm[0], norm[1], norm[2])
                    fig.axes[inorm].set_title(title, y=1.05)

                for ip, p in enumerate(self.property):
                    if 'shear' in p or 'young' in p: ptitle = '{} (GPa)'.format(p)
                    else: ptitle = p
                    fig.axes[int(ip*nplane)].text(
                        -0.1, 0.5, ptitle, rotation='vertical', horizontalalignment='right',
                        verticalalignment='center', transform=fig.axes[int(ip*nplane)].transAxes,
                        fontsize=fig.axes[0].title._fontproperties._size
                    )
            return fig

    def plot_3D(self, property, nphi=90, ntheta=90, nchi=180,
                scale_radius=True, range_cbar=None, range_x=None, range_y=None,
                range_z=None, u=None, utext=None, use_cartesian=True,
                plot_lattice=False,  colormap='jet', title='default',
                figsize=[6.4, 4.8], plot_lib='matplotlib', **kwargs):
        """
        Plot 3D crystal elastic properties on the basis of the elastic tensor.
        The 3D surface is colored by the elastic property.

        Options:

        * "young": Young's modulus in GPa.
        * "comp": Linear compressibility.
        * "shear avg": Averaged shear modulus in the normal plane, in GPa.
        * "shear min": Minimum shear modulus in the normal plane, in GPa.
        * "shear max": Maximum shear modulus in the normal plane, in GPa.
        * "poisson avg": Averaged Poisson ratio in the normal plane.
        * "poisson min": Minimum Poisson ratio in the normal plane.
        * "poisson max": Maximum Poisson ratio in the normal plane.

        Args:
            property (str): Property to compute. See options above.
            nphi (int): Resolution of azimuth angle :math:`\\phi` on xy plane,
                in radian :math:`[0, 2\\pi)`.
            ntheta (int): Resolution of polar angle :math:`\\theta` in radian
                :math:`[0, 2\\pi)`. In practice only half of the points defined
                in :math:`[0, \\pi)` are used.
            nchi (int): Resolution of auxiliary angle  :math:`\\chi`, in radian
                :math:`[0, 2\\pi)`. Shear modulus and Poisson ratio only.
            scale_radius (bool): To scale the radius by values of the elastic
                property, or plot a sphere with radius = 1.
            range_cbar, range_x, range_y, range_z (list[float,float]): *Not
                suggested* Explicitly specifying the ranges of colorbar, x, y
                and z axes.
            u (numpy.ndarray): 3\*1 or nu\*3 array of vectors to be added into
                the figure. A line segment with doubled radius is plotted to
                indicate the vector. Or 'max' / 'min' / 'bothends', plot
                vectors corresponding to the max and min of elastic properties.
                For 'bothends', min first.
            utext (list[str] | str): A string or a list of string. Used to mark
                the vector. If ``None``, the input ``u`` and the corresponding
                value is annotated. If 'value', the value is annotated.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)
            plot_lattice (bool): Draw the lattice box around the 3D surface.
            colormap (str): Colormap name.
            title (str|None): Plot title. 'default' for default values and
                'None' for no title.
            figsize (list): Figure size in inches. For plotly, 300 ppi is
                assumed.
            plot_lib (bool): 'matplotlib' or 'plotly'. By default use
                `matplotlib <https://matplotlib.org/>`_. Alternatively use
                `plotly <https://plotly.com/>`_ to enable interactive
                inspections in Jupyter Notebook.
            \*\*kwargs: Parameters passed to ``Axes3D.view_init`` (matplotlib)
                or ``fig.update_layout()`` (plotly). Only camera position
                keywords are suggested.

        Returns:
            fig (Figure): Matplotlib or plotly ``Figure`` object. It ``axes``
                is a 2\*1 list of matplotlib ``Axis`` object. The first element
                is the axis of colorbar. The second is the axis of 3D surface.
                ``None`` if ``plot_lib = 'plotly'``.

        .. note::

            When using ``u='max'``, only one vector will be plotted if the same
            value is repeated along different directions. Same for ``min`` and
            ``bothends``.

        """
        import warnings

        # sanity check
        property_list = ['young', 'comp', 'shear avg', 'shear min', 'shear max',
                         'poisson avg', 'poisson min', 'poisson max']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        if np.all(self.lattice==None):
            if use_cartesian != True or plot_lattice == True:
                raise Exception('Lattice matrix not available.')

        if np.all(plot_lib!=None) \
        and (plot_lib.lower() != 'matplotlib' and plot_lib.lower() != 'plotly'):
            raise ValueError("Unknown plot library : '{}'. You must use 'matplotlib' or 'plotly'.".format(plot_lib))

        # Generate mesh in polar coordinate
        phi = np.linspace(0, 2*np.pi, nphi+1)
        theta = np.linspace(0, np.pi, ntheta//2+1)
        self.X = np.cos([phi]).transpose() @ np.sin([theta])
        self.Y = np.sin([phi]).transpose() @ np.sin([theta])
        self.Z = (np.zeros([nphi+1,1])+1) @ np.cos([theta])
        # Get 3D plot
        v0 = np.vstack([self.X.flatten(), self.Y.flatten(), self.Z.flatten()]).transpose()
        self.R = np.zeros([len(v0),], dtype=float)
        if 'shear' in property or 'poisson' in property:
            v1 = np.zeros([len(v0), nchi+1, 3], dtype=float)
            chi = np.linspace(0, 2*np.pi, nchi+1)
            sinp = np.sin(phi)
            cosp = np.cos(phi)
            sint = np.sin(theta)
            cost = np.cos(theta)
            sinc = np.sin(chi)
            cosc = np.cos(chi)
            for i in range(len(v0)):
                iphi = i // len(theta)
                itheta = i % len(theta)
                v1[i, :, :] = np.vstack([
                    cosp[iphi]*cost[itheta]*cosc - sinp[iphi]*sinc,
                    sinp[iphi]*cost[itheta]*cosc + cosp[iphi]*sinc,
                    -sint[itheta]*cosc
                ]).transpose()
            del sinp, cosp, sint, cost, sinc, cosc # save memory
        else:
            v1 = None; chi = None

        if property == 'young':
            self.R = young(self._Scart, v0, 3)
        elif property == 'comp':
            self.R = comp(self._Scart, v0, 3)
        elif property == 'shear avg':
            self.R = np.average(shear(self._Scart, v0, v1, 3), axis=1)
        elif property == 'shear min':
            self.R = np.min(shear(self._Scart, v0, v1, 3), axis=1)
        elif property == 'shear max':
            self.R = np.max(shear(self._Scart, v0, v1, 3), axis=1)
        elif property == 'poisson avg':
            self.R = np.average(poisson(self._Scart, v0, v1, 3), axis=1)
        elif property == 'poisson min':
            self.R = np.min(poisson(self._Scart, v0, v1, 3), axis=1)
        elif property == 'poisson max':
            self.R = np.max(poisson(self._Scart, v0, v1, 3), axis=1)

        del v0, v1, phi, theta, chi # save memory
        self.R = np.reshape(self.R, self.X.shape)
        if scale_radius == True:
            self.X = self.R * self.X
            self.Y = self.R * self.Y
            self.Z = self.R * self.Z
        self.scale_radius = scale_radius

        # Get 1D plot
        if np.all(u!=None):
            if isinstance(u, str):
                if u.lower() == 'max':
                    irow = np.argmax(self.R)//np.shape(self.R)[1]
                    icol = np.argmax(self.R)%np.shape(self.R)[1]
                    uplt = np.array([
                        self.X[irow, icol], self.Y[irow, icol], self.Z[irow, icol],
                    ], dtype=float)
                elif u.lower() == 'min':
                    irow = np.argmin(self.R)//np.shape(self.R)[1]
                    icol = np.argmin(self.R)%np.shape(self.R)[1]
                    uplt = np.array([
                        self.X[irow, icol], self.Y[irow, icol], self.Z[irow, icol],
                    ], dtype=float)
                elif u.lower() == 'bothends':
                    irow1 = np.argmin(self.R)//np.shape(self.R)[1]
                    icol1 = np.argmin(self.R)%np.shape(self.R)[1]
                    irow2 = np.argmax(self.R)//np.shape(self.R)[1]
                    icol2 = np.argmax(self.R)%np.shape(self.R)[1]
                    uplt = np.array([
                        [self.X[irow1, icol1], self.Y[irow1, icol1], self.Z[irow1, icol1]],
                        [self.X[irow2, icol2], self.Y[irow2, icol2], self.Z[irow2, icol2]],
                    ], dtype=float)
                else:
                    raise ValueError("Unknown u value: '{}'".format(u))
                use_cartesian = True
            else:
                uplt = np.array(u, dtype=float) # cartesian u, scaled

            if len(np.shape(uplt)) == 1:
                uplt = np.array([uplt])

            if use_cartesian != True:
                for i in range(len(uplt)):
                    uplt[i, :] = self.lattice.get_cartesian_coords(uplt[i, :])
                    uplt[i, :] = uplt[i, :] / np.linalg.norm(uplt[i, :])
            else:
                for i in range(len(uplt)):
                    uplt[i, :] = uplt[i, :] / np.linalg.norm(uplt[i, :])

            Ru = self.get_1D(property, uplt, nchi, use_cartesian=True)
            if len(uplt) == 1:
                Ru = [Ru]
            # Plot setups
            if scale_radius == True:
                for i in range(len(uplt)):
                    uplt[i, :] = 2 * Ru[i] * uplt[i, :]
            else:
                for i in range(len(uplt)):
                    uplt[i, :] = 2 * uplt[i, :]
            if np.all(utext!=None):
                if isinstance(utext, str):
                    if utext.lower() == 'value':
                        utext = ['{:.2f}'.format(R) for R in Ru]
                    else:
                        utext = [utext for i in range(len(uplt))]
                else:
                    if len(utext) != len(uplt):
                        warnings.warn('Length of vector annotations should be equal to length of vectors. Default annotations are used',
                                      stacklevel=2)
                        utext = None

            if np.all(utext==None):
                if len(uplt) == 1:
                    u = np.array([u], dtype=float)
                utext = []
                for v, r in zip(u, Ru):
                    utext.append('[{:<4.1f}{:^4.1f}{:>4.1f}]  {:.2f}'.format(
                        v[0], v[1], v[2], r))
        else:
            uplt = None; utext=[]

        self.u = uplt
        self.utext = utext

        # Get lattice plot
        if plot_lattice == False:
            self.lattice_plot = None
        else:
            self.lattice_plot = self.lattice.matrix

        # set title
        if np.all(title!=None):
            if title.lower() == 'default':
                if 'comp' in property or 'poisson' in property:
                    title = property
                else:
                    title = '{} (GPa)'.format(property)

        # Select plot lib
        if np.all(plot_lib!=None):
            if plot_lib.lower() == 'matplotlib':
                fig = _plot3D_mplib(self, range_cbar, range_x, range_y, range_z,
                                    colormap, figsize, fig=None, ax_index=None,
                                    Rmax=None, **kwargs)
                if np.all(title!=None):
                    fig.axes[0].set_title(title)
            else:
                fig = _plot3D_plotly(self, range_cbar, range_x, range_y, range_z,
                                     colormap, figsize, **kwargs)
                if np.all(title!=None):
                    fig.update_layout(title=title)
            return fig
        else: # Developer settings. return all the variables required by matplotlib.
            return self

    def transform(self, new_lattice):
        """
        Transform (rotate) the lattice and change matrix elements. Useful when
        lattice information is available.

        Args:
            new_lattice (numpy.ndarray) 3\*3 array of new lattice.

        Returns:
            self (Tensor3D)
        """
        from scipy.spatial.transform import Rotation as Rot
        from pymatgen.core.lattice import Lattice
        import copy

        # sanity check
        if np.all(self.lattice==None):
            raise Exception('Lattice information not available.')
        ## if it is a rotation
        nlatt = np.array(new_lattice, dtype=float)
        for i in range(3):
            unorm = np.linalg.norm(self.lattice.matrix[i, :])
            for j in range(i, 3):
                if unorm - np.linalg.norm(nlatt[j, :]) < 1e-6:
                    nlatt[[i, j], :] = nlatt[[j, i], :]
                    break
                else:
                    continue

        rmx = nlatt.transpose() @ np.linalg.inv(self.lattice.matrix.transpose())
        if np.abs(np.linalg.det(rmx))-1 > 1e-6:
            raise ValueError('Determinant of rotation matrix does not equal to 1. Not a rotation. Check your lattice input.')

        # Rotation of lattice ---> rotation of cartesian base vectors
        rot = Rot.from_matrix(rmx)
        newbase = rot.inv().apply(np.eye(3)) # New base represented by the old coordinate

        voigt = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]], dtype=int)
        newC = np.zeros([3, 3, 3, 3], dtype=float)
        for m in range(3):
            for n in range(3):
                for o in range(3):
                    for p in range(3):
                        C = 0.
                        for i in range(3):
                            for j in range(3):
                                for k in range(3):
                                    for l in range(3):
                                        C += newbase[m, i]*newbase[n, j]*newbase[o, k]*newbase[p, l]*self._C[voigt[i, j], voigt[k, l]]
                        newC[m, n, o, p] = C

        self._Ccart = copy.deepcopy(newC)
        self._C = cart2voigt(newC)
        self._S = np.linalg.inv(self._C)
        self._Scart = voigt2cart(self._S)
        self.lattice = Lattice(new_lattice)
        return self


class Tensor2D():
    """
    2D elastic tensor and related properties. Periodic boundary conditions are
    consistent with CRYSTAL definitions, i.e., xy-periodic and slab vertical to
    z. Basic units: GPa, m

    Args:
        matrix (numpy.ndarray): 3\*3 compliance or stiffness matrix.
        lattice (numpy.ndarray): 2\*2 lattice matrix.
        is_compliance (bool): Compliance or stiffness matrix used as input.
    """
    def __init__(self, matrix, lattice=None, is_compliance=True):
        import warnings
        from pymatgen.core.lattice import Lattice

        if np.all(lattice==None):
            warnings.warn("Lattice information not available, Miller indices and fractional coordinates disabled.",
                          stacklevel=2)
            self.lattice = None
        else:
            latt = np.eye(3)*500
            latt[0:2, 0:2] = lattice
            self.lattice = Lattice(np.array(latt, dtype=float))

        if is_compliance == True:
            self._S = np.array(matrix, dtype=float)
            self._C = np.linalg.inv(self._S)
        else:
            self._C = np.array(matrix, dtype=float)
            self._S = np.linalg.inv(self._C)

        # convert matrices from Voigt notation to Cartesian notation
        self._Ccart = voigt2cart(self._C)
        self._Scart = voigt2cart(self._S)

    @property
    def stiffness(self):
        """
        3\*3 stiffness matrix in Voigt annotation. Basic units: GPa, m.
        """
        return self._C

    @property
    def compliance(self):
        """
        3\*3 compliance matrix in Voigt annotation. Basic units: GPa, m.
        """
        return self._S

    @classmethod
    def from_file(cls, output, thickness):
        """
        Read elastic tensor from CRYSTAL output file and generate ``Tensor2D``
        object. Calls the ``crystal_io.Crystal_output.get_elatensor()`` method.
        Lattice information is obtained.

        Args:
            output (str): CRYSTAL output file.
            thickness (float): Effective thickness of low dimensional materials,
                in :math:`\\AA`.
        Returns:
            cls (Tensor2D)
        """
        from CRYSTALpytools.crystal_io import Crystal_output

        out = Crystal_output(output)
        if out.get_dimensionality() != 2:
            raise Exception("Dimensionality error. Input system is not a 3D system.")

        return cls(matrix=out.get_elatensor(thickness),
                   lattice=out.get_lattice(initial=False)[0:2,0:2],
                   is_compliance=False)

    def get_1D(self, property, u, use_cartesian=True):
        """
        Compute properties along the given vector. Available properties are
        defined by the ``property`` variable. To get shear modulus / Poisson
        ratio of 2 specified vectors, use ``get_pair``.

        Options:

        * "young": Young's modulus in GPa.
        * "comp": Linear compressibility.
        * "shear": Shear modulus between the vector and its normal vector in plane, in GPa.
        * "poisson": Poisson ratio between the vector and its normal vector in plane.

        Args:
            property (str): Property to compute. See options above.
            u (numpy.ndarray): 2\*1 or nu\*2 array of vectors.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)

        Returns:
            e (float | numpy.ndarray): float or nu\*1, some sort of property.
        """
        # sanity check
        property_list = ['young', 'comp', 'shear', 'poisson']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        if use_cartesian != True and np.all(self.lattice==None):
            raise Exception('Lattice matrix not available.')

        u = np.array(u, dtype=float, ndmin=2)
        u3d = np.hstack([u[:, 0:2], [[0] for i in range(len(u))]])

        if use_cartesian != True:
            for i in range(len(u)):
                u[i, :] = self.lattice.get_cartesian_coords(u3d[i])[0:2]
                u[i, :] = u[i, :] / np.linalg.norm(u[i, :])
        else:
            for i in range(len(u)):
                u[i, :] = u[i, :] / np.linalg.norm(u[i, :])

        # get the normal vector
        if 'shear' in property or 'poisson' in property:
            v = np.cross(u3d, [0, 0, 1])[:, 0:2]
            for i in range(len(v)):
                v[i, :] = v[i, :] / np.linalg.norm(v[i, :])
            # add 1 extra dimension to v (nu*1*3)
            v = np.expand_dims(v, axis=1)

        if property == 'young':
            e = young(self._Scart, u, 2)
        elif property == 'comp':
            e = comp(self._Scart, u, 2)
        elif property == 'shear':
            e = shear(self._Scart, u, v, 2)
            e = np.array([i[0] for i in e])
        elif property == 'poisson':
            e = poisson(self._Scart, u, v, 2)
            e = np.array([i[0] for i in e])

        if len(u) == 1:
            e = e[0]
        return e

    def get_pair(self, property, u, v, use_cartesian=True):
        """
        Compute properties between given pairs of in-plane vectors. Used for
        shear modulus and Poisson ratio. For properties along 1 given vector,
        use ``get_1D``.

        Options:

        * "shear": Shear modulus between in-plane ``u`` and ``v``, in GPa.
        * "poisson": Poisson ratio between in-plane ``u`` and ``v``.

        Args:
            property (str): Property to compute. See options above.
            u (numpy.ndarray): 2\*1 or nu\*2 array of vector 1.
            v (numpy.ndarray): 2\*1 or nu\*2 array of vector 2.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)

        Returns:
            e (float | numpy.ndarray): float or nu\*1, some sort of property.
        """
        # sanity check
        property_list = ['shear', 'poisson']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        if use_cartesian != True and np.all(self.lattice==None):
            raise Exception('Lattice matrix not available.')

        u = np.array(u, dtype=float, ndmin=2)
        v = np.array(v, dtype=float, ndmin=2)
        if u.shape[0] != v.shape[0]:
            raise ValueError("Inconsistent length of u and v.")
        u3d = np.hstack([u[:, 0:2], [[0] for i in range(len(u))]])
        v3d = np.hstack([v[:, 0:2], [[0] for i in range(len(v))]])

        if use_cartesian != True:
            for i in range(np.shape(u)[0]):
                u[i, :] = self.lattice.get_cartesian_coords(u3d[i])[0:2]
                v[i, :] = self.lattice.get_cartesian_coords(v3d[i])[0:2]
                u[i, :] = u[i, :] / np.linalg.norm(u[i, :])
                v[i, :] = v[i, :] / np.linalg.norm(v[i, :])
        else:
            for i in range(len(u)):
                u[i, :] = u[i, :] / np.linalg.norm(u[i, :])
                v[i, :] = v[i, :] / np.linalg.norm(v[i, :])
        # orthogonality
        for i in range(len(u)):
            if np.abs(np.dot(u[i, :], v[i, :])) > 1e-6:
                raise ValueError("Vector pair {:d} (from 0) is non-orthogonal.".format(i))

        # add 1 extra dimension to v (nu*1*3)
        v = np.expand_dims(v, axis=1)
        if property == 'shear':
            e = shear(self._Scart, u, v, 2)
        elif property == 'poisson':
            e = poisson(self._Scart, u, v, 2)

        if len(u) == 1:
            e = e[0, 0]
        else:
            e = e[:, 0]
        return e

    def plot_2D(self, property, ntheta=90, u=None, utext=None,
                use_cartesian=True, plot_lattice=True, loop_color=None,
                loop_linestyle=None, loop_linewidth=None, layout=None,
                figsize=[6.4, 4.8], get_data=False, **kwargs):
        """
        Plot 2D crystal elastic properties of the system.

        Options:

        * "young": Young's modulus in GPa.
        * "comp": Linear compressibility.
        * "shear": In-plane shear modulus between 2 vectors in the plane, in GPa.
        * "poisson": In-plane Poisson ratio between 2 vectors in the plane.

        This method can be used for multi-property plots. Subplot titles are
        added automatically. The layout of subplots is, by default, 1\*nProp.

        .. note::

            ``u`` is defined by 2D vectors.

        Args:
            property (str|list): Property(ies) to compute. See options above.
            plane (numpy.ndarray): 3\*1 or nplane\*3 array of planes.
            ntheta (int): Resolution of azimuth angle :math:`\\theta` on plane,
                in radian :math:`[0, 2\\pi)`.
            nchi (int): Resolution of auxiliary angle  :math:`\\chi`, in radian
                :math:`[0, 2\\pi)`. Shear modulus and Poisson ratio only.
            plane_definition (str): The method used to define the plane. See
                options above.
            u (numpy.ndarray): 2\*1 or nu\*2 array of vectors to be added into
                the figure. A line segment with doubled radius is plotted to
                indicate the vector. Or 'max' / 'min' / 'bothends', plot
                vectors corresponding to the max and min of elastic properties.
                For 'bothends', min first.
            utext (list[str] | str): A string or a list of string. Used to mark
                The vector. If ``None`` or 'value', the value is annotated.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)
            plot_lattice (bool): Draw lattice base vectors to indicate
                orientation of the plane.
            loop_color (str|list): Color of 2D loops. If only one string is
                given, apply it to all loops. 'None' for default values
                (matplotlib Tableau Palette).
            loop_linestyle (str|list): Linestyle of 2D loops. If only one
                string is given, apply it to all plots. 'None' for default
                values ('-').
            loop_linewidth (str|list): Linewidth of band structure. If only one
                number is given, apply it to all plots. For the 'multi' mode,
                1\*nSystem list of linewidth numbers. 'None' for default values
                (1.0).
            layout (list|tuple): 2\*1 list. The layout of subplots. The first
                element is nrows, the second is ncolumns. By default, 1\*nProp.
            figsize (list): The figure size specified as \[width, height\].
            get_data (bool): *Developer only* Return data or figure.
            \*\*kwargs: Parameters passed to ``PolarAxes.plot`` to customize
                the loops. Applied to all the loops.

        Returns:
            fig (Figure): Matplotlib ``Figure`` object.
        """
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        import warnings

        # sanity check
        property_list = ['young', 'comp', 'shear', 'poisson']
        self.property = np.array(property, ndmin=1)
        nprop = len(self.property)
        for i in range(nprop):
            if self.property[i].lower() not in property_list:
                raise ValueError("Unknown property input: '{}'".format(self.property[i]))
            self.property[i] = self.property[i].lower()

        if np.all(self.lattice==None):
            if use_cartesian != True or plot_lattice == True:
                raise Exception('Lattice matrix not available.')

        # get properties
        self.r = np.zeros([nprop, 1, ntheta+1], dtype=float) # to be commensurate with 3D cases
        chi = np.linspace(0, 2*np.pi, ntheta+1) # pole angle is actually chi.
        self.chi = chi
        v0 = np.array([0, 0, 1], dtype=float)
        theta = 0.
        phi = 0.
        for ip, p in enumerate(self.property):
            # get pole angles and 2D vectors. use v1
            v1 = np.vstack([
                np.cos(phi) * np.cos(theta) * np.cos(chi) - np.sin(phi) * np.sin(chi),
                np.sin(phi) * np.cos(theta) * np.cos(chi) + np.cos(phi) * np.sin(chi)
            ]).transpose()
            # get pole angles and 2D vectors. use v1, v2
            if p == 'shear' or p == 'poisson':
                # v2 is just an in-plane norm of v1
                v2 = np.zeros([ntheta+1, 1, 2], dtype=float) # another chi.
                for i in range(ntheta+1):
                    v2[i, 0, :] = np.cross([v1[i, 0], v1[i, 1], 0],
                                           [0, 0, 1])[0:2]
                    v2[i, 0, :] = v2[i, 0, :] / np.linalg.norm(v2[i, 0, :])
            # compute
            if p == 'young':
                self.r[ip, 0, :] = young(self._Scart, v1, 2)
            elif p == 'comp':
                self.r[ip, 0, :] = comp(self._Scart, v1, 2)
            elif p == 'shear':
                self.r[ip, 0, :] = shear(self._Scart, v1, v2, 2)[:, 0]
            elif p == 'poisson':
                self.r[ip, 0, :] = poisson(self._Scart, v1, v2, 2)[:, 0]

        # Plot
        if plot_lattice == False:
            self.lattice_plot = None
        else:
            self.lattice_plot = self.lattice.matrix

        if np.all(u!=None):
            if isinstance(u, str):
                ureal = u
            else:
                u = np.array(u, ndmin=2)
                ureal = [[i[0], i[1], 0.] for i in u]
        else:
            ureal = u

        self.norm = np.array([[0, 0, 1]], dtype=float)
        if get_data == True:
            return self
        else:
            fig = _plot2D(
                self, None, use_cartesian, ureal, utext, None, loop_color,
                loop_linestyle, loop_linewidth, layout, figsize, None, None,
                **kwargs
            )
            # set title
            for ip, p in enumerate(self.property):
                if 'shear' in p or 'young' in p: ptitle = '{} (GPa)'.format(p)
                else: ptitle = p
                fig.axes[ip].set_title(ptitle, y=1.05)
            return fig

    def transform(self, new_lattice):
        """
        Transform (rotate) the lattice and change matrix elements. Useful when
        lattice information is available.

        Args:
            new_lattice (numpy.ndarray) 2\*2 array of new lattice.

        Returns:
            self (Tensor2D)
        """
        from scipy.spatial.transform import Rotation as Rot
        from pymatgen.core.lattice import Lattice
        import copy

        # sanity check
        if np.all(self.lattice==None):
            raise Exception('Lattice information not available.')
        ## if it is a rotation
        nlatt = np.array(new_lattice, dtype=float)[0:2, 0:2]
        olatt = self.lattice.matrix[0:2, 0:2]
        for i in range(2):
            unorm = np.linalg.norm(olatt[i, :])
            for j in range(i, 2):
                if unorm - np.linalg.norm(nlatt[j, :]) < 1e-6:
                    nlatt[[i, j], :] = nlatt[[j, i], :]
                    break
                else:
                    continue

        latt3d = np.eye(3)*500
        latt3d[0:2, 0:2] = nlatt
        nlatt3d = np.eye(3)*500
        nlatt3d[0:2, 0:2] = np.array(new_lattice, dtype=float)[0:2, 0:2]
        rmx = latt3d.transpose() @ np.linalg.inv(self.lattice.matrix.transpose())
        if np.abs(np.linalg.det(rmx))-1 > 1e-6:
            raise ValueError('Determinant of rotation matrix does not equal to 1. Not a rotation. Check your lattice input.')

        # Rotation of lattice ---> rotation of cartesian base vectors
        rot = Rot.from_matrix(rmx)
        newbase = rot.inv().apply(np.eye(3)) # New base represented by the old coordinate

        voigt = np.array([[0, 2], [2, 1]], dtype=int)
        newC = np.zeros([2, 2, 2, 2], dtype=float)
        for m in range(2):
            for n in range(2):
                for o in range(2):
                    for p in range(2):
                        C = 0.
                        for i in range(2):
                            for j in range(2):
                                for k in range(2):
                                    for l in range(2):
                                        C += newbase[m, i]*newbase[n, j]*newbase[o, k]*newbase[p, l]*self._C[voigt[i, j], voigt[k, l]]
                        newC[m, n, o, p] = C

        self._Ccart = copy.deepcopy(newC)
        self._C = cart2voigt(newC)
        self._S = np.linalg.inv(self._C)
        self._Scart = voigt2cart(self._S)
        self.lattice = Lattice(nlatt3d)
        return self

# ----
# Basic elastic property functions.
# ----
def young(S, u, ndimen):
    """
    The equation of Young's modulus. Using ``Tensor*`` objects is recommended.

    .. math::

        E = \\frac{1}{u_{i}u_{j}u_{k}u_{l}S_{ijkl}}

    Einstein's notation is used. :math:`i,j,k,l=1,ndimen; j>i; l>k`.

    Args:
        S (numpy.ndarray): ndimen\*ndimen\*ndimen\*ndimen compliance matrix in
            Cartesian form. Base units: GPa, m.
        u (numpy.ndarray): nu\*ndimen array. Must be normalized Cartesian
            vectors.
        ndimen (int): Dimensionality, either 2 or 3.

    Returns:
        e (numpy.ndarray): nu\*1 Youngs's modulus. Base units: GPa, m.
    """
    rf = 1 / np.array([
        [1., 2., 2.],
        [2., 1., 2.],
        [2., 2., 1.]
    ])
    rf = rf[0:ndimen, 0:ndimen]
    e = 1 / np.einsum('...i,...j,...k,...l,ij,kl,ijkl->...', u, u, u, u, rf, rf, S)
    return e


def comp(S, u, ndimen):
    """
    The equation of linear compressibility. Using ``Tensor*`` objects is recommended.

    .. math::

        \\beta = u_{i}u_{j}S_{ijkk}

    Einstein's notation is used. :math:`i,j,k=1,ndimen; j>i`.

    Args:
        S (numpy.ndarray): ndimen\*ndimen\*ndimen\*ndimen compliance matrix in
            Cartesian form. Base units: GPa, m.
        u (numpy.ndarray): nu\*ndimen array. Must be normalized Cartesian
            vectors.
        ndimen (int): Dimensionality, either 2 or 3.

    Returns:
        beta (numpy.ndarray): nu\*1 Linear compressibility.
    """
    rf = 1 / np.array([
        [1., 2., 2.],
        [2., 1., 2.],
        [2., 2., 1.]
    ])
    rf = rf[0:ndimen, 0:ndimen]
    beta = np.einsum('...i,...j,ij,ijkk->...', u, u, rf, S)
    return beta


def shear(S, u, v, ndimen):
    """
    The equation of shear modulus. Using ``Tensor*`` objects is recommended.

    .. math::

        G = \\frac{1}{u_{i}v_{j}u_{k}v_{l}S_{ijkl}}

    Einstein's notation is used. :math:`i,j,k,l=1,ndimen`.

    Args:
        S (numpy.ndarray): ndimen\*ndimen\*ndimen\*ndimen compliance matrix in
            Cartesian form. Base units: GPa, m.
        u (numpy.ndarray): nu\*ndimen array. Must be normalized Cartesian
            vectors.
        v (numpy.ndarray): nu\*nv\*ndimen array. Must be normalized Cartesian
            vectors.
        ndimen (int): Dimensionality, either 2 or 3.

    Returns:
        g (numpy.ndarray): nu\*nv Shear modulus. Base units: GPa, m.
    """
    rf = 1 / np.array([
        [1., 2., 2.],
        [2., 1., 2.],
        [2., 2., 1.]
    ]) * 2
    rf = rf[0:ndimen, 0:ndimen]
    g = np.zeros([np.shape(v)[0], np.shape(v)[1]], dtype=float)
    for iu in range(len(u)):
        g[iu, :] = 1 / np.einsum('i,...j,k,...l,ij,kl,ijkl->...', u[iu], v[iu], u[iu], v[iu], rf, rf, S)

    return g


def poisson(S, u, v, ndimen):
    """
    The equation of Poisson ratio. Using ``Tensor*`` objects is recommended.

    .. math::

        \\nu = -\\frac{u_{i}u_{j}v_{k}v_{l}S_{ijkl}}{u_{i}u_{j}u_{k}u_{l}S_{ijkl}}

    Einstein's notation is used. :math:`i,j,k=1,ndimen; j>i`.

    Args:
        S (numpy.ndarray): ndimen\*ndimen\*ndimen\*ndimen compliance matrix in
            Cartesian form. Base units: GPa, m.
        u (numpy.ndarray): nu\*ndimen array. Must be normalized Cartesian
            vectors.
        v (numpy.ndarray): nu\*nv\*ndimen array. Must be normalized Cartesian
            vectors.
        ndimen (int): Dimensionality, either 2 or 3.

    Returns:
        nu (numpy.ndarray): nu\*nv Poisson ratio.
    """
    rf = 1 / np.array([
        [1., 2., 2.],
        [2., 1., 2.],
        [2., 2., 1.]
    ])
    rf = rf[0:ndimen, 0:ndimen]
    nu = np.zeros([np.shape(v)[0], np.shape(v)[1]], dtype=float)
    for iu in range(len(u)):
        nu[iu, :] = -np.einsum('i,j,...k,...l,ij,kl,ijkl->...', u[iu], u[iu], v[iu], v[iu], rf, rf, S)
        nu[iu, :] =  nu[iu, :] / np.einsum('i,j,k,l,ij,kl,ijkl', u[iu], u[iu], u[iu], u[iu], rf, rf, S)
    return nu

# ----
# auxiliary functions
# ----
def tensor_from_file(output, conventional_lattice=True, *thickness):
    """
    A wrapper function to get tensor objects of corresponding dimensionality.

    Args:
        output (str): CRYSTAL output file.
        conventional_lattice (bool): 3D only. Get conventional lattice.
        thickness (float): 2D only. Effective thickness of low dimensional
            materials, in :math:`\\AA`.

    Returns:
        t (Tensor3D | Tensor2D): Tensor objects.
    """
    from CRYSTALpytools.crystal_io import Crystal_output

    out = Crystal_output(output)
    if out.get_dimensionality() == 3:
        t = Tensor3D.from_file(output, conventional_lattice=conventional_lattice)
    elif out.get_dimensionality() == 2:
        if len(thickness) > 0:
            t = Tensor2D.from_file(output, thickness[0])
        else:
            t = Tensor2D.from_file(output, [])
    else:
        raise ValueError('1D or 0D systems. No available tensor objects.')
    return t


def voigt2cart(V):
    """
    Convert stiffness / compliance matrix in Voigt representation into
    Cartesian representation. For 3D: 6\*6 to 3\*3\*3\*3. For 2D: 3\*3 to
    2\*2\*2\*2.

    Args:
        V (numpy.ndarray): Voigt represented matrix.

    Returns:
        C (numpy.ndarray): Cartesian represented matrix.
    """
    dimen = int(np.shape(V)[0] / 3 + 1)
    if dimen == 3:
        voigt = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]], dtype=int)
    elif dimen == 2:
        voigt = np.array([[0, 2], [2, 1]], dtype=int)
    C = np.zeros([dimen,dimen,dimen,dimen], dtype=float)

    for i in range(dimen):
        for j in range(dimen):
            for k in range(dimen):
                for l in range(dimen):
                    C[i, j, k, l] = V[voigt[i, j], voigt[k,l]]
    return C


def cart2voigt(C):
    """
    Convert 3\*3\*3\*3  stiffness / compliance matrix in Cartesian
    representation into 6\*6 Voigt representation.

    Args:
        C (numpy.ndarray): Cartesian represented matrix.

    Returns:
        V (numpy.ndarray): Voigt represented matrix.
    """
    dimen = np.shape(C)[0]
    if dimen == 3:
        voigt = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]], dtype=int)
    elif dimen == 2:
        voigt = np.array([[0, 2], [2, 1]], dtype=int)
    V = np.zeros([int((dimen-1)*3),int((dimen-1)*3)], dtype=float)

    for i in range(dimen):
        for j in range(dimen):
            for k in range(dimen):
                for l in range(dimen):
                    V[voigt[i, j], voigt[k,l]] = C[i, j, k, l]
    return V


# ----
# Tensor plot functions. Private.
# ----
def _plot3D_mplib(tensor, range_cbar, range_x, range_y, range_z, colormap,
                  figsize, fig, ax_index, Rmax, **kwargs):
    """
    Get 3D plot using matplotlib. Private.

    Args:
        tensor3d (Tensor3D): 3D tensor object modified by the ``plot_3D``
            method. Must has 'R', 'X', 'Y', 'Z', 'scale_radius', 'u', 'utext'
            and 'lattice_plot' attributes.
        range_cbar, range_x, range_y, range_z (list[float,float]): *Not
            suggested* Explicitly specifying the ranges of colorbar, x, y and z
            axes.
        colormap (str): Colormap name.
        figsize (list): Figure size in inches.
        fig (Figure): *Developer Only*, matplotlib Figure class.
        ax_index (int): *Developer Only*, index of the axis in ``fig.axes``.
        Rmax (numpy.ndarray): *Developers Only* Auxiliary data to plot lattices
            of the same scale. The length of lattice vectors are arbitrary and
            no intersection with plot is ensured.
        \*\*kwargs: Parameters passed to ``Axes3D.view_init`` Only camera
            position keywords are suggested.

    Returns:
        fig (Figure): Matplotlib figure object
    """
    import matplotlib.pyplot as plt
    from matplotlib import colormaps, colors, cm
    from mpl_toolkits.mplot3d import Axes3D, axes3d

    if np.all(fig==None):
        fig, ax = plt.subplots(1, 1, figsize=figsize, layout='tight',
                               subplot_kw={'projection' : '3d'})
    else:
        if np.all(ax_index==None):
            raise ValueError('For subplots, the axes index must be set.')
        ax = fig.axes[int(ax_index)]

    if np.all(range_cbar==None):
        rmin = np.min(tensor.R)
        rmax = np.max(tensor.R)
    else:
        rmin = np.min(range_cbar[0])
        rmax = np.max(range_cbar[1])
    norm = colors.Normalize(vmin=rmin, vmax=rmax, clip=False)

    ax.plot_surface(
        tensor.X,
        tensor.Y,
        tensor.Z,
        rstride=1,
        cstride=1,
        facecolors=colormaps[colormap](norm(tensor.R)),
        antialiased=True,
        alpha=0.75,
    )

    m = cm.ScalarMappable(cmap=colormaps[colormap], norm=norm)
    m.set_array(np.linspace(rmin, rmax, 100))
    colorbar = fig.colorbar(m, ax=ax, location="right", shrink=0.45, pad=0.1)
    # colorbar.set_ticks(np.linspace(rmin, rmax, 5))

    # Plot 1D vectors
    u = tensor.u; utext = tensor.utext
    if np.all(u!=None):
        for v, text in zip(u, utext):
            ax.plot([0, v[0]], [0, v[1]], [0, v[2]], color='k', linewidth=1)
            ax.text(v[0], v[1], v[2], text)
    else:
        u = np.zeros([1, 3])

    # Plot lattice
    if np.all(tensor.lattice_plot!=None):
        if tensor.scale_radius == True:
            if np.all(Rmax==None):
                r = np.max(np.abs(tensor.R))
            else:
                r = Rmax
        else:
            r = 1
        lens = [np.linalg.norm(i) for i in tensor.lattice_plot]
        shrink = 2*r / np.min(lens)
        plot_latt = tensor.lattice_plot * shrink
        points = np.array([[-0.5, -0.5, -0.5], [0.5, -0.5, -0.5],
                           [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
                           [-0.5, -0.5, -0.5], [-0.5, -0.5, 0.5],
                           [0.5, -0.5, 0.5], [0.5, 0.5, 0.5],
                           [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5],
                           [0.5, -0.5, 0.5], [0.5, -0.5, -0.5],
                           [0.5, 0.5, -0.5], [0.5, 0.5, 0.5],
                           [-0.5, 0.5, 0.5], [-0.5, 0.5, -0.5]])
        points = points @ plot_latt
        ax.plot(points[:, 0], points[:, 1], points[:, 2], color='tab:gray', linewidth=1)
    else:
        points = np.zeros([1,3])

    # Make the planes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # Make the grid lines transparent
    # ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    # ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    # ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    # Fixing limits
    if tensor.scale_radius == True:
        xmx = np.max(np.abs(np.concatenate([points[:, 0], [rmax], u[:, 0]])))
        ymx = np.max(np.abs(np.concatenate([points[:, 1], [rmax], u[:, 1]])))
        zmx = np.max(np.abs(np.concatenate([points[:, 2], [rmax], u[:, 2]])))
    else:
        xmx = np.max(np.abs(np.concatenate([points[:, 0], [1], u[:, 0]])))
        ymx = np.max(np.abs(np.concatenate([points[:, 1], [1], u[:, 1]])))
        zmx = np.max(np.abs(np.concatenate([points[:, 2], [1], u[:, 2]])))
    if np.all(range_x==None):
        range_x = [-xmx, xmx]
    if np.all(range_y==None):
        range_y = [-ymx, ymx]
    if np.all(range_z==None):
        range_z = [-zmx, zmx]

    ax.set_xlim(np.min(range_x), np.max(range_x))
    ax.set_ylim(np.min(range_y), np.max(range_y))
    ax.set_zlim3d(np.min(range_z), np.max(range_z))

    ax.locator_params(nbins=5)  # tight=True,
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    # Fix aspect ratio
    ax.set_aspect('equal')
    if len(kwargs.keys()) > 0:
        ax.view_init(**kwargs)

    return fig


def _plot3D_plotly(tensor, range_cbar, range_x, range_y, range_z, colormap,
                   figsize, **kwargs):
    """
    Get 3D plot using plotly. Private.

    Args:
        tensor3d (Tensor3D): 3D tensor object modified by the ``plot_3D``
            method. Must has 'R', 'X', 'Y', 'Z', 'scale_radius', 'u', 'utext'
            and 'lattice_plot' attributes.
        range_cbar, range_x, range_y, range_z (list[float,float]): *Not
            suggested* Explicitly specifying the ranges of colorbar, x, y and z
            axes.
        colormap (str): Colormap name.
        figsize (list): Figure size in inches. For plotly, 300 ppi is assumed.
        \*\*kwargs: Parameters passed to ``fig.update_layout()`` Only camera
            position keywords are suggested.

    Returns:
        fig (Figure): Plotly Figure object.
    """
    import plotly.graph_objects as go

    if np.all(range_cbar==None):
        rmin = np.min(tensor.R)
        rmax = np.max(tensor.R)
    else:
        rmin = np.min(range_cbar[0])
        rmax = np.max(range_cbar[1])

    surface = go.Surface(x=tensor.X, y=tensor.Y, z=tensor.Z, surfacecolor=tensor.R,
                         colorscale=colormap, cmin=rmin,cmax=rmax)
    layout = go.Layout(
        title=''
    )
    fig = go.Figure(data=[surface], layout=layout)

    # Plot 1D vectors
    u = tensor.u; utext = tensor.utext
    if np.all(u!=None):
        for v, vlabel in zip(u, utext):
            fig.add_trace(
                go.Scatter3d(
                    x = [0, v[0]],
                    y = [0, v[1]],
                    z = [0, v[2]],
                    mode = 'lines+text',
                    line_color = 'black',
                    line_width = 2,
                    text = ['', vlabel],
                    textposition = "top center",
                    showlegend = False,
                )
            )
    else:
        u = np.zeros([1, 3])

    # Plot lattice
    if np.all(tensor.lattice_plot!=None):
        if tensor.scale_radius == True:
            r = np.max(np.abs(tensor.R))
        else:
            r = 1
        lens = [np.linalg.norm(i) for i in tensor.lattice_plot]
        shrink = 2*r / np.min(lens)
        plot_latt = tensor.lattice_plot * shrink
        points = np.array([[-0.5, -0.5, -0.5], [0.5, -0.5, -0.5],
                           [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
                           [-0.5, -0.5, -0.5], [-0.5, -0.5, 0.5],
                           [0.5, -0.5, 0.5], [0.5, 0.5, 0.5],
                           [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5],
                           [0.5, -0.5, 0.5], [0.5, -0.5, -0.5],
                           [0.5, 0.5, -0.5], [0.5, 0.5, 0.5],
                           [-0.5, 0.5, 0.5], [-0.5, 0.5, -0.5]])
        points = points @ plot_latt
        fig.add_trace(
            go.Scatter3d(
                x = points[:, 0],
                y = points[:, 1],
                z = points[:, 2],
                mode='lines',
                name='',
                line_color = 'grey',
                line_width = 2,
            )
        )
    else:
        points = np.zeros([1,3])

    # Fixing limits. Fix text display issue
    if tensor.scale_radius == True:
        xmx = np.max(np.abs(np.concatenate([points[:, 0], [rmax], u[:, 0]])))
        ymx = np.max(np.abs(np.concatenate([points[:, 1], [rmax], u[:, 1]])))
        zmx = np.max(np.abs(np.concatenate([points[:, 2], [rmax], u[:, 2]])))
    else:
        xmx = np.max(np.abs(np.concatenate([points[:, 0], [1], u[:, 0]])))
        ymx = np.max(np.abs(np.concatenate([points[:, 1], [1], u[:, 1]])))
        zmx = np.max(np.abs(np.concatenate([points[:, 2], [1], u[:, 2]])))
    if np.all(range_x==None):
        range_x = [-xmx, xmx]
    if np.all(range_y==None):
        range_y = [-ymx, ymx]
    if np.all(range_z==None):
        range_z = [-zmx, zmx]

    fig.update_layout(
        width=int(figsize[0]*300),
        height=int(figsize[1]*300),
        scene=dict(
            xaxis=dict(
                showbackground=True,
                backgroundcolor='white',
                gridcolor='lightgrey',
                zerolinecolor='black',
                range=[np.min(range_x), np.max(range_x)],),
            xaxis_title='X',
            yaxis=dict(
                showbackground=True,
                backgroundcolor='white',
                gridcolor='lightgrey',
                zerolinecolor='black',
                range=[np.min(range_y), np.max(range_y)],),
            yaxis_title='Y',
            zaxis=dict(
                showbackground=True,
                backgroundcolor='white',
                gridcolor='lightgrey',
                zerolinecolor='black',
                range=[np.min(range_z), np.max(range_z)],
            ),
            zaxis_title='Z',
        ),
        scene_aspectmode='manual',
        scene_aspectratio=dict(x=xmx, y=ymx, z=zmx),
    )
    if len(kwargs.keys()) > 0:
        fig.update_layout(**kwargs)
    return fig


def _plot2D(tensor, rmax, use_cartesian, u, utext, loop_label, loop_color,
            loop_linestyle, loop_linewidth, layout, figsize, fig, ax_index,
            **kwargs):
    """
    Get 2D plots of properties with matplotlib. Private.

    Args:
        tensor (Tensor3D|Tensor2D): Tensor classes with data, attribute 'r',
            shape: nProp\*nPlane\*nTheta. 'norm', norm vectors, shape nPlane\*3.
            'theta', pole angle. 'lattice', lattice.
        rmax (list|None): 1\*nProp list. Explicitly set the plot radius of
            every property.
        use_cartesian (bool): Vector is defined as cartesian or fractional
            coordinates. (*Only when lattice information is available*.)
        u (numpy.ndarray): 3\*1 or nu\*3 array of vectors to be added into the
            figure. A line segment with doubled radius is plotted to indicate
            the vector. Or 'max' / 'min' / 'bothends', plot vectors
            corresponding to the max and min of elastic properties. For
            'bothends', min first.
        utext (list[str] | str): A string or a list of string. Used to mark the
            vector. If ``None`` or 'value', the value is annotated.
        loop_label (str|list): Label of 2D loops. Used only for multi-system
            plots.
        loop_color (str|list): Color of 2D loops. If only one string is given,
            apply it to all loops. 'None' for default values (matplotlib
            Tableau Palette).
        loop_linestyle (str|list): Linestyle of 2D loops. If only one string is
            given, apply it to all plots. 'None' for default values ('-').
        loop_linewidth (str|list): Linewidth of band structure. If only one
            number is given, apply it to all plots. For the 'multi' mode,
            1\*nSystem list of linewidth numbers. 'None' for default values (1.0).
        layout (list|tuple): 2\*1 list. The layout of subplots. The first element
            is nrows, the second is ncolumns. By default, 1\*nProp or
            nProp\*nPlane.
        figsize (list): The figure size specified as \[width, height\].
        fig (Figure): *Developer Only*, matplotlib Figure class.
        ax_index (list[int]): *Developer Only*, 1D list of indices of axes in
            ``fig.axes``.
        \*\*kwargs: Parameters passed to ``PolarAxes.plot`` to customize the
            loops. Applied to all the loops.

    Returns:
        fig (Figure): Matplotlib figure object
    """
    import warnings
    import numpy as np
    import matplotlib.pyplot as plt
    from CRYSTALpytools.base.plotbase import _plot_label_preprocess

    nprop = tensor.r.shape[0]
    nplane = tensor.r.shape[1]
    ntheta = tensor.r.shape[2]

    # plot radius
    if np.all(rmax!=None):
        rmax = np.array(rmax, ndmin=1, dtype=float)
        if len(rmax) != nprop:
            warnings.warn('The specified plot radius must have the same length as plot properties.',
                          stacklevel=2)
            rmax = None
    if np.all(rmax==None):
        rmax = []
        for value in tensor.r: rmax.append(np.max(value))
        rmax = np.array(rmax, dtype=float)

    # plot layout
    if np.all(layout!=None):
        if layout[0] * layout[1] < nprop * nplane:
            warnings.warn('Inconsistent plane number and figure layout. Default layout is used.',
                          stacklevel=2)
            layout = None
        nstyle = layout[1] # nstyle is ncols
    if np.all(layout==None):
        if nplane == 1:
            nstyle = nprop
            layout = [1, nstyle]
        else:
            nstyle = nplane
            layout = [nprop, nstyle]

    if np.all(fig==None):
        fig, ax = plt.subplots(nrows=layout[0], ncols=layout[1], figsize=figsize,
                               subplot_kw={'projection' : 'polar'}, layout='tight')
        ax_index = [i for i in range(len(fig.axes))]
    else:
        if np.all(ax_index==None):
            raise ValueError("Indices of axes must be set for subplots.")
        ax_index = np.array(ax_index, ndmin=1)

    # plot styles
    doss = np.zeros([nstyle, 1, 1]) # a dummy doss input
    commands = _plot_label_preprocess(doss, loop_label, loop_color, loop_linestyle, loop_linewidth)
    keywords = ['label', 'color', 'linestyle', 'linewidth']

    # plot and deal with 1D vectors
    for iax in ax_index:
        ax = fig.axes[iax]
        irow = iax // nstyle; icol = iax % nstyle
        iprop = iax // nplane; iplane = iax % nplane
        rplt = rmax[iprop]
        theta = np.arccos(tensor.norm[iplane,2]) # polar angle of surface norm
        phi = np.arctan2(tensor.norm[iplane,1], tensor.norm[iplane,0]) # azimuth angle of surface norm
        # plot lattice vectors
        if np.all(tensor.lattice_plot!=None):
            a = tensor.lattice_plot[0, :] / np.linalg.norm(tensor.lattice_plot[0, :])
            b = tensor.lattice_plot[1, :] / np.linalg.norm(tensor.lattice_plot[1, :])
            c = tensor.lattice_plot[2, :] / np.linalg.norm(tensor.lattice_plot[2, :])
            rmarker = rplt / 20
            for v, vlabel in zip([a,b,c], ['a','b','c']):
                if np.abs(np.abs(np.dot(v, tensor.norm[iplane])) - 1) < 1e-6: # parallel
                    if np.linalg.norm(v+tensor.norm[iplane]) > 1e-6: # same direction, dot
                        ax.scatter(0, 0, c='tab:gray', marker='o')
                    else: # opposite direction, cross
                        ax.scatter(0, 0, c='tab:gray', marker='x')
                    ax.text(0, rmarker, vlabel, color='tab:gray',
                            fontstyle='italic', ha='center', va='center')
                else:
                    # get the projection of base vector on plane
                    proj = np.cross(tensor.norm[iplane], np.cross(v, tensor.norm[iplane]))
                    proj = proj / np.linalg.norm(proj)
                    cossin_chi = np.matmul(
                        np.linalg.inv(np.array([[np.cos(phi)*np.cos(theta), -np.sin(phi)],
                                                [np.sin(phi)*np.cos(theta), np.cos(phi)]])),
                        [proj[0], proj[1]]
                    )
                    chi = np.arctan2(cossin_chi[1], cossin_chi[0])
                    ax.vlines(chi, 0, rplt, colors='tab:gray', linewidth=1.5)
                    ax.text(chi, rplt+rmarker, vlabel, color='tab:gray',
                            fontstyle='italic', ha='center', va='center')

        # plot loop
        for icmd in range(4):
            if np.all(commands[icmd]==None): continue
            kwargs[keywords[icmd]] = commands[icmd][icol][0]
        ax.plot(tensor.chi, tensor.r[iprop, iplane], **kwargs)

        # deal with 1D vectors and plot
        if np.all(u!=None):
            # get cartesian u vectors
            if not isinstance(u, str):
                utmp = np.array(u, dtype=float, ndmin=2) # cartesian u, normalised
                if use_cartesian != True:
                    for i in range(len(utmp)):
                        utmp[i, :] = tensor.lattice.get_cartesian_coords(utmp[i, :])
                # normalize
                for i in range(len(utmp)):
                    utmp[i, :] = utmp[i, :] / np.linalg.norm(utmp[i, :])

                ## text
                if np.all(utext!=None):
                    utext = np.array(utext, ndmin=1)
                    if len(utext)!=len(utmp) and len(utext)!=1:
                        warnings.warn('Length of vector annotations should be equal to length of vectors. Default annotations are used',
                                      stacklevel=2)
                        utext = ['value' for i in range(len(utmp))]
                    elif len(utext)==1 and isinstance(utext[0], str): # text or 'values'
                        if utext[0].lower() == 'value':
                            utext = ['value' for i in range(len(utmp))]
                        else:
                            utext = [utext[0] for i in range(len(utmp))]
                else:
                    utext = ['value' for i in range(len(utmp))]

                # orthogonality
                uplt = []; utextplt = []; ru = []
                for iv, v in enumerate(utmp):
                    if np.abs(np.dot(v, tensor.norm[iplane])) > 1e-6: # non-orthogonal
                        continue
                    uplt.append(v)
                    # 2D system
                    if tensor.compliance.shape[0] == 3:
                        vreal = v[:2]
                    else:
                        vreal = v
                    ## values
                    if (tensor.property[iprop] == 'shear' or tensor.property[iprop] == 'poisson') \
                    and tensor.compliance.shape[0] == 6: # 2D system only compute in-plane props
                        ru.append(tensor.get_pair(
                            tensor.property[iprop], vreal, nreal, use_cartesian=True
                        ))
                    else:
                        ru.append(tensor.get_1D(
                            tensor.property[iprop], vreal, use_cartesian=True
                        ))
                    if utext[iv] == 'value':
                        utextplt.append('{:.2f}'.format(ru[-1]))
                    else:
                        utextplt.append(utext[iv])

                uplt = np.array(uplt, dtype=float)

                # convert u into polar angles
                uchi = []
                for v in uplt:
                    cossin_chi = np.matmul(
                        np.linalg.inv(np.array([[np.cos(phi)*np.cos(theta), -np.sin(phi)],
                                                [np.sin(phi)*np.cos(theta), np.cos(phi)]])),
                        [v[0], v[1]]
                    )
                    uchi.append(np.arctan2(cossin_chi[1], cossin_chi[0]))
            # text u
            else:
                ithetamx = np.argmax(tensor.r[iprop, iplane])
                ithetami = np.argmin(tensor.r[iprop, iplane])
                if u.lower() == 'max':
                    uchi = np.array([tensor.chi[ithetamx]])
                    # u text
                    if np.all(utext!=None):
                        utextplt = np.array(utext, ndmin=1)
                        if len(utextplt) != 1:
                            warnings.warn('Length of vector annotations should be equal to length of vectors. Default annotations are used',
                                           stacklevel=2)
                            utextplt[0] = 'value'
                        if utextplt[0].lower() == 'value':
                            utextplt[0] = '{:.2f}'.format(tensor.r[iprop, iplane, ithetamx])

                elif u.lower() == 'min':
                    uchi = np.array([tensor.chi[ithetami]])
                    # u text
                    if np.all(utext!=None):
                        utextplt = np.array(utext, ndmin=1)
                        if len(utextplt) != 1:
                            warnings.warn('Length of vector annotations should be equal to length of vectors. Default annotations are used',
                                           stacklevel=2)
                            utextplt[0] = 'value'
                        if utextplt[0].lower() == 'value':
                            utextplt[0] = '{:.2f}'.format(tensor.r[iprop, iplane, ithetami])

                elif u.lower() == 'bothends':
                    uchi = np.array([tensor.chi[ithetami], tensor.chi[ithetammx]])
                    # u text
                    if np.all(utext!=None):
                        utexpltt = np.array(utext, ndmin=1)
                        if len(utextplt) != 2:
                            warnings.warn('Length of vector annotations should be equal to length of vectors. Default annotations are used',
                                           stacklevel=2)
                            utextplt[0] = 'value'
                            utextplt[1] = 'value'
                        if utextplt[0].lower() == 'value':
                            utextplt[0] = '{:.2f}'.format(tensor.r[iprop, iplane, ithetami])
                        if utextplt[1].lower() == 'value':
                            utextplt[1] = '{:.2f}'.format(tensor.r[iprop, iplane, ithetamx])

                else:
                    raise ValueError("Unknown u value: '{}'.".format(u))

            # plot
            for v, vlabel in zip(uchi, utextplt):
                ax.vlines(v, 0, rplt, colors='k', linewidth=1)
                ax.text(v, rplt*1.1, vlabel, color='k', ha='center', va='center')

        # axis and title
        ax.set_rmax(rplt)
        ax.set_theta_zero_location('E')
        ax.set_rlabel_position(0.0)
        ax.grid(True)
        ax.set_xticklabels([])
        ax.set_yticks(np.round(np.linspace(0, rplt, 4), 2))
        ax.tick_params(axis='y', labelrotation=-45)
    return fig
