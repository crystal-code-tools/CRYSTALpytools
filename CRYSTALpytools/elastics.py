#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A post-processing module for elastic properties.
"""
from CRYSTALpytools import units
import numpy as np

__all__ = ['Tensor3D']


class Tensor3D():
    """
    3D elastic tensor and related properties.

    Args:
        matrix (numpy.ndarray): 6\*6 compliance or stiffness matrix. Unit: GPa:math:`^{-1}` or GPa.
        lattice (numpy.ndarray): lattice matrix.
        is_compliance (bool): Compliance or stiffness matrix used as input.
    """
    def __init__(self, matrix, lattice=None, is_compliance=True):
        import warnings
        from pymatgen.core.lattice import Lattice

        if lattice is None:
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
        voigt = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]], dtype=int)
        self._Ccart = np.zeros([3, 3, 3, 3], dtype=float)
        self._Scart = np.zeros([3, 3, 3, 3], dtype=float)
        for i in range(3):
            for j in range(3):
                v = voigt[i, j]
                for k in range(3):
                    for l in range(3):
                        self._Ccart[i, j, k, l] = self._C[v, voigt[k,l]]
                        self._Scart[i, j, k, l] = self._S[v, voigt[k,l]]

    @property
    def stiffness(self):
        """
        6\*6 stiffness matrix. Unit: GPa.
        """
        return self._C

    @property
    def compliance(self):
        """
        6\*6 compliance matrix. Unit: GPa:math:`^{-1}`.
        """
        return self._S

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
            nchi (int): Resolution of auxiliary angle  :math:`\\chi`, in radian
                :math:`[0, 2\\pi)`. Shear modulus and Poisson ratio only.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)

        Returns:
            e (float | numpy.ndarray): Some sort of property.
        """
        property_list = ['young', 'comp', 'shear', 'shear avg', 'shear min',
                         'shear max', 'poisson', 'poisson avg', 'poisson min',
                         'poisson max']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        if use_cartesian != True and self.lattice is None:
            raise Exception('Lattice matrix not available.')

        u = np.array(u, dtype=float)
        if len(np.shape(u)) == 1:
            u = np.array([u])

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
                for j in range(nchi+1):
                    v[i, j, :] = np.array([
                        np.cos(phi[i]) * np.cos(theta[i]) * np.cos(chi[j]) - np.sin(phi[i]) * np.sin(chi[j]),
                        np.sin(phi[i]) * np.cos(theta[i]) * np.cos(chi[j]) + np.cos(phi[i]) * np.sin(chi[j]),
                        -np.sin(theta[i]) * np.cos(chi[j])
                    ])

        if property == 'young':
            e = young(self._Scart, u)
        elif property == 'comp':
            e = comp(self._Scart, u)
        elif property == 'shear':
            e = shear(self._Scart, u, v)
        elif property == 'shear avg':
            e = np.average(shear(self._Scart, u, v), axis=1)
        elif property == 'shear min':
            e = np.min(shear(self._Scart, u, v), axis=1)
        elif property == 'shear max':
            e = np.max(shear(self._Scart, u, v), axis=1)
        elif property == 'poisson':
            e = poisson(self._Scart, u, v)
        elif property == 'poisson avg':
            e = np.average(poisson(self._Scart, u, v), axis=1)
        elif property == 'poisson min':
            e = np.min(poisson(self._Scart, u, v), axis=1)
        elif property == 'poisson max':
            e = np.max(poisson(self._Scart, u, v), axis=1)

        if len(u) == 1:
            e = e[0]
        return e

    def get_pair(self, property, u, v, use_cartesian=True):
        """
        Compute properties between given pairs of vectors. Used for shear
        modulus and Poisson ratio. For properties along one vector, use ``get_1D``.

        Options:

        * "shear": Shear modulus on the 'loop' normal to vector, in GPa.
        * "poisson": Poisson ratio on the 'loop' normal to vector.

        Args:
            property (str): Property to compute. See options above.
            u (numpy.ndarray): 3\*1 or nu\*3 array of vector 1.
            v (numpy.ndarray): 3\*1 or nu\*3 array of vector 2.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)

        Returns:
            e (float | numpy.ndarray): Some sort of property.
        """
        property_list = ['shear', 'poisson']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        if use_cartesian != True and self.lattice is None:
            raise Exception('Lattice matrix not available.')

        u = np.array(u, dtype=float)
        v = np.array(v, dtype=float)
        if len(np.shape(u)) == 1:
            u = np.array([u])
            v = np.array([v])
        if np.shape(u)[0] != np.shape(v)[0]:
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
            if np.abs(np.dot(u[i, :], v[i, :])) > 1e-4:
                raise ValueError("Vector pair {:d} (from 0) is non-orthogonal.".format(i))

        # add 1 extra dimension to v (nu*1*3)
        v = np.expand_dims(v, axis=1)
        if property == 'shear':
            e = shear(self._Scart, u, v)
        elif property == 'poisson':
            e = poisson(self._Scart, u, v)

        if len(u) == 1:
            e = e[0, 0]
        else:
            e = e[:, 0]
        return e

    def plot_2D(self, property, plane, ntheta=90, nchi=180,
                plane_definition='miller', u=None, utext=None,
                use_cartesian=True, layout=None, add_title=True,
                uniform_scale=True,plot_lattice=True):
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

        Args:
            property (str): Property to compute. See options above.
            plane (numpy.ndarray): 3\*1 or nplane\*3 array of planes.
            ntheta (int): Resolution of polar angle :math:`\\theta` on plane,
                in radian :math:`[0, 2\\pi)`.
            nchi (int): Resolution of auxiliary angle  :math:`\\chi`, in radian
                :math:`[0, 2\\pi)`. Shear modulus and Poisson ratio only.
            plane_definition (str): The method used to define the plane. See
                options above.
            u (numpy.ndarray): 3\*1 or nu\*3 array of vectors to be added into
                the figure. When the given u is in the plane, a line segment to
                the limit of axis is plotted to indicate the vector.
            utext (list[str] | str): A string or a list of string. Used to mark
                The vector. If ``None``, the input ``u`` and the corresponding
                value is annotated. If 'value', the value is annotated.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)
            layout (list|tuple): 2\*1 list. The layout of subplots. The first
                element is nrows, the second is ncolumns. By default, at most 3
                subplots per row is used.
            add_title (bool): Add title to subplot. Only entries of ``plane``
                can be used.
            uniform_scale (bool): Use the same radial scale (elastic properties)
                for all the plots.
            plot_lattice (bool): Draw lattice base vectors to indicate
                orientation of the plane.

        Returns:
            fig (Figure): Matplotlib ``Figure`` object.
            ax (Axis): Matplotlib ``Axis`` object or a list of it if multiple
                planes are plotted.
        """
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        import warnings

        # sanity check
        property_list = ['young', 'comp', 'shear', 'shear avg', 'shear min',
                         'shear max', 'poisson', 'poisson avg', 'poisson min',
                         'poisson max']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        planedef_list = ['miller', 'cartesian', 'fractional']
        plane_definition = plane_definition.lower()
        if plane_definition not in planedef_list:
            raise ValueError("Unknown plane definition: '{}'".format(plane_definition))

        if self.lattice is None:
            if use_cartesian != True or plot_lattice == True:
                raise Exception('Lattice matrix not available.')

        # plane definition: use normalised cartesian for plot
        planeplt = np.array(plane, dtype=float)
        if len(np.shape(plane)) == 1:
            plane = [plane]
            planeplt = np.array([planeplt])

        if plane_definition == 'miller':
            for i in range(len(planeplt)):
                planeplt[i, :] = self.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(planeplt[i, :])
                planeplt[i, :] = planeplt[i, :] / np.linalg.norm(planeplt[i, :])
        elif plane_definition == 'fractional':
            for i in range(len(planeplt)):
                planeplt[i, :] = self.lattice.get_cartesian_coords(planeplt[i, :])
                planeplt[i, :] = planeplt[i, :] / np.linalg.norm(planeplt[i, :])
        else:
            for i in range(len(planeplt)):
                planeplt[i, :] = planeplt[i, :] / np.linalg.norm(planeplt[i, :])

        # Compute each 2D plane separately
        r = np.zeros([len(planeplt), ntheta+1], dtype=float)
        for iplot, v0 in enumerate(planeplt):
            # get pole angles and 3D vectors. use (v0) v1
            theta = np.arccos(v0[2])
            phi = np.arctan2(v0[1], v0[0])
            chi = np.linspace(0, 2*np.pi, ntheta+1) # pole angle is actually chi.
            v1 = np.zeros([1, ntheta+1, 3], dtype=float)
            v1[0, :, :] = np.vstack([
                np.cos(phi) * np.cos(theta) * np.cos(chi) - np.sin(phi) * np.sin(chi),
                np.sin(phi) * np.cos(theta) * np.cos(chi) + np.cos(phi) * np.sin(chi),
                -np.sin(theta) * np.cos(chi)
            ]).transpose()

            # get pole angles and 3D vectors. use v1, v2
            if 'avg' in property or 'min' in property or 'max' in property:
                v2 = np.zeros([ntheta+1, nchi+1, 3], dtype=float) # another chi.
                chi2 = np.linspace(0, 2*np.pi, nchi+1)
                for i in range(ntheta+1):
                    theta2 = np.arccos(v[0, i, 2])
                    phi2 = np.arctan2(v[0, i, 1], v[0, i, 0])
                    v2[i, :, :] = np.vstack([
                        np.cos(phi2) * np.cos(theta2) * np.cos(chi2) - np.sin(phi2) * np.sin(chi2),
                        np.sin(phi2) * np.cos(theta2) * np.cos(chi2) + np.cos(phi2) * np.sin(chi2),
                        -np.sin(theta2) * np.cos(chi2)
                    ]).transpose()

            # compute
            if property == 'young':
                r[iplot, :] = young(self._Scart, v1)
            elif property == 'comp':
                r[iplot, :] = comp(self._Scart, v1)
            elif property == 'shear':
                r[iplot, :] = shear(self._Scart, [v0], v1)
            elif property == 'shear avg':
                r[iplot, :] = np.average(shear(self._Scart, v1, v2), axis=1)
            elif property == 'shear min':
                r[iplot, :] = np.min(shear(self._Scart, v1, v2), axis=1)
            elif property == 'shear max':
                r[iplot, :] = np.max(shear(self._Scart, v1, v2), axis=1)
            elif property == 'poisson':
                r[iplot, :] = poisson(self._Scart, [v0], v1)
            elif property == 'poisson avg':
                r[iplot, :] = np.average(poisson(self._Scart, v1, v2), axis=1)
            elif property == 'poisson min':
                r[iplot, :] = np.min(poisson(self._Scart, v1, v2), axis=1)
            elif property == 'poisson max':
                r[iplot, :] = np.max(poisson(self._Scart, v1, v2), axis=1)

        # Plot
        # Set plot framework
        if layout is not None:
            if layout[0] + layout[1] != len(planeplt):
                warnings.warn('Inconsistent plane number and figure layout. Default layout is used.',
                              stacklevel=2)
                layout = None
        if layout is None:
            layout = [len(planeplt)//3+1, 0]
            if len(planeplt) < 3:
                layout[1] = len(planeplt)
            else:
                layout[1] = 3

        fig, ax = plt.subplots(nrows=layout[0], ncols=layout[1],
                               subplot_kw={'projection' : 'polar'})
        if len(planeplt) == 1:
            ax = [ax]

        # Prepare 1D plot
        if u is not None:
            utmp = np.array(u, dtype=float) # cartesian u, normalised
            if len(np.shape(utmp)) == 1:
                utmp = np.array([utmp])

            if use_cartesian != True:
                for i in range(len(utmp)):
                    utmp[i, :] = self.lattice.get_cartesian_coords(utmp[i, :])
                    utmp[i, :] = utmp[i, :] / np.linalg.norm(utmp[i, :])
            else:
                for i in range(len(utmp)):
                    utmp[i, :] = utmp[i, :] / np.linalg.norm(utmp[i, :])

            if utext is not None:
                if isinstance(utext, str):
                    if utext.lower() == 'value':
                        utext = ['value' for i in range(len(utmp))]
                    else:
                        utext = [utext]
                if len(utext) != len(utmp):
                     warnings.warn('Length of vector annotations should be equal to length of vectors. Default annotations are used',
                                   stacklevel=2)

        # Prepare lattice
        if plot_lattice == False:
            platt = None
        else:
            platt = self.lattice.matrix

        # Set max radius
        if uniform_scale == True:
            rmax = np.max(r)
        else:
            rmax = None

        pcolor = list(mcolors.TABLEAU_COLORS.keys())
        for iplot, v0 in enumerate(planeplt):
            # plot 1D vectors
            if u is not None:
                uplt = []
                utextplt = []
                ru = []
                for ivp, vp in enumerate(utmp):
                    if np.abs(np.dot(vp, v0)) > 1e-4: # non-orthogonal
                        continue
                    uplt.append(vp)
                    if property == 'shear' or property == 'poisson':
                        ru.append(self.get_pair(property, v0, vp, use_cartesian=True))
                    else:
                        ru.append(self.get_1D(property, vp, nchi, use_cartesian=True))

                    if utext is not None:
                        if utext[ivp].lower() == 'value':
                            utextplt.append('{:.2f}'.format(ru[-1]))
                        else:
                            utextplt.append(utext[ivp])
                    else:
                        if len(utmp) == 1:
                            u = np.array([u], dtype=float)
                        else:
                            u = np.array(u, dtype=float)
                        utextplt.append('[{:<4.1f}{:^4.1f}{:>4.1f}]  {:.2f}'.format(
                            u[ivp, 0], u[ivp, 1], u[ivp, 2], ru[-1]
                        ))
                uplt = np.array(uplt, dtype=float)
            else:
                uplt=None; utextplt=[]

            # plot figures
            if add_title == True:
                if plane_definition == 'miller':
                    title = '({:<2d} {:^2d} {:>2d})'.format(
                        plane[iplot][0], plane[iplot][1], plane[iplot][2])
                elif plane_definition == 'cartesian':
                    title = 'Cartesian Norm: ({:<4.1f} {:^4.1f} {:>4.1f})'.format(
                        plane[iplot][0], plane[iplot][1], plane[iplot][2])
                else:
                    title = 'Fractional Norm: ({:<4.1f} {:^4.1f} {:>4.1f})'.format(
                        plane[iplot][0], plane[iplot][1], plane[iplot][2])
            else:
                title = ''

            if layout[0] > 1:
                irow = iplot // layout[1]
                icol = iplot % layout[1]
                ax[irow, icol] = _plot2D_single(
                    ax[irow, icol], chi, r[iplot, :], v0, pcolor[iplot%10],
                     title, rmax, uplt, utextplt, platt
                )
            else:
                ax[iplot] = _plot2D_single(
                    ax[iplot], chi, r[iplot, :], v0, pcolor[iplot%10], title,
                    rmax, uplt, utextplt, platt
                )


        if len(planeplt) == 1:
            ax = ax[0]

        return fig, ax

    def plot_3D(self, property, nphi=90, ntheta=90, nchi=180, scale_radius=True,
                u=None, utext=None, use_cartesian=True, plot_lattice=False,
                use_plotly=False):
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
            u (numpy.ndarray): 3\*1 or nu\*3 array of vectors to be added into
                the figure. A line segment with doubled radius is plotted to
                indicate the vector.
            utext (list[str] | str): A string or a list of string. Used to mark
                The vector. If ``None``, the input ``u`` and the corresponding
                value is annotated. If 'value', the value is annotated.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)
            plot_lattice (bool): Draw the lattice box around the 3D surface.
            use_plotly (bool): By default use `matplotlib <https://matplotlib.org/>`_.
                Alternatively use `plotly <https://plotly.com/>`_ to enable
                interactive inspections in Jupyter Notebook.

        Returns:
            fig (Figure): Matplotlib or plotly ``Figure`` object.
            ax (list[Axis] | None): 2\*1 list of matplotlib ``Axis`` object. The
                first element is the axis of colorbar. The second is the axis of
                3D surface. ``None`` if ``use_plotly = True``.
        """
        import warnings

        # sanity check
        property_list = ['young', 'comp', 'shear avg', 'shear min', 'shear max',
                         'poisson avg', 'poisson min', 'poisson max']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        if self.lattice is None:
            if use_cartesian != True or plot_lattice == True:
                raise Exception('Lattice matrix not available.')

        # Generate mesh in polar coordinate
        phi = np.linspace(0, 2*np.pi, nphi+1)
        theta = np.linspace(0, np.pi, ntheta//2+1)
        X = np.cos([phi]).transpose() @ np.sin([theta])
        Y = np.sin([phi]).transpose() @ np.sin([theta])
        Z = (np.zeros([nphi+1,1])+1) @ np.cos([theta])
        # Get 3D plot
        v0 = np.vstack([X.flatten(), Y.flatten(), Z.flatten()]).transpose()
        R = np.zeros([len(v0),], dtype=float)
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
                for j in range(nchi+1):
                    v1[i, j, :] = np.array([
                        cosp[iphi]*cost[itheta]*cosc[j] - sinp[iphi]*sinc[j],
                        sinp[iphi]*cost[itheta]*cosc[j] + cosp[iphi]*sinc[j],
                        -sint[itheta]*cosc[j]
                    ])
            del sinp, cosp, sint, cost, sinc, cosc # save memory
        else:
            v1 = None; chi = None

        if property == 'young':
            R = young(self._Scart, v0)
        elif property == 'comp':
            R = comp(self._Scart, v0)
        elif property == 'shear avg':
            R = np.average(shear(self._Scart, v0, v1), axis=1)
        elif property == 'shear min':
            R = np.min(shear(self._Scart, v0, v1), axis=1)
        elif property == 'shear max':
            R = np.max(shear(self._Scart, v0, v1), axis=1)
        elif property == 'poisson avg':
            R = np.average(poisson(self._Scart, v0, v1), axis=1)
        elif property == 'poisson min':
            R = np.min(poisson(self._Scart, v0, v1), axis=1)
        elif property == 'poisson max':
            R = np.max(poisson(self._Scart, v0, v1), axis=1)

        del v0, v1, phi, theta, chi # save memory
        R = np.reshape(R, X.shape)
        if scale_radius == True:
            X = R * X
            Y = R * Y
            Z = R * Z

        # Get 1D plot
        if u is not None:
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
            if utext is not None:
                if isinstance(utext, str):
                    if utext.lower() == 'value':
                        utext = ['value' for i in range(len(uplt))]
                    else:
                        utext = [utext]
                if len(utext) != len(uplt):
                    warnings.warn('Length of vector annotations should be equal to length of vectors. Default annotations are used',
                                  stacklevel=2)
                    utext = None
                else:
                    for itxt, txt in enumerate(utext):
                        if txt.lower() == 'value':
                            utext[itxt] = '{:.2f}'.format(Ru[itxt])
            if utext is None:
                if len(uplt) == 1:
                    u = np.array([u], dtype=float)
                utext = []
                for v, r in zip(u, Ru):
                    utext.append('[{:<4.1f}{:^4.1f}{:>4.1f}]  {:.2f}'.format(
                        v[0], v[1], v[2], r))
        else:
            uplt = None; utext=[]

        # Get lattice plot
        if plot_lattice == False:
            platt = None
        else:
            platt = self.lattice.matrix

        # Select plot lib
        if use_plotly == False:
            fig, ax = _plot3D_mplib(R, X, Y, Z, scale_radius, uplt, utext, platt)
        else:
            fig = _plot3D_plotly(R, X, Y, Z, scale_radius, uplt, utext, platt)
            ax = None

        return fig, ax


# ----
# Basic elastic property functions. Currently private.
# ----
def young(S, u):
    """
    The equation of Young's modulus. Private, use ``Tensor`` object.

    Args:
        S (numpy.ndarray): Compliance matrix in Cartesian form. Unit: GPa:math:`^{-1}`.
        u (numpy.ndarray): nu\*3 array. Must be normalized Cartesian vectors.

    Returns:
        e (numpy.ndarray): nu\*1 Youngs's modulus. Unit: GPa
    """
    rf = 1 / np.array([
        [1., 2., 2.],
        [2., 1., 2.],
        [2., 2., 1.]
    ])
    e = 1 / np.einsum('...i,...j,...k,...l,ij,kl,ijkl->...', u, u, u, u, rf, rf, S)
    return e


def comp(S, u):
    """
    The equation of linear compressibility. Private, use ``Tensor`` object.

    Args:
        S (numpy.ndarray): Compliance matrix in Cartesian form. Unit: GPa:math:`^{-1}`
        u (numpy.ndarray): nu\*3 array. Must be normalized Cartesian vectors.

    Returns:
        beta (numpy.ndarray): nu\*1 Linear compressibility.
    """
    rf = 1 / np.array([
        [1., 2., 2.],
        [2., 1., 2.],
        [2., 2., 1.]
    ])
    beta = np.einsum('...i,...j,ij,ijkk->...', u, u, rf, S)
    return beta


def shear(S, u, v):
    """
    The equation of shear modulus. Private, use ``Tensor`` object.

    Args:
        S (numpy.ndarray): Compliance matrix in Cartesian form.. Unit: GPa:math:`^{-1}`.
        u (numpy.ndarray): nu\*3 array. Must be normalized Cartesian vectors.
        v (numpy.ndarray): nu\*nv\*3 array. Must be normalized Cartesian vectors.

    Returns:
        g (numpy.ndarray): nu\*nv Shear modulus. Unit: GPa.
    """
    rf = 1 / np.array([
        [1., 2., 2.],
        [2., 1., 2.],
        [2., 2., 1.]
    ]) * 2
    g = np.zeros([np.shape(v)[0], np.shape(v)[1]], dtype=float)
    for iu in range(len(u)):
        g[iu, :] = 1 / np.einsum('i,...j,k,...l,ij,kl,ijkl->...', u[iu], v[iu], u[iu], v[iu], rf, rf, S)

    return g


def poisson(S, u, v):
    """
    The equation of Poisson ratio. Private, use ``Tensor`` object.

    Args:
        S (numpy.ndarray): Compliance matrix in Cartesian form. Unit: GPa:math:`^{-1}`.
        u (numpy.ndarray): nu\*3 array. Must be normalized Cartesian vectors.
        v (numpy.ndarray): nu\*nv\*3 array. Must be normalized Cartesian vectors.

    Returns:
        nu (numpy.ndarray): nu\*nv Poisson ratio.
    """
    rf = 1 / np.array([
        [1., 2., 2.],
        [2., 1., 2.],
        [2., 2., 1.]
    ])
    nu = np.zeros([np.shape(v)[0], np.shape(v)[1]], dtype=float)
    for iu in range(len(u)):
        nu[iu, :] = -np.einsum('i,j,...k,...l,ij,kl,ijkl->...', u[iu], u[iu], v[iu], v[iu], rf, rf, S)
        nu[iu, :] =  nu[iu, :] / np.einsum('i,j,k,l,ij,kl,ijkl', u[iu], u[iu], u[iu], u[iu], rf, rf, S)
    return nu


# ----
# Tensor plot functions. Private.
# ----
def _plot3D_mplib(R, X, Y, Z, scale_radius, u, utext, lattice):
    """
    Get 3D plot using matplotlib. Private.

    Args:
        R (numpy.ndarray): Computed elastic proeprties.
        X (numpy.ndarray): X coordinate mesh.
        Y (numpy.ndarray): Y coordinate mesh.
        Z (numpy.ndarray): Z coordinate mesh.
        scale_radius (bool): Whether to scale the spherical radius by R.
        u (numpy.ndarray): If not None, plot properties along specified vectors.
            nu\*3 array in the Cartesian coordinate.
        utext (list): nvu\*1 list of strings. Legends of vectors.
        lattice (numpy.ndarray): Lattice matrix. If not None, plot lattice
            around the 3D plot. The length of lattice vectors are arbitrary and
            no intersection with plot is ensured.

    Returns:
        fig (Figure): Matplotlib figure object
        ax (Axes): Matplotlib axes object
    """
    import matplotlib.pyplot as plt
    from matplotlib import animation, cm, colors
    from mpl_toolkits.mplot3d import Axes3D, axes3d

    rmin = np.min(R)
    rmax = np.max(R)
    norm = colors.Normalize(vmin=rmin, vmax=rmax, clip=False)
    fig = plt.figure(figsize=[5.5, 5])
    grid = plt.GridSpec(1, 11) # 5 for plot, 0.5 for colorbar
    ax = [fig.add_subplot(grid[0]),
          fig.add_subplot(grid[1:], projection='3d')]

    ax[1].plot_surface(
        X,
        Y,
        Z,
        rstride=1,
        cstride=1,
        facecolors=cm.jet(norm(R)),
        antialiased=True,
        alpha=0.75,
    )

    m = cm.ScalarMappable(cmap=cm.jet, norm=norm)
    m.set_array(np.linspace(rmin, rmax, 100))
    colorbar = plt.colorbar(m, cax=ax[0], shrink=0.7, location="left")
    colorbar.set_ticks(np.linspace(rmin, rmax, 5))

    # Plot 1D vectors
    if u is not None:
        for v, text in zip(u, utext):
            ax[1].plot([0, v[0]], [0, v[1]], [0, v[2]], color='k', linewidth=1)
            ax[1].text(v[0], v[1], v[2], text)
    else:
        u = np.zeros([1, 3])

    # Plot lattice
    if lattice is not None:
        if scale_radius == True:
            r = np.max(np.abs(R))
        else:
            r = 1
        lens = [np.linalg.norm(i) for i in lattice]
        shrink = 2*r / np.min(lens)
        plot_latt = lattice * shrink
        points = np.array([[-0.5, -0.5, -0.5], [0.5, -0.5, -0.5],
                           [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
                           [-0.5, -0.5, -0.5], [-0.5, -0.5, 0.5],
                           [0.5, -0.5, 0.5], [0.5, 0.5, 0.5],
                           [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5],
                           [0.5, -0.5, 0.5], [0.5, -0.5, -0.5],
                           [0.5, 0.5, -0.5], [0.5, 0.5, 0.5],
                           [-0.5, 0.5, 0.5], [-0.5, 0.5, -0.5]])
        points = np.matmul(points, plot_latt)
        ax[1].plot(points[:, 0], points[:, 1], points[:, 2], color='tab:gray', linewidth=1)
    else:
        points = np.zeros([1,3])

    # Make the planes transparent
    ax[1].xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax[1].yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax[1].zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # Make the grid lines transparent
    # ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    # ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    # ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    # Fixing limits
    if scale_radius == True:
        xmx = np.max(np.concatenate([points[:, 0], [rmax], u[:, 0]]))
        ymx = np.max(np.concatenate([points[:, 1], [rmax], u[:, 1]]))
        zmx = np.max(np.concatenate([points[:, 2], [rmax], u[:, 2]]))
    else:
        xmx = np.max(np.concatenate([points[:, 0], [1], u[:, 0]]))
        ymx = np.max(np.concatenate([points[:, 1], [1], u[:, 1]]))
        zmx = np.max(np.concatenate([points[:, 2], [1], u[:, 2]]))
    ax[1].set_xlim(-xmx, xmx)
    ax[1].set_ylim(-ymx, ymx)
    ax[1].set_zlim3d(-zmx, zmx)
    ax[1].locator_params(nbins=5)  # tight=True,
    ax[1].set_xlabel("X")
    ax[1].set_ylabel("Y")
    ax[1].set_zlabel("Z")
    # Fix aspect ratio
    # ax[1].set_box_aspect(aspect=(1, 1, 1))
    ax[1].set_aspect('equal')
    return fig, ax


def _plot3D_plotly(R, X, Y, Z, scale_radius, u, utext, lattice):
    """
    Get 3D plot using plotly. Private.

    Args:
        R (numpy.ndarray): Computed elastic proeprties.
        X (numpy.ndarray): X coordinate mesh.
        Y (numpy.ndarray): Y coordinate mesh.
        Z (numpy.ndarray): Z coordinate mesh.
        scale_radius (bool): Whether to scale the spherical radius by R.
        u (numpy.ndarray): If not None, plot properties along specified vectors.
            nu\*3 array in the Cartesian coordinate.
        utext (list): nvu\*1 list of strings. Legends of vectors.
        lattice (numpy.ndarray): Lattice matrix. If not None, plot lattice
            around the 3D plot. The length of lattice vectors are arbitrary and
            no intersection with plot is ensured.

    Returns:
        fig (Figure): Plotly Figure object.
    """
    import plotly.graph_objects as go

    rmin = np.min(R)
    rmax = np.max(R)

    surface = go.Surface(x=X, y=Y, z=Z, surfacecolor=R, colorscale='Jet')
    layout = go.Layout(
        title=''
    )
    fig = go.Figure(data=[surface], layout=layout)

    # Plot 1D vectors
    if u is not None:
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
    if lattice is not None:
        if scale_radius == True:
            r = np.max(np.abs(R))
        else:
            r = 1
        lens = [np.linalg.norm(i) for i in lattice]
        shrink = 2*r / np.min(lens)
        plot_latt = lattice * shrink
        points = np.array([[-0.5, -0.5, -0.5], [0.5, -0.5, -0.5],
                           [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
                           [-0.5, -0.5, -0.5], [-0.5, -0.5, 0.5],
                           [0.5, -0.5, 0.5], [0.5, 0.5, 0.5],
                           [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5],
                           [0.5, -0.5, 0.5], [0.5, -0.5, -0.5],
                           [0.5, 0.5, -0.5], [0.5, 0.5, 0.5],
                           [-0.5, 0.5, 0.5], [-0.5, 0.5, -0.5]])
        points = np.matmul(points, plot_latt)
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

    # Fixing limits
    if scale_radius == True:
        xmx = np.max(np.concatenate([points[:, 0], [rmax], u[:, 0]]))
        ymx = np.max(np.concatenate([points[:, 1], [rmax], u[:, 1]]))
        zmx = np.max(np.concatenate([points[:, 2], [rmax], u[:, 2]]))
    else:
        xmx = np.max(np.concatenate([points[:, 0], [1], u[:, 0]]))
        ymx = np.max(np.concatenate([points[:, 1], [1], u[:, 1]]))
        zmx = np.max(np.concatenate([points[:, 2], [1], u[:, 2]]))
    fig.update_layout(
        width=720,
        height=720,
        scene=dict(
            xaxis=dict(
                showbackground=True,
                backgroundcolor='white',
                gridcolor='lightgrey',
                zerolinecolor='black',
                range=[-xmx, xmx],),
            xaxis_title='X',
            yaxis=dict(
                showbackground=True,
                backgroundcolor='white',
                gridcolor='lightgrey',
                zerolinecolor='black',
                range=[-ymx, ymx],),
            yaxis_title='Y',
            zaxis=dict(
                showbackground=True,
                backgroundcolor='white',
                gridcolor='lightgrey',
                zerolinecolor='black',
                range=[-zmx, zmx],
            ),
            zaxis_title='Z',
        ),
        scene_aspectmode='manual',
        scene_aspectratio=dict(x=xmx, y=ymx, z=zmx),
    )
    return fig


def _plot2D_single(ax, theta, r, norm, color, title, rmax, u, utext, lattice):
    """
    Get a single 2D plot using matplotlib. Private.

    Args:
        ax (Axes): Matplotlib axes object.
        theta (numpy.ndarray): ntheta\*1 polar angles.
        r (numpy.ndarray): ntheta\*1 values.
        norm (numpy.ndarray): 3\*1 array of normalised Cartesian surface norm
            vector.
        color (str): Color of the loop.
        title (str): Subplot title.
        rmax (float): Maximum radius. 
        u (numpy.ndarray): If not None, plot properties along specified vectors.
            nu\*3 **normalized** array in the Cartesian coordinate.
        utext (list): nvu\*1 list of strings. Legends of vectors.
        lattice (numpy.ndarray): Lattice matrix. If not None, indicate the
            orientations of lattice base vectors in the 2D plot.

    Returns:
        ax (Axes): Matplotlib axes object
    """
    ax.plot(theta, r, color=color, linewidth=2)
    if rmax is None:
        rmax = np.max(r)

    # plot lattice vectors
    theta = np.arccos(norm[2])
    phi = np.arctan2(norm[1], norm[0])
    if lattice is not None:
        a = lattice[0, :] / np.linalg.norm(lattice[0, :])
        b = lattice[1, :] / np.linalg.norm(lattice[1, :])
        c = lattice[2, :] / np.linalg.norm(lattice[2, :])
        rmarker = rmax/20
        for v, vlabel in zip([a,b,c], ['a','b','c']):
            if np.abs(np.abs(np.dot(v, norm)) - 1) < 1e-4: # parallel
                if np.linalg.norm(v+norm) > 1e-4: # same direction, dot
                    ax.scatter(0, 0, c='tab:gray', marker='o')
                else: # opposite direction, cross
                    ax.scatter(0, 0, c='tab:gray', marker='x')
                ax.text(0, rmarker, vlabel, color='tab:gray',
                        fontstyle='italic', ha='center', va='center')
            else:
                # get the projection of base vector on plane
                proj = np.cross(norm, np.cross(v, norm))
                proj = proj / np.linalg.norm(proj)
                cossin_chi = np.matmul(
                    np.linalg.inv(np.array([[np.cos(phi)*np.cos(theta), -np.sin(phi)],
                                            [np.sin(phi)*np.cos(theta), np.cos(phi)]])),
                    [proj[0], proj[1]]
                )
                chi = np.arctan2(cossin_chi[1], cossin_chi[0])
                ax.vlines(chi, 0, rmax, colors='tab:gray', linewidth=1.5)
                ax.text(chi, rmax+rmarker, vlabel, color='tab:gray',
                        fontstyle='italic', ha='center', va='center')

    # plot 1D vectors
    if u is not None:
        for v, vlabel in zip(u, utext):
            if np.abs(np.dot(v, norm)) > 1e-4: # non-orthogonal
                continue
            cossin_chi = np.matmul(
                np.linalg.inv(np.array([[np.cos(phi)*np.cos(theta), -np.sin(phi)],
                                        [np.sin(phi)*np.cos(theta), np.cos(phi)]])),
                [v[0], v[1]]
            )
            chi = np.arctan2(cossin_chi[1], cossin_chi[0])
            ax.vlines(chi, 0, rmax, colors='k', linewidth=2)
            ax.text(chi, rmax*0.75, vlabel, color='k', va='center')

    ax.set_rmax(rmax)
    ax.set_theta_zero_location('E')
    ax.set_rlabel_position(0.0)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_yticks(np.round(np.linspace(0, rmax, 4), 2))
    ax.tick_params(axis='y', labelrotation=-45)
    ax.set_title(title, y=1.05)
    return ax





