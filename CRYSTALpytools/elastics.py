#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The module for elastic properties.
"""
from CRYSTALpytools import units
import numpy as np

# __all__ = ['Tensor3D', 'Tensor2D', 'Tensor1D',
#            'young', 'comp', 'shear', 'poisson',
#            'tensor_from_file', 'voigt2cart', 'cart2voigt']


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
        6\*6 compliance matrix in Voigt annotation. Unit: GPa:math:`^{-1}`.
        """
        return self._S

    @classmethod
    def from_file(cls, output):
        """
        Read elastic tensor from CRYSTAL output file and generate ``Tensor3D``
        object. Calls the ``crystal_io.Crystal_output.get_elatensor()`` method.
        Lattice information is obtained.

        Args:
            output (str): CRYSTAL output file.

        Returns:
            cls (Tensor3D)
        """
        from CRYSTAlpytools.crystal_io import Crystal_output

        out = Crystal_output(output)
        if out.get_dimensionality() != 3:
            raise Exception("Dimensionality error. Input system is not a 3D system.")

        return cls(matrix=out.get_elatensor(),
                   lattice=out.get_lattice(initial=False),
                   is_compliance=False)

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

        if use_cartesian != True and np.all(self.lattice==None):
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

        if use_cartesian != True and np.all(self.lattice==None):
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
            if np.abs(np.dot(u[i, :], v[i, :])) > 1e-6:
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
                uniform_scale=True, plot_lattice=True):
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
            ax (Axes): Matplotlib ``Axis`` object or a list of it if multiple
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

        if np.all(self.lattice==None):
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
        v1 = np.zeros([len(planeplt), ntheta+1, 3], dtype=float)
        for iplot, v0 in enumerate(planeplt):
            # get pole angles and 3D vectors. use (v0) v1
            theta = np.arccos(v0[2])
            phi = np.arctan2(v0[1], v0[0])
            chi = np.linspace(0, 2*np.pi, ntheta+1) # pole angle is actually chi.
            v1[iplot, :, :] = np.vstack([
                np.cos(phi) * np.cos(theta) * np.cos(chi) - np.sin(phi) * np.sin(chi),
                np.sin(phi) * np.cos(theta) * np.cos(chi) + np.cos(phi) * np.sin(chi),
                -np.sin(theta) * np.cos(chi)
            ]).transpose()

            # get pole angles and 3D vectors. use v1, v2
            if 'avg' in property or 'min' in property or 'max' in property:
                v2 = np.zeros([ntheta+1, nchi+1, 3], dtype=float) # another chi.
                chi2 = np.linspace(0, 2*np.pi, nchi+1)
                for i in range(ntheta+1):
                    theta2 = np.arccos(v1[iplot, i, 2])
                    phi2 = np.arctan2(v1[iplot, i, 1], v1[iplot, i, 0])
                    v2[i, :, :] = np.vstack([
                        np.cos(phi2) * np.cos(theta2) * np.cos(chi2) - np.sin(phi2) * np.sin(chi2),
                        np.sin(phi2) * np.cos(theta2) * np.cos(chi2) + np.cos(phi2) * np.sin(chi2),
                        -np.sin(theta2) * np.cos(chi2)
                    ]).transpose()

            # compute
            if property == 'young':
                r[iplot, :] = young(self._Scart, v1[iplot])
            elif property == 'comp':
                r[iplot, :] = comp(self._Scart, v1[iplot])
            elif property == 'shear':
                r[iplot, :] = shear(self._Scart, [v0], [v1[iplot]])
            elif property == 'shear avg':
                r[iplot, :] = np.average(shear(self._Scart, v1[iplot], v2), axis=1)
            elif property == 'shear min':
                r[iplot, :] = np.min(shear(self._Scart, v1[iplot], v2), axis=1)
            elif property == 'shear max':
                r[iplot, :] = np.max(shear(self._Scart, v1[iplot], v2), axis=1)
            elif property == 'poisson':
                r[iplot, :] = poisson(self._Scart, [v0], [v1[iplot]])
            elif property == 'poisson avg':
                r[iplot, :] = np.average(poisson(self._Scart, v1[iplot], v2), axis=1)
            elif property == 'poisson min':
                r[iplot, :] = np.min(poisson(self._Scart, v1[iplot], v2), axis=1)
            elif property == 'poisson max':
                r[iplot, :] = np.max(poisson(self._Scart, v1[iplot], v2), axis=1)

        # Plot
        # Set plot framework
        if np.all(layout!=None):
            if layout[0] + layout[1] != len(planeplt):
                warnings.warn('Inconsistent plane number and figure layout. Default layout is used.',
                              stacklevel=2)
                layout = None
        if np.all(layout==None):
            layout = [len(planeplt)//3+1, 0]
            if len(planeplt) < 3:
                layout[1] = len(planeplt)
            else:
                layout[1] = 3

        fig, ax = plt.subplots(nrows=layout[0], ncols=layout[1],
                               subplot_kw={'projection' : 'polar'}, layout='tight')
        if len(planeplt) == 1:
            ax = [ax]

        # Prepare 1D plot
        if np.all(u!=None):
            if not isinstance(u, str):
                utmp0 = np.array(u, dtype=float) # cartesian u, normalised
                if len(np.shape(utmp0)) == 1:
                    utmp0 = np.array([utmp0])

                if use_cartesian != True:
                    for i in range(len(utmp0)):
                        utmp0[i, :] = self.lattice.get_cartesian_coords(utmp0[i, :])
                        utmp0[i, :] = utmp0[i, :] / np.linalg.norm(utmp0[i, :])
                else:
                    for i in range(len(utmp0)):
                        utmp0[i, :] = utmp0[i, :] / np.linalg.norm(utmp0[i, :])

                # nplt*nu array
                utmp = np.zeros([len(planeplt), len(utmp0), 3], dtype=float)
                for i in range(len(planeplt)):
                    utmp[i, :, :] = utmp0
            else:
                use_cartesian = True
                if u.lower() == 'bothends':
                    utmp = np.zeros([len(planeplt), 2, 3], dtype=float)
                else:
                    utmp = np.zeros([len(planeplt), 1, 3], dtype=float)
                for iplot in range(len(planeplt)):
                    if u.lower() == 'max':
                        utmp[iplot, 0, :] = v1[iplot, np.argmax(r[iplot, :]), :] / np.linalg.norm(v1[iplot, np.argmax(r[iplot, :]), :])
                    elif u.lower() == 'min':
                        utmp[iplot, 0, :] = v1[iplot, np.argmin(r[iplot, :]), :] / np.linalg.norm(v1[iplot, np.argmin(r[iplot, :]), :])
                    elif u.lower() == 'bothends':
                        utmp[iplot, 0, :] = v1[iplot, np.argmax(r[iplot, :]), :] / np.linalg.norm(v1[iplot, np.argmax(r[iplot, :]), :])
                        utmp[iplot, 1, :] = v1[iplot, np.argmin(r[iplot, :]), :] / np.linalg.norm(v1[iplot, np.argmin(r[iplot, :]), :])
                    else:
                        raise ValueError("Unknown u value: '{}'.".format(u))

            if np.all(utext!=None):
                if isinstance(utext, str):
                    if utext.lower() == 'value':
                        utext = [['value' for i in range(np.shape(utmp)[1])] for j in range(np.shape(utmp)[0])]
                    else:
                        raise ValueError("Unknown utext value: '{}'.".format(utext))
                else:
                    utext = [utext for j in range(np.shape(utmp)[0])]
                    if len(utext) != np.shape(utmp)[1]:
                        warnings.warn('Length of vector annotations should be equal to length of vectors. Default annotations are used',
                                      stacklevel=2)
                        utext = None

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
            if np.all(u!=None):
                uplt = []
                utextplt = []
                ru = []
                for ivp, vp in enumerate(utmp[iplot]):
                    if np.abs(np.dot(vp, v0)) > 1e-6: # non-orthogonal
                        continue
                    uplt.append(vp)
                    if property == 'shear' or property == 'poisson':
                        ru.append(self.get_pair(property, v0, vp, use_cartesian=True))
                    else:
                        ru.append(self.get_1D(property, vp, nchi, use_cartesian=True))

                    if np.all(utext!=None):
                        if utext[iplot][ivp].lower() == 'value':
                            utextplt.append('{:.2f}'.format(ru[-1]))
                        else:
                            utextplt.append(utext[iplot][ivp])
                    else:
                        if not isinstance(u, str):
                            if len(vp) == 1:
                                u = np.array([u], dtype=float)
                            else:
                                u = np.array(u, dtype=float)
                        else:
                            u = utmp[iplot]
                        utextplt.append('[{:<4.1f}{:^4.1f}{:>4.1f}]  {:.2f}'.format(
                            u[ivp, 0], u[ivp, 1], u[ivp, 2], ru[-1]
                        ))
                uplt = np.array(uplt, dtype=float)
            else:
                uplt=None; utextplt=[]

            # plot figures
            if add_title == True:
                if plane_definition == 'miller':
                    title = '{}\n({:<2d} {:^2d} {:>2d})'.format(
                        property, plane[iplot][0], plane[iplot][1], plane[iplot][2])
                elif plane_definition == 'cartesian':
                    title = '{}\nCartesian Norm: ({:<4.1f} {:^4.1f} {:>4.1f})'.format(
                        property, plane[iplot][0], plane[iplot][1], plane[iplot][2])
                else:
                    title = '{}\nFractional Norm: ({:<4.1f} {:^4.1f} {:>4.1f})'.format(
                        property, plane[iplot][0], plane[iplot][1], plane[iplot][2])
            else:
                title = ''

            ax.flat[iplot] = _plot2D_single(
                    ax.flat[iplot], chi, r[iplot, :], v0, pcolor[iplot%10],
                    title, rmax, uplt, utextplt, platt
            )

        if len(planeplt) == 1:
            ax = ax[0]

        return fig, ax

    def plot_3D(self, property, nphi=90, ntheta=90, nchi=180,
                scale_radius=True, range_cbar=None, range_x=None, range_y=None,
                range_z=None, u=None, utext=None, use_cartesian=True,
                plot_lattice=False, plot_lib='matplotlib', **kwargs):
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
                For 'bothends', max first.
            utext (list[str] | str): A string or a list of string. Used to mark
                The vector. If ``None``, the input ``u`` and the corresponding
                value is annotated. If 'value', the value is annotated.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)
            plot_lattice (bool): Draw the lattice box around the 3D surface.
            plot_lib (bool): 'matplotlib' or 'plotly'. By default use
                `matplotlib <https://matplotlib.org/>`_. Alternatively use
                `plotly <https://plotly.com/>`_ to enable interactive
                inspections in Jupyter Notebook.
            \*\*kwargs: Parameters passed to ``Axes3D.view_init`` (matplotlib)
                or ``fig.update_layout()`` (plotly). Only camera position
                keywords are suggested.

        Returns:
            fig (Figure): Matplotlib or plotly ``Figure`` object.
            ax (list[Axes] | None): 2\*1 list of matplotlib ``Axis`` object. The
                first element is the axis of colorbar. The second is the axis of
                3D surface. ``None`` if ``plot_lib = 'plotly'``.

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
        if np.all(u!=None):
            if isinstance(u, str):
                if u.lower() == 'max':
                    irow = np.argmax(R)//np.shape(R)[1]
                    icol = np.argmax(R)%np.shape(R)[1]
                    uplt = np.array([
                        X[irow, icol], Y[irow, icol], Z[irow, icol],
                    ], dtype=float)
                elif u.lower() == 'min':
                    irow = np.argmin(R)//np.shape(R)[1]
                    icol = np.argmin(R)%np.shape(R)[1]
                    uplt = np.array([
                        X[irow, icol], Y[irow, icol], Z[irow, icol],
                    ], dtype=float)
                elif u.lower() == 'bothends':
                    irow1 = np.argmax(R)//np.shape(R)[1]
                    icol1 = np.argmax(R)%np.shape(R)[1]
                    irow2 = np.argmin(R)//np.shape(R)[1]
                    icol2 = np.argmin(R)%np.shape(R)[1]
                    uplt = np.array([
                        [X[irow1, icol1], Y[irow1, icol1], Z[irow1, icol1]],
                        [X[irow2, icol2], Y[irow2, icol2], Z[irow2, icol2]],
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
            if np.all(utext==None):
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
        if np.all(plot_lib!=None):
            if plot_lib.lower() == 'matplotlib':
                fig, ax = _plot3D_mplib(R, X, Y, Z, scale_radius, uplt, utext, platt,
                                        range_cbar, range_x, range_y, range_z, **kwargs)
            else:
                fig = _plot3D_plotly(R, X, Y, Z, scale_radius, uplt, utext, platt,
                                     range_cbar, range_x, range_y, range_z, **kwargs)
                ax = None
            return fig, ax
        else: # Developer settings. return all the variables required by matplotlib.
            return R, X, Y, Z, scale_radius, uplt, utext, platt

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


class Tensor2D(Tensor3D):
    """
    2D elastic tensor and related properties. Periodic boundary conditions are
    consistent with CRYSTAL definitions, i.e., xy-periodic and slab vertical to
    z.

    Args:
        matrix (numpy.ndarray): 3\*3 compliance or stiffness matrix. Unit: GPa:math:`^{-1}` or GPa.
        lattice (numpy.ndarray): 2\*2 lattice matrix.
        is_compliance (bool): Compliance or stiffness matrix used as input.
    """
    def __init__(self, matrix, lattice=None, is_compliance=True):
        import warnings
        from pymatgen.core.lattice import Lattice

        if np.all(lattice==None):
            warnings.warn("Lattice information not available, fractional coordinates disabled.",
                          stacklevel=2)
            self.lattice = None
        else:
            latt = np.eye(3)*500
            latt[0:2, 0:2] = lattice
            self.lattice = Lattice(latt)

        mx = np.zeros([6,6], dtype=float)
        mx[0:2, 0:2] = matrix[0:2, 0:2]
        mx[0, 5] = matrix[0, 2]
        mx[5, 0] = mx[0, 5]
        mx[1, 5] = matrix[1, 2]
        mx[5, 1] = mx[1, 5]
        mx[5, 5] = matrix[2, 2]
        if is_compliance == True:
            self._S = np.array(mx, dtype=float)
            self._C = np.linalg.inv(mx)
        else:
            self._C = np.array(mx, dtype=float)
            self._S = np.linalg.inv(mx)

        # convert matrices from Voigt notation to Cartesian notation
        self._Ccart = voigt2cart(self._C)
        self._Scart = voigt2cart(self._S)

    @property
    def stiffness(self):
        """
        3\*3 stiffness matrix in Voigt annotation. Unit: GPa.
        """
        C = np.zeros([3,3], dtype=float)
        C[0:2, 0:2] = self._C[0:2, 0:2]
        C[0, 2] = self._C[0, 5]
        C[2, 0] = C[0, 2]
        C[1, 2] = self._C[1, 5]
        C[2, 1] = C[1, 2]
        C[2, 2] = self._C[5, 5]
        return C

    @property
    def compliance(self):
        """
        6\*6 compliance matrix in Voigt annotation. Unit: GPa:math:`^{-1}`.
        """
        S = np.zeros([3,3], dtype=float)
        S[0:2, 0:2] = self._S[0:2, 0:2]
        S[0, 2] = self._S[0, 5]
        S[2, 0] = S[0, 2]
        S[1, 2] = self._S[1, 5]
        S[2, 1] = S[1, 2]
        S[2, 2] = self._S[5, 5]
        return S

    @classmethod
    def from_file(cls, output):
        """
        Read elastic tensor from CRYSTAL output file and generate ``Tensor2D``
        object. Calls the ``crystal_io.Crystal_output.get_elatensor()`` method.
        Lattice information is obtained.

        Args:
            output (str): CRYSTAL output file.

        Returns:
            cls (Tensor2D)
        """
        from CRYSTAlpytools.crystal_io import Crystal_output

        out = Crystal_output(output)
        if out.get_dimensionality() != 2:
            raise Exception("Dimensionality error. Input system is not a 3D system.")

        return cls(matrix=out.get_elatensor(),
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
            e (float | numpy.ndarray): Some sort of property.
        """
        property_list = ['young', 'comp', 'shear', 'poisson']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        if use_cartesian != True and np.all(self.lattice==None):
            raise Exception('Lattice matrix not available.')

        u = np.array(u, dtype=float)
        if len(np.shape(u)) == 1:
            u = np.array([u])

        u = np.hstack([u, [[0] for i in range(len(u))]])
        if use_cartesian != True:
            for i in range(len(u)):
                u[i, :] = self.lattice.get_cartesian_coords(u[i, :])
                u[i, :] = u[i, :] / np.linalg.norm(u[i, :])
        else:
            for i in range(len(u)):
                u[i, :] = u[i, :] / np.linalg.norm(u[i, :])

        # define the norm vector of u
        if 'shear' in property or 'poisson' in property:
            v = np.cross(u, [0, 0, 1])
            for iv in range(len(v)):
                v[iv, :] = v[iv, :] / np.linalg.norm(v[iv, :])

        if property == 'young':
            e = young(self._Scart, u)
        elif property == 'comp':
            e = comp(self._Scart, u)
        elif property == 'shear':
            e = shear(self._Scart, u, v)[:, 0]
        elif property == 'poisson':
            e = poisson(self._Scart, u, v)[:, 0]

        if len(u) == 1:
            e = e[0]
        return e

    def plot_2D(self, property, ntheta=90, u=None, utext=None,
                use_cartesian=True, add_title=True, uniform_scale=True,
                plot_lattice=True):
        """
        Plot 2D crystal elastic properties of the system.

        Options:

        * "young": Young's modulus in GPa.
        * "comp": Linear compressibility.
        * "shear": Shear modulus between plane norm and vectors in the plane, in GPa.
        * "shear plane": Shear modulus between 2 vectors in the plane, in GPa.
        * "poisson": Poisson ratio between plane norm and vector in the plane.
        * "poisson plane": Poisson ratio between 2 vectors in the plane.

        Args:
            property (str): Property to compute. See options above.
            ntheta (int): Resolution of polar angle :math:`\\theta` on plane,
                in radian :math:`[0, 2\\pi)`.
            u (numpy.ndarray): 2\*1 or nu\*2 array of vectors to be added into
                the figure. A line segment to the limit of axis is plotted to
                indicate the vector.
            utext (list[str] | str): A string or a list of string. Used to mark
                The vector. If ``None``, the input ``u`` and the corresponding
                value is annotated. If 'value', the value is annotated.
            use_cartesian (bool): Vector is defined as cartesian or fractional
                coordinates. (*Only when lattice information is available*.)
            add_title (bool): Add title to subplot. Only entries of ``plane``
                can be used.
            uniform_scale (bool): Use the same radial scale (elastic properties)
                for all the plots.
            plot_lattice (bool): Draw lattice base vectors to indicate
                orientation of the plane.

        Returns:
            fig (Figure): Matplotlib ``Figure`` object.
            ax (Axes): Matplotlib ``Axis`` object.
        """
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        import warnings

        # sanity check
        property_list = ['young', 'comp', 'shear', 'shear plane', 'poisson', 'poisson plane']
        property = property.lower()
        if property not in property_list:
            raise ValueError("Unknown property input: '{}'".format(property))

        if np.all(self.lattice==None):
            if use_cartesian != True or plot_lattice == True:
                raise Exception('Lattice matrix not available.')

        # plane definition: use normalised cartesian for plot
        planeplt = np.array([0, 0, 1], dtype=float)
        v0 = planeplt

        # Compute and plot the 2D plane
        # get pole angles and 3D vectors. use (v0) v1
        theta = 0.
        phi = 0.
        chi = np.linspace(0, 2*np.pi, ntheta+1) # pole angle is actually chi.
        v1 = np.zeros([ntheta+1, 3], dtype=float)
        v1[:, :] = np.vstack([
            np.cos(phi) * np.cos(theta) * np.cos(chi) - np.sin(phi) * np.sin(chi),
            np.sin(phi) * np.cos(theta) * np.cos(chi) + np.cos(phi) * np.sin(chi),
            -np.sin(theta) * np.cos(chi)
        ]).transpose()

        # get pole angles and 3D vectors. use v1, v2
        if 'plane' in property:
            v2 = np.zeros([1, ntheta+1, 3], dtype=float) # another chi.
            v2[0, :, :] = np.cross(v1, [0, 0, 1])
            for i in range(ntheta+1):
                v2[0, i, :] = v2[0, i, :] / np.linspace(v2[0, i, :])
            v2 = np.transpose(v2, axes=[1, 0, 2])

        # compute
        if property == 'young':
            r = young(self._Scart, v1)
        elif property == 'comp':
            r = comp(self._Scart, v1)
        elif property == 'shear':
            r = shear(self._Scart, [v0], [v1])
        elif property == 'shear plane':
            r = shear(self._Scart, v1, v2)
            r = r.transpose()[0]
        elif property == 'poisson':
            r = poisson(self._Scart, [v0], [v1])
        elif property == 'poisson plane':
            r = poisson(self._Scart, v1, v2)
            r = r.transpose()[0]

        # Plot
        fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection' : 'polar'})

        # Prepare 1D plot
        if np.all(u!=None):
            if not isinstance(u, str):
                utmp = np.array(u, dtype=float) # cartesian u, normalised
                if len(np.shape(uplt)) == 1:
                    utmp = np.array([utmp])
                if np.shape(utmp)[1] != 2:
                    raise ValueError('Vector u must have only 2 dimensions.')

                uplt = np.zeros([np.shape(utmp)[0], 3], dtype=float)
                uplt[:, 0:2] = utmp
                if use_cartesian != True:
                    for i in range(len(uplt)):
                        uplt[i, :] = self.lattice.get_cartesian_coords(uplt[i, :])
                        uplt[i, :] = uplt[i, :] / np.linalg.norm(uplt[i, :])
                else:
                    for i in range(len(uplt)):
                        uplt[i, :] = uplt[i, :] / np.linalg.norm(uplt[i, :])
            else:
                use_cartesian = True
                if u.lower() == 'max':
                    uplt = np.zeros([1, 3], dtype=float)
                    uplt[0, :] = v1[np.argmax(r), :] / np.linalg.norm(v1[np.argmax(r), :])
                elif u.lower() == 'min':
                    uplt = np.zeros([1, 3], dtype=float)
                    uplt[0, :] = v1[np.argmin(r), :] / np.linalg.norm(v1[np.argmin(r), :])
                elif u.lower() == 'bothends':
                    uplt = np.zeros([2, 3], dtype=float)
                    uplt[0, :] = v1[np.argmax(r), :] / np.linalg.norm(v1[np.argmax(r), :])
                    uplt[1, :] = v1[np.argmin(r), :] / np.linalg.norm(v1[np.argmin(r), :])
                else:
                    raise ValueError("Unknown u value: '{}'.".format(u))

            utextplt = []
            for ivp, vp in enumerate(uplt):
                if property == 'shear' or property == 'poisson':
                    ru = super().get_pair(property, v0, vp, use_cartesian=True)
                else:
                    ru = self.get_1D(property, vp, use_cartesian=True)

                if np.all(utext!=None):
                    if isinstance(utext, str):
                        if utext.lower() == 'value':
                            utextplt.append('{:.2f}'.format(ru))
                        else:
                            raise ValueError("Unknown utext value: '{}'.".format(utext))
                    else:
                        utextplt.append(utext[ivp])
                else:
                    if not isinstance(u, str):
                        if len(vp) == 1:
                            u = np.array([u], dtype=float)
                        else:
                            u = np.array(u, dtype=float)
                    else:
                        u = uplt
                    utextplt.append('[{:<4.1f}{:^4.1f}{:>4.1f}]  {:.2f}'.format(
                        u[ivp, 0], u[ivp, 1], u[ivp, 2], ru
                    ))
        else:
            uplt=None; utextplt=[]

        # Prepare lattice
        if plot_lattice == False:
            platt = None
        else:
            platt = self.lattice.matrix

        pcolor = list(mcolors.TABLEAU_COLORS.keys())

        # plot figures
        if add_title == True:
            title = property
        else:
            title = ''

        rmax = np.max(r)
        ax = _plot2D_single(ax, chi, r, v0, pcolor[0], title, rmax, uplt, utextplt, platt)
        return fig, ax

# class Tensor1D():


    
# ----
# Basic elastic property functions.
# ----
def young(S, u):
    """
    The equation of Young's modulus. Using ``Tensor3D`` object is recommended.

    .. math::

        E = \\frac{1}{u_{i}u_{j}u_{k}u_{l}S_{ijkl}}

    Einstein's notation is used. :math:`i,j,k,l=1,3; j>i; l>k`.

    Args:
        S (numpy.ndarray): 3\*3\*3\*3 compliance matrix in Cartesian form.
            Unit: GPa:math:`^{-1}`.
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
    The equation of linear compressibility. Using ``Tensor3D`` object is recommended.

    .. math::

        \\beta = u_{i}u_{j}S_{ijkk}

    Einstein's notation is used. :math:`i,j,k=1,3; j>i`.

    Args:
        S (numpy.ndarray): 3\*3\*3\*3 compliance matrix in Cartesian form.
            Unit: GPa:math:`^{-1}`
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
    The equation of shear modulus. Using ``Tensor3D`` object is recommended.

    .. math::

        G = \\frac{1}{u_{i}v_{j}u_{k}v_{l}S_{ijkl}}

    Einstein's notation is used. :math:`i,j,k,l=1,3`.

    Args:
        S (numpy.ndarray): 3\*3\*3\*3 compliance matrix in Cartesian form.
            Unit: GPa:math:`^{-1}`.
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
    The equation of Poisson ratio. Using ``Tensor3D`` object is recommended.

    .. math::

        \\nu = -\\frac{u_{i}u_{j}v_{k}v_{l}S_{ijkl}}{u_{i}u_{j}u_{k}u_{l}S_{ijkl}}

    Einstein's notation is used. :math:`i,j,k=1,3; j>i`.

    Args:
        S (numpy.ndarray): 3\*3\*3\*3 compliance matrix in Cartesian form.
            Unit: GPa:math:`^{-1}`.
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
# auxiliary functions
# ----
# def tensor_from_file():


def voigt2cart(V):
    """
    Convert 6\*6 stiffness / compliance matrix in Voigt representation into
    3\*3\*3\*3 Cartesian representation.

    Args:
        V (numpy.ndarray): Voigt represented matrix.

    Returns:
        C (numpy.ndarray): Cartesian represented matrix.
    """
    voigt = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]], dtype=int)
    C = np.zeros([3,3,3,3], dtype=float)

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
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
    voigt = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]], dtype=int)
    V = np.zeros([6,6], dtype=float)

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    V[voigt[i, j], voigt[k,l]] = C[i, j, k, l]
    return V


# ----
# Tensor plot functions. Private.
# ----
def _plot3D_mplib(R, X, Y, Z, scale_radius, u, utext, lattice, range_cbar,
                  range_x, range_y, range_z, Rlatt=None, **kwargs):
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
        range_cbar, range_x, range_y, range_z (list[float,float]): *Not
            suggested* Explicitly specifying the ranges of colorbar, x, y and z
            axes.
        \*\*kwargs: Parameters passed to ``Axes3D.view_init`` Only camera
            position keywords are suggested.s
        Rlatt (numpy.ndarray): *Developers Only* Auxiliary data set to plot
            lattices of the same scale.

    Returns:
        fig (Figure): Matplotlib figure object
        ax (Axes): Matplotlib axes object
    """
    import matplotlib.pyplot as plt
    from matplotlib import animation, cm, colors
    from mpl_toolkits.mplot3d import Axes3D, axes3d

    fig = plt.figure(figsize=[5.5, 5], layout='tight')
    grid = plt.GridSpec(1, 12) # 5 for plot, 0.5 for colorbar, 0.5 for gap
    ax = [fig.add_subplot(grid[0]), fig.add_subplot(grid[2:], projection='3d')]

    if np.all(range_cbar==None):
        rmin = np.min(R)
        rmax = np.max(R)
    else:
        rmin = np.min(range_cbar[0])
        rmax = np.max(range_cbar[1])
    norm = colors.Normalize(vmin=rmin, vmax=rmax, clip=False)

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
    colorbar = plt.colorbar(m, cax=ax[0], shrink=0.5, location="left")
    # colorbar.set_ticks(np.linspace(rmin, rmax, 5))

    # Plot 1D vectors
    if np.all(u!=None):
        for v, text in zip(u, utext):
            ax[1].plot([0, v[0]], [0, v[1]], [0, v[2]], color='k', linewidth=1)
            ax[1].text(v[0], v[1], v[2], text)
    else:
        u = np.zeros([1, 3])

    # Plot lattice
    if np.all(lattice!=None):
        if scale_radius == True:
            if np.all(Rlatt==None):
                r = np.max(np.abs(R))
            else:
                r = np.max(np.abs(Rlatt))
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
        points = points @ plot_latt
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

    ax[1].set_xlim(np.min(range_x), np.max(range_x))
    ax[1].set_ylim(np.min(range_y), np.max(range_y))
    ax[1].set_zlim3d(np.min(range_z), np.max(range_z))

    ax[1].locator_params(nbins=5)  # tight=True,
    ax[1].set_xlabel("X")
    ax[1].set_ylabel("Y")
    ax[1].set_zlabel("Z")
    # Fix aspect ratio
    ax[1].set_aspect('equal')
    if len(kwargs.keys()) > 0:
        ax[1].view_init(**kwargs)
    if np.all(fig!=None):
        ax = fig.axes

    return fig, ax


def _plot3D_plotly(R, X, Y, Z, scale_radius, u, utext, lattice, range_cbar,
                   range_x, range_y, range_z, **kwargs):
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
        range_cbar, range_x, range_y, range_z (list[float,float]): *Not
            suggested* Explicitly specifying the ranges of colorbar, x, y and z
            axes.
        \*\*kwargs: Parameters passed to ``fig.update_layout()`` Only camera
            position keywords are suggested.

    Returns:
        fig (Figure): Plotly Figure object.
    """
    import plotly.graph_objects as go

    if np.all(range_cbar==None):
        rmin = np.min(R)
        rmax = np.max(R)
    else:
        rmin = np.min(range_cbar[0])
        rmax = np.max(range_cbar[1])

    surface = go.Surface(x=X, y=Y, z=Z, surfacecolor=R, colorscale='Jet',
                         cmin=rmin,cmax=rmax)
    layout = go.Layout(
        title=''
    )
    fig = go.Figure(data=[surface], layout=layout)

    # Plot 1D vectors
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
    if np.all(lattice!=None):
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
    if scale_radius == True:
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
        width=720,
        height=720,
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
    if np.all(rmax==None):
        rmax = np.max(r)

    # plot lattice vectors
    theta = np.arccos(norm[2])
    phi = np.arctan2(norm[1], norm[0])
    if np.all(lattice!=None):
        a = lattice[0, :] / np.linalg.norm(lattice[0, :])
        b = lattice[1, :] / np.linalg.norm(lattice[1, :])
        c = lattice[2, :] / np.linalg.norm(lattice[2, :])
        rmarker = rmax/20
        for v, vlabel in zip([a,b,c], ['a','b','c']):
            if np.abs(np.abs(np.dot(v, norm)) - 1) < 1e-6: # parallel
                if np.linalg.norm(v+norm) > 1e-6: # same direction, dot
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
    if np.all(u!=None):
        for v, vlabel in zip(u, utext):
            if np.abs(np.dot(v, norm)) > 1e-6: # non-orthogonal
                continue
            cossin_chi = np.matmul(
                np.linalg.inv(np.array([[np.cos(phi)*np.cos(theta), -np.sin(phi)],
                                        [np.sin(phi)*np.cos(theta), np.cos(phi)]])),
                [v[0], v[1]]
            )
            chi = np.arctan2(cossin_chi[1], cossin_chi[0])
            ax.vlines(chi, 0, rmax, colors='k', linewidth=2)
            ax.text(chi, rmax*1.1, vlabel, color='k', ha='center', va='center')

    ax.set_rmax(rmax)
    ax.set_theta_zero_location('E')
    ax.set_rlabel_position(0.0)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_yticks(np.round(np.linspace(0, rmax, 4), 2))
    ax.tick_params(axis='y', labelrotation=-45)
    ax.set_title(title, y=1.05)
    return ax





