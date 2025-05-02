#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods to parse files used by `XCrySDen <http://www.xcrysden.org>`_.
Developed for  XCrySDen <= 1.6.2, the latest version by 2024.
"""
class XSF():
    """
    A collection of methods to parse and write XCrySDen XSF files.

    .. note::

        At most 1 grid data is saved with this class.

    Args:
        struc (CStructure): Extended Pymatgen Structure object.
        forces (array): nForce\*4 array of forces (or any vector associated
            with atoms). The first element is the index of atom (from 0) and
            the last 3 are Cartesian components along x, y and z directions.
            Unit: Hartree/:math:`\\AA`.
        grid_base (array): 3(4)\*3 array of Cartesian coordinates of points O,
            A, B(, C) to define a 2D(3D) grid. Vectors OA, OB (and OC) are used.
            Unit: :math:`\\AA`.
        grid_data (arary): (nZ\*)nY\*nX array of 2D(3D) data grid.
    """
    def __init__(self, struc, forces=[], grid_base=[], grid_data=[]):
        import numpy as np
        import warnings

        self.structure = struc
        self.ndim = struc.lattice.pbc.count(True)

        if len(forces) > 0:
            self.forces = np.array(forces, ndmin=2, dtype=float)
            if self.forces.shape[1] != 4:
                raise ValueError("Input forces must be a nForce*4 array. The first element is the index of atom and the others are Cartesian components.")
        else:
            self.forces = []

        if len(grid_base) > 0:
            self.grid_base = np.array(grid_base, ndmin=2, dtype=float)
            self.grid_ndim = self.grid_base.shape[0]-1
        else:
            self.grid_ndim = 0
            self.grid_base = []

        if len(grid_data) > 0:
            self.grid_data = np.array(grid_data, ndmin=2, dtype=float)
            if self.grid_ndim != self.grid_data.ndim or np.any(self.grid_data.shape==1):
                raise ValueError("Inconsistent dimensionalities of grid base vector and data grid")
        else:
            self.grid_data = []
            if self.grid_ndim != 0:
                warnings.warn("Empty grid data! Grid base vectors are removed.", stacklevel=2)
                self.grid_ndim = 0
                self.grid_base = []

    @classmethod
    def read(cls, filename, grid_index=1):
        """
        Generate a XSF class from XSF file. Only the grid specified by
        ``grid_index`` (starting from 1) is saved.

        Args:
            filename (str): Iutput name.
            grid_index (int): The index of grid to be saved, if grid exists.
                Otherwise return to a class with geometry only.
        Returns:
            cls (XSF)
        """
        import os, re
        import pandas as pd
        import numpy as np
        from CRYSTALpytools.geometry import CStructure

        if not os.path.exists(filename):
            raise FileNotFoundError("Input file '{}' does not exist.".format(filename))

        df = pd.DataFrame(open(filename, 'r'))
        # remove empty lines and comments
        df = df.replace(r'^\s*#*\s*$', np.nan, regex=True).dropna()
        df.reset_index(drop=True, inplace=True)
        # geometry keyword
        if len(df[df[0].str.contains(r'^\s*ATOMS\s*$')].index) == 1:
            ndim = 0; read_cell = False
        elif len(df[df[0].str.contains(r'^\s*MOLECULE\s*$')].index) == 1:
            ndim = 0; read_cell = True
        elif len(df[df[0].str.contains(r'^\s*POLYMER\s*$')].index) == 1:
            ndim = 1; read_cell = True
        elif len(df[df[0].str.contains(r'^\s*SLAB\s*$')].index) == 1:
            ndim = 2; read_cell = True
        elif len(df[df[0].str.contains(r'^\s*CRYSTAL\s*$')].index) == 1:
            ndim = 3; read_cell = True
        else:
            raise Exception("No or more than 1 geometry found in '{}'.".format(filename))
        # read cell
        if read_cell == True:
            iat = df[df[0].str.contains(r'^\s*PRIMCOORD\s*$')].index
            icel = df[df[0].str.contains(r'^\s*PRIMVEC\s*$')].index
            if len(iat) != 1 or len(icel) != 1:
                raise Exception("No or more than 1 geometry is found in '{}'.".format(filename))
            latt = df[0].loc[icel[0]+1:icel[0]+3].map(lambda x: x.strip().split()).tolist()
            latt = np.array(latt, dtype=float)
        else:
            iat = df[df[0].str.contains(r'^\s*ATOMS\s*$')].index
            latt = np.eye(3) * 500.
        # read atom and force
        natom = int(df[0].loc[iat[0]+1].strip().split()[0])
        species = []; coords = []; forces = []
        for i in range(natom):
            ldata = df[0].loc[iat[0]+2+i].strip().split()
            species.append(ldata[0])
            coords.append(ldata[1:4])
            if len(ldata) == 7:
                forces.append([i+1, ldata[4], ldata[5], ldata[6]])
        if re.match(r'[0-9]+', ldata[0]):
            species = np.array(species, dtype=int)
        else:
            species = np.array(species, dtype=str)
        coords = np.array(coords, dtype=float)
        if len(forces) > 0:
            forces = np.array(forces, dtype=float)
        # get structure
        pbc = {0 : (False, False, False),
               1 : (True, False, False),
               2 : (True, True, False),
               3 : (True, True, True)}
        struc = CStructure(
            latt, species, coords, pbc=pbc[ndim], standarize=False,
            coords_are_cartesian=True
        )

        # grid keyword
        bgrd = df[df[0].str.contains(r'^\s*BEGIN_DATAGRID')].index
        egrd = df[df[0].str.contains(r'^\s*END_DATAGRID')].index
        if len(bgrd) != len(egrd):
            raise Exception("Data grid 'BEGIN' and 'END' keywords are not paired. File might have been broken.")
        if len(bgrd) > 0:
            if grid_index-1 >= len(bgrd):
                raise Exception("Only {:d} grids are found. Index {:d} exceeds that limit.".format(len(bgrd), grid_index))

            bgrd = bgrd[grid_index-1]; egrd = egrd[grid_index-1]
            gsize = np.array(df[0].loc[bgrd+1].strip().split(), dtype=int)
            gdim = len(gsize)
            gbase = df[0].loc[bgrd+2:bgrd+2+gdim].map(lambda x: x.strip().split()).tolist()
            gbase = np.array(gbase, dtype=float)
            for i in range(1, gdim+1):
                gbase[i] = gbase[i] + gbase[0]
            ndata = np.prod(gsize)
            nperline = len(df[0].loc[bgrd+3+gdim].strip().split())
            nline = ndata // nperline
            gdata = df[0].loc[bgrd+3+gdim:bgrd+2+gdim+nline].map(
                lambda x: x.strip().split()
            ).tolist()
            gdata = np.array(gdata, dtype=float).flatten()
            if ndata % nperline != 0:
                lastl = df[0].loc[egrd-1].strip().split()
                gdata = np.hstack([gdata, np.array(lastl, dtype=float)])
            gdata = gdata.reshape(gsize[::-1])
        else:
            gbase = []; gdata = []

        return cls(struc, forces=forces, grid_base=gbase, grid_data=gdata)

    def write(self, filename, geomonly=False, grid_name='UNKNOWN'):
        """
        Write data into a new XSF file.

        Args:
            filename (str): Output name.
            geomoly (bool): Writes geometry and forces (if present) only into XSF file.
            grid_name (str): Name of the grid, valid only if ``self.grid_ndim``
                is not 0.
        """
        if geomonly == True:
            obj = XSF(self.structure)
            obj._write_geom(filename)
        else:
            self._write_geom(filename)
            if self.grid_ndim > 0:
                self._write_grid(filename, grid_name)
        return

    def _write_geom(self, filename):
        """
        Write geometry and forces into a new XSF file.

        Args:
            filename (str): Output name.
        """
        import numpy as np
        import os, warnings

        if os.path.exists(filename):
            warnings.warn(
                "File '{}' exists! It will be overwritten.".format(filename),
                stacklevel=3)
        file = open(filename, 'w')

        # write geometry
        header = '# Generated by CRYSTALpytools\n'
        dimen_key = {3 : 'CRYSTAL', 2 : 'SLAB', 1 : 'POLYMER', 0 : 'MOLECULE'}
        if self.ndim > 0:
            header += ' %s\n' % dimen_key[self.ndim]
            header += ' %s\n' % 'PRIMVEC'
            lattmx = self.structure.lattice.matrix
            for i in range(3):
                header += ' %15.9f%15.9f%15.9f\n' % (lattmx[i,0], lattmx[i,1], lattmx[i,2])
            header += ' %s\n' % 'PRIMCOORD'
            header += ' %10i%10i\n' % (self.structure.num_sites, 1)
        else:
            header += ' %s\n' % 'ATOMS'

        if len(self.forces)>0:
            idx = self.forces[:,0]
            for i in range(self.structure.num_sites):
                fidx = np.where(idx-i==0)[0]
                if len(fidx) > 0:
                    header += ' %-4s%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f\n' % \
                          (self.structure.species_symbol[i], self.structure.cart_coords[i,0],
                           self.structure.cart_coords[i,1], self.structure.cart_coords[i,2],
                           self.forces[fidx[0], 1], self.forces[fidx[0], 2], self.forces[fidx[0], 3])
                else:
                    header += ' %-4s%15.9f%15.9f%15.9f\n' % \
                          (self.structure.species_symbol[i], self.structure.cart_coords[i,0],
                           self.structure.cart_coords[i,1], self.structure.cart_coords[i,2])
        else:
            for i in range(self.structure.num_sites):
                header += ' %-4s%15.9f%15.9f%15.9f\n' % \
                          (self.structure.species_symbol[i], self.structure.cart_coords[i,0],
                           self.structure.cart_coords[i,1], self.structure.cart_coords[i,2])
        file.write("%s\n" % header)
        file.close()
        return

    def _write_grid(self, filename, grid_name='UNKNOWN'):
        """
        Append scalar field data into XSF file. For every 'block', only 1 data
        grid is written.

        Args:
            filename (str): Output name.
            grid_name (str): Name of the grid.
        """
        import numpy as np

        file = open(filename, 'r+')
        header = file.read()
        file.close()

        # write 3D data
        if self.grid_ndim == 3:
            header += ' %s\n' % 'BEGIN_BLOCK_DATAGRID_3D'
            header += '   %s\n' % grid_name
            header += '   %s_%s\n' % ('BEGIN_DATAGRID_3D', grid_name)
            header += '   %8i%8i%8i\n' % (self.grid_data.shape[2],
                                          self.grid_data.shape[1],
                                          self.grid_data.shape[0]) # [[nX] nY] nZ
            header += '   %15.9f%15.9f%15.9f\n' % (self.grid_base[0,0],
                                                   self.grid_base[0,1],
                                                   self.grid_base[0,2])
            va = self.grid_base[1] - self.grid_base[0]
            header += '   %15.9f%15.9f%15.9f\n' % (va[0], va[1], va[2])
            vb = self.grid_base[2] - self.grid_base[0]
            header += '   %15.9f%15.9f%15.9f\n' % (vb[0], vb[1], vb[2])
            vc = self.grid_base[3] - self.grid_base[0]
            header += '   %15.9f%15.9f%15.9f' % (vc[0], vc[1], vc[2])
            ## write grid
            grid = self.grid_data.flatten(order='C') # [[nX] nY] nZ
            ## footer
            footer = ''
            left = grid.shape[0] % 5
            if left > 0:
                for i in range(grid.shape[0]-left, grid.shape[0]):
                    footer += '%15.6e' % grid[i]
                footer += '\n'
            footer += '   %s_%s\n' % ('END_DATAGRID_3D', grid_name)
            footer += ' %s\n' % 'END_BLOCK_DATAGRID_3D'
            grid = grid[:grid.shape[0]-left].reshape([-1, 5], order='C')
            np.savetxt(filename, grid, fmt='%15.6e', header=header, footer=footer, comments='')
            del grid
        # write 2D data
        else:
            header += ' %s\n' % 'BEGIN_BLOCK_DATAGRID_2D'
            header += '   %s\n' % grid_name
            header += '   %s_%s\n' % ('BEGIN_DATAGRID_2D', grid_name)
            header += '   %8i%8i\n' % (self.grid_data.shape[1], self.grid_data.shape[0]) # [nX] nY
            header += '   %15.9f%15.9f%15.9f\n' % (self.grid_base[0,0],
                                                   self.grid_base[0,1],
                                                   self.grid_base[0,2])
            va = self.grid_base[1] - self.grid_base[0]
            header += '   %15.9f%15.9f%15.9f\n' % (va[0], va[1], va[2])
            vb = self.grid_base[2] - self.grid_base[0]
            header += '   %15.9f%15.9f%15.9f\n' % (vb[0], vb[1], vb[2])
            ## write grid
            grid = self.grid_data.flatten(order='C') # [nX] nY
            ## footer
            footer = ''
            left = grid.shape[0] % 5
            if left > 0:
                last_line = ''
                for i in range(grid.shape[0]-left, grid.shape[0]):
                    footer += '%15.6e' % grid[i]
                footer += '\n'
            footer += '   %s_%s\n' % ('END_DATAGRID_2D', grid_name)
            footer += ' %s\n' % 'END_BLOCK_DATAGRID_2D'
            grid = grid[:grid.shape[0]-left].reshape([-1, 5], order='C')
            np.savetxt(filename, grid, fmt='%15.6e', header=header, footer=footer, comments='')
        return


class BXSF():
    """
    The collection of methods to parse and write XCrySDen BXSF files (Fermi
    surface).

    .. note::

        3D systems only.

    Args:
        relattice (array|CStructure): Reciprocal lattice matrix or structure
            object.
        bands (array): nBand\*nZ\*nY\*nX\*nSpin array. Unit: eV. Aligned to
            :math:`E_{F}=0`.
        efermi (float): Fermi energy. Unit: eV.
        band_index (list|str|int): Indices of bands. Starting from 1. For spin-
            polarized cases, one can specify '4a' for the :math:`\\alpha` state
            of band 4. If only band number is given for spin polarized bands,
            both spin states are saved. Use 'vb' and 'cb' for the highest
            valance or the lowest conduction band. ``[]`` for all bands. 
    """
    def __init__(self, rlattice, bands, efermi=0., band_index=[]):
        import warnings, re
        import numpy as np
        from pymatgen.core.structure import Structure
        from CRYSTALpytools.electronics import FermiSurface

        if isinstance(rlattice, Structure):
            self.rlattice = rlattice.lattice.reciprocal_lattice
        else:
            self.rlattice = np.array(rlattice, dtype=float)

        self.efermi = efermi

        bands = np.array(bands, dtype=float)
        pdirs = np.where(np.array(bands.shape[1:-1])>1)[0]
        if bands.ndim != 5:
            raise ValueError("Input band must be in the shape of nBand*nZ*nY*nX*nSpin.")
        if len(pdirs) != 3:
            raise Exception('BXSF format only applies to 3D structures. Read XCrySDen manual.')

        # Band index
        if len(band_index) == 0:
            band_index = [i+1 for i in range(bands.shape[0])]
        iband, ispin = FermiSurface._get_band_index(bands, band_index)

        # XCrySDen defines a general mesh
        self.bands = np.zeros([len(iband), bands.shape[1], bands.shape[2], bands.shape[3]])
        self.band_labels = np.zeros_like(np.zeros([len(iband),]), dtype=str)
        for i in range(len(iband)):
            self.bands[i, :, :, :] = bands[iband[i], :, :, :, ispin[i]]
            self.band_labels[i] = '{:d}0{:d}'.format(iband[i]+1, ispin[i]+1)
        self.bands = self.bands + self.efermi

    def write(self, filename, grid_name='UNKNOWN'):
        """
        Write Fermi surface data into a new XSF file. Band label is a string,
        with indices:

        * \[0\:-2\]: Band index strating from 1;  
        * -2: Useless, 0;  
        * -1: Spin. Alpha / No spin, 1; Beta, 2

        Args:
            filename (str): Output name.
            grid_name (str): Name of the grid.
        """
        import numpy as np
        import os, warnings

        if os.path.exists(filename):
            warnings.warn(
                "File '{}' exists! It will be overwritten.".format(filename),
                stacklevel=3)
        file = open(filename, 'w')

        # write 3D data
        file.write('# Generated by CRYSTALpytools\n')
        file.write(' BEGIN_INFO\n')
        file.write('   Fermi Energy: %12.6f\n' % self.efermi)
        file.write(' END_INFO\n\n')
        file.write(' BEGIN_BLOCK_BANDGRID_3D\n')
        file.write('   Fermi_surfaces_of_%s\n' % filename.replace(' ', '_'))
        file.write('   BEGIN_BANDGRID_3D_%s\n' % grid_name)
        file.write('     %i\n' % self.bands.shape[0])
        file.write('     %-6i%-6i%-6i\n' % (self.bands.shape[1],
                                            self.bands.shape[2],
                                            self.bands.shape[3]))
        file.write('     {:< 12.6f}{:< 12.6f}{:< 12.6f}\n'.format(0., 0., 0.))
        file.write('     {:< 12.6f}{:< 12.6f}{:< 12.6f}\n'.format(
            self.rlattice[0, 0], self.rlattice[0, 1], self.rlattice[0, 2]
        ))
        file.write('     {:< 12.6f}{:< 12.6f}{:< 12.6f}\n'.format(
            self.rlattice[1, 0], self.rlattice[1, 1], self.rlattice[1, 2]
        ))
        file.write('     {:< 12.6f}{:< 12.6f}{:< 12.6f}\n'.format(
            self.rlattice[2, 0], self.rlattice[2, 1], self.rlattice[2, 2]
        ))
        for label, band in zip(self.band_labels, self.bands):
            file.write('   BAND:%6s\n' % label)
            band = band.flatten(order='F')
            for i in range(band.shape[0] // 4):
                idx = int(i*4)
                file.write('       {:< 12.6f}{:< 12.6f}{:< 12.6f}{:< 12.6f}\n'.format(
                    band[idx], band[idx+1], band[idx+2], band[idx+3]
                ))
            if band.shape[0]%4 != 0:
                lastl = '       '
                for i in range(band.shape[0]%4, 0, -1):
                    lastl += '{:< 12.6f}'.format(band[-i])
                file.write('%s\n' % lastl)
        file.write('   END_BANDGRID_3D\n')
        file.write(' END_BLOCK_BANDGRID_3D\n\n')
        return

