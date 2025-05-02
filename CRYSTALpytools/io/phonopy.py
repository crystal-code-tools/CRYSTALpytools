#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods to parse files used by `Phonopy <https://phonopy.github.io/phonopy/>`_.
Currently only YAML and FORCE_CONSTANTS files are supported.
"""
import numpy as np
from yaml import safe_load
from warnings import warn

try:
    from phonopy.structure.atoms import PhonopyAtoms
    from phonopy import Phonopy
    from phonopy import units as punits
    from spglib import find_primitive
except ModuleNotFoundError:
    raise ModuleNotFoundError("Phonopy is required for this module.")

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

from CRYSTALpytools import units
from CRYSTALpytools.geometry import CStructure



def _phonon_rep(dumper, value):
    """Formatting output phonon info."""
    return dumper.represent_scalar('tag:yaml.org,2002:float', '{0:.10f}'.format(value))
def _header_rep(dumper, value):
    """Formatting output header info."""
    return dumper.represent_scalar('tag:yaml.org,2002:float', '{0:.15f}'.format(value))


class YAML():
    """
    Read and write phonopy YAML files. `Phonopy python API <https://phonopy.github.io/phonopy/phonopy-module.html>`_
    is used for writing. The ``Phonopy`` instance is saved in ``self._phonopy``
    attribute, with basic geometry and calculator information.

    .. note::

        Phonopy default units (:math:`\\AA`, AMU, THz, eV) are used for all
        calculators, which leads to no error as far as developers are aware of.

    Args:
        struc (Structure|CStructure): Pymatgen Structure class, unit cell.
        dim (list[int]): 1\*3 or 1\*9 list of supercell expansion matrix, i.e.,
            the ``--dim`` option of phonopy.
        calculator (str): Name of calculator. Will be used to determine the
            conversion factors.
        primitive (str|array): 9\*1 primitive matrix in phonopy convention, or
            'auto' to automatically identify the primitive cell.
        \*\*kwargs: Other attributes. Listed below.
        qpoint (array): nQpoint\*4, Fractional coordinates and weight.
        frequency (array): nQpoint\*nMode, In THz.
        mode_symm (array): nQpoint\*nMode str, in Mulliken symbols.
        eigenvector (array): nQpoint\*nMode\*nAtom\*3 complex, Mass-weighted and phased, normalized to 1.
    """
    def __init__(self, struc, dim, calculator='crystal', primitive='auto', **kwargs):
        if not isinstance(struc, Structure):
            raise TypeError("A pymatgen Structure or CRYSTALpytools CStructure object must be used.")

        freqfactor = {
            'vasp'     : 'punits.VaspToTHz',
            'wien2k'   : 'punits.Wien2kToTHz',
            'qe'       : 'punits.PwscfToTHz',
            'abinit'   : 'punits.AbinitToTHz',
            'siesta'   : 'punits.SiestaToTHz',
            'elk'      : 'punits.ElkToTHz',
            'crystal'  : 'punits.CrystalToTHz',
            'turbomole': 'punits.TurbomoleToTHz',
            'cp2k'     : 'punits.CP2KToTHz',
            'fhi-aims' : 'punits.VaspToTHz',
            'fleur'    : 'punits.FleurToTHz',
            'castep'   : 'punits.CastepToTHz',
            'abacus'   : 'punits.AbinitToTHz',
            'lammps'   : 'punits.VaspToTHz'
        }
        if calculator.lower() not in freqfactor.keys():
            raise ValueError("Calculator not supported: '{}'.".format(calculator))
        freqfac = eval(freqfactor[calculator.lower()])

        # supercell
        dim = np.array(dim, dtype=int, ndmin=1)
        s_matrix = np.zeros([3,3], dtype=int) # Phonopy convention
        if dim.shape[0] == 3:
            for i in range(3): s_matrix[i,i] = dim[i]
        elif dim.shape[0] == 9:
            for i in range(9): s_matrix[i//3, i%3] = dim[i]
        elif dim.shape[0] == 1:
            for i in range(3): s_matrix[i,i] = dim[0]
        else:
            raise ValueError('Dimensionality must have the length of 1, 3, or 9.')

        # primitive cell
        if isinstance(primitive, str) and primitive.lower() == 'auto':
            cell = (struc.lattice.matrix, struc.frac_coords, [i.Z for i in struc.species])
            p_cell, p_coords, p_species = find_primitive(cell)
            p_matrix = np.linalg.inv(struc.lattice.matrix.T) @ p_cell.T # Phonopy convention
        else:
            p_matrix = np.array(primitive, ndmin=1, dtype=float)
            if p_matrix.shape[0] != 9: raise ValueError('Primitive axes must be a 1*9 1D list.')
            p_matrix = p_matrix.reshape([3,3])

        atom = PhonopyAtoms(symbols=[i.symbol for i in struc.species],
                            cell=struc.lattice.matrix,
                            scaled_positions=struc.frac_coords)
        self._phonopy = Phonopy(atom,
                                supercell_matrix=s_matrix,
                                primitive_matrix=p_matrix,
                                factor=freqfac,
                                calculator=calculator.lower())
        self._phonopy._build_primitive_cell()
        self.natom = self._phonopy.primitive.scaled_positions.shape[0]
        self.structure = CStructure(self._phonopy.primitive.cell,
                                    self._phonopy.primitive.symbols,
                                    self._phonopy.primitive.scaled_positions)

        # Other attributes
        for key, value in zip(kwargs.keys(), kwargs.values()):
            setattr(self, key, value)
            if key == 'qpoint': self.nqpoint = len(self.qpoint)
            elif key == 'frequency': self.nmode = len(self.frequency[0])

    @classmethod
    def read(cls, struc, phonon=''):
        """
        Read data from YAML. Currently 'phonopy', 'phononpy_disp' (structure only),
        'mesh', 'band', 'qpoints' and 'irreps' are supported.

        Args:
            struc (str): Geometry information. Only for 'phonopy', 'phononpy_disp', 'mesh' and 'band'
            phonon (str): Frequency information, including q points, frequency,
                eigenvector and irreducible representations. For 'mesh', 'band',
                'qpoints' and 'irreps'.
        Returns:
            cls
        """
        struc_file= open(struc, 'r')
        data = safe_load(struc_file)
        struc_file.close()

        # Structure
        ## file type, phonopy and phonopy-disp are not distinguished
        ftype = ''
        for key, file in zip(['space_group', 'segment_nqpoint', 'mesh'],
                             ['phonopy', 'band', 'mesh']):
            try:
                _ = data[key]
                ftype = file
                break
            except KeyError:
                continue
        if ftype == '': raise Exception("Unknown file format for structure. Only 'phonopy', 'phonopy_disp', 'mesh' and 'band' are read.")

        try:
            ## unit
            if ftype == 'phonopy':
                len_unit = data['physical_unit']['length']
            elif ftype == 'band':
                len_unit = data['length_unit']
            elif ftype == 'mesh':
                warn("Unknown length unit. 'angstrom' is assumed.", stacklevel=2)
                len_unit = 'angstrom'

            if len_unit == 'angstrom':
                unit_len = 1.0
            elif len_unit == 'au':
                unit_len = units.au_to_angstrom(1.0)
            else:
                raise Exception("Unknown length unit. Available options: au, angstrom.")
            ## geometry
            spec = []; coord = []
            if ftype == 'band' or ftype == 'mesh':
                latt = np.array(data['lattice'], dtype=float) * unit_len
                for idx_a, atom in enumerate(data['points']):
                    spec.append(atom['symbol'])
                    coord.append(atom['coordinates'])
                calculator = 'crystal'
                smatrix = [1,1,1]
                pmatrix = np.eye(3).flatten()
            elif ftype == 'phonopy':
                latt = np.array(data['primitive_cell']['lattice'], dtype=float) * unit_len
                for idx_a, atom in enumerate(data['primitive_cell']['points']):
                        spec.append(atom['symbol'])
                        coord.append(atom['coordinates'])
                calculator = data['phonopy']['calculator']
                try:
                    smatrix = np.array(data['supercell_matrix'], dtype=int).flatten()
                except KeyError:
                    smatrix = [1,1,1]
                try:
                    pmatrix = np.array(data['primitive_matrix'], dtype=float).flatten()
                except KeyError:
                    pmatrix = np.eye(3).flatten()

            struc = CStructure(latt, spec, coord, )
        except KeyError:
            raise Exception("Geometry file: '{}' is broken. Check your input file.".format(struc))

        if phonon == '' and ftype == 'phonopy':
            return cls(struc, smatrix, calculator=calculator, primitive=pmatrix)

        # Phonon
        if phonon != '':
            phonon_file= open(phonon, 'r')
            data = safe_load(phonon_file)
            phonon_file.close()

        ## file type, qpoints must put at the end
        ftype = ''
        for key, file in zip(['segment_nqpoint', 'mesh', 'point_group', 'phonon'],
                             ['band', 'mesh', 'irreps', 'qpoints']):
            try:
                _ = data[key]
                ftype = file
                break
            except KeyError:
                continue
        if ftype == '': raise Exception("Unknown file format for phonons. Only 'band', 'mesh', 'irreps' and 'qpoints' are read.")

        try:
            ## phonon
            if ftype == 'band' or ftype == 'mesh' or ftype == 'qpoints':
                nqpoint = data['nqpoint']
                nmode = int(struc.num_sites*3)
                natom = struc.num_sites
                try:
                    _ = data['phonon'][0]['band'][0]['eigenvector']
                    read_eigvec = True
                except KeyError:
                    read_eigvec = False

                try:
                    _ = data['phonon'][0]['weight']
                    read_weight = True
                except KeyError:
                    read_weight = False

                qpoint = np.zeros([nqpoint, 4], dtype=float)
                frequency = np.zeros([nqpoint, nmode], dtype=float)
                mode_symm = np.array([['' for i in range(nmode)] for j in range(nqpoint)], dtype=str)
                if read_eigvec == True:
                    eigenvector = np.zeros([nqpoint, nmode, natom, 3], dtype=complex)
                    for iq in range(nqpoint):
                        qpoint[iq, 0:3] = data['phonon'][iq]['q-position']
                        if read_weight == True:
                            qpoint[iq, 3] = data['phonon'][iq]['weight']
                        else:
                            qpoint[iq, 3] = 1
                        for im in range(nmode):
                            frequency[iq, im] = data['phonon'][iq]['band'][im]['frequency']
                            eigvt = np.array(data['phonon'][iq]['band'][im]['eigenvector'], dtype=float)
                            eigenvector[iq, im] = eigvt[:, :, 0] + eigvt[:, :, 1]*1j
                else:
                    eigenvector = np.array([[[[] for i in range(natom)] for j in range(nmode)] for k in range(nqpoint)])
                    for iq in range(nqpoint):
                        qpoint[iq, 0:3] = data['phonon'][iq]['q-position']
                        if read_weight == True:
                            qpoint[iq, 3] = data['phonon'][iq]['weight']
                        else:
                            qpoint[iq, 3] = 1
                        for im in range(nmode):
                            frequency[iq, im] = data['phonon'][iq]['band'][im]['frequency']
            elif ftype == 'irreps':
                nqpoint = 1
                nmode = int(struc.num_sites*3)
                natom = struc.num_sites

                qpoint = np.zeros([1, 4], dtype=float) + 1
                qpoint[0, 0:3] = data['q-position']
                frequency = np.zeros([nqpoint, nmode], dtype=float)
                mode_symm = []
                eigenvector = np.array([[[[] for i in range(natom)] for j in range(nmode)] for k in range(nqpoint)])

                im = 0
                for normode in data['normal_modes']:
                    neq = len(normode['band_indices'])
                    freq = normode['frequency']
                    symm = str(normode['ir_label'])
                    if symm == 'None': symm = ''
                    for i in range(neq):
                        frequency[0, im+i] = freq
                        mode_symm.append(symm)
                    im += neq
                mode_symm = np.array(mode_symm, dtype=str).reshape([nqpoint, nmode])

        except KeyError:
            raise Exception("Phonon file: '{}' is broken. Check your input file.".format(phonon))

        return cls(struc, smatrix, calculator=calculator, primitive=pmatrix,
                   qpoint=qpoint, frequency=frequency, mode_symm=mode_symm, eigenvector=eigenvector)

    def write_phonopy(self, filename='phonopy.yaml'):
        """
        Save computational setups and structure into 'phonopy.yaml'.

        Args:
            filename (str): The YAML file name.
        Returns:
            None
        """
        self._phonopy.save(filename=filename)
        return

    def write_qpoints(self, filename='qpoints.yaml', write_eigenvector=False):
        """
        Write vibration data into 'qpoints.yaml'.

        Args:
            filename (str): 'qpoints' formatted file name.
            write_eigenvector (bool): Whether to write eigenvector if present.
        Returns:
            None
        """
        qcoords = self.qpoint[:, 0:3]

        lattmx = self._phonopy.primitive.cell
        rlattmx = Lattice(lattmx).reciprocal_lattice_crystallographic.matrix

        # header
        header = dict(nqpoint=int(self.nqpoint),
                      natom=int(self.natom),
                      reciprocal_lattice=rlattmx.tolist())

        file = open(filename, 'w')
        yaml.add_representer(float, _header_rep)
        yaml.dump(header, file, sort_keys=False, default_flow_style=None)
        file.write('\n')
        del header

        # phonons
        phonon = self._write_phonon(rlattmx, 'qpoints', write_eigenvector)

        yaml.add_representer(float, _phonon_rep)
        yaml.dump({'phonon' : phonon}, file, sort_keys=False, default_flow_style=None)
        file.write('\n')
        file.close()
        return

    def write_mesh(self, filename='mesh.yaml', write_eigenvector=False):
        """
        Write vibration data into 'mesh.yaml'. The mesh size is inferred from
        coordinates of qpoints, so it is important to use data obtained from
        the regular mesh grid.

        Args:
            filename (str): 'mesh' formatted file name.
            write_eigenvector (bool): Whether to write eigenvector if present.
        Returns:
            None
        """
        qcoords = self.qpoint[:, 0:3]
        qweight = self.qpoint[:, 3].flatten()
        # infer mesh size
        idx = np.where(qcoords!=0)
        if len(idx) == 0:
            mesh = [1, 1, 1]
        else:
            inrc = np.min(np.abs(qcoords[idx[0], idx[1], idx[2]]), axis=0)
            mesh = np.array(np.round(0.5/inrc), dtype=int)

        lattmx = self._phonopy.primitive.cell
        rlattmx = Lattice(lattmx).reciprocal_lattice_crystallographic.matrix

        # header
        header = dict(mesh=mesh.tolist(),
                      nqpoint=int(self.nqpoint),
                      reciprocal_lattice=rlattmx.tolist(),
                      natom=int(self.natom),
                      lattice=lattmx.tolist())
        points = []
        for crd, ele, mas in zip(self._phonopy.primitive.scaled_positions,
                                 self._phonopy.primitive.symbols,
                                 self._phonopy.primitive.masses):
            points.append(dict(symbol=str(ele), coordinates=crd.tolist(), mass=float(mas)))

        file = open(filename, 'w')
        yaml.add_representer(float, _header_rep)
        yaml.dump(header, file, sort_keys=False, default_flow_style=None)
        yaml.dump({'points' : points}, file, sort_keys=False, default_flow_style=None)
        file.write('\n')
        del points, header

        # phonons
        phonon = self._write_phonon(rlattmx, 'mesh', write_eigenvector)

        yaml.add_representer(float, _phonon_rep)
        yaml.dump({'phonon' : phonon}, file, sort_keys=False, default_flow_style=None)
        file.write('\n')
        file.close()
        return

    def write_band(self, filename='band.yaml'):
        """
        Write phonon band structure into 'band.yaml'. The line segment is
        inferred from coordinates of qpoints, so it is important to use data
        sampled by a path.

        Args:
            filename (str): 'band' formatted file name.
        Returns:
            None
        """
        qcoords = self.qpoint[:, 0:3]
        if self.nqpoint <= 3:
            raise Exception("Only {:d} q points are present in this object. Are you sure it is a phonon band output?".format(self.nqpoint))

        lattmx = self._phonopy.primitive.cell
        rlattmx = Lattice(lattmx).reciprocal_lattice_crystallographic.matrix

        # infer path
        path = []; a_path = [0, 1]; org = qcoords[0]
        for i in range(1, self.nqpoint-1): # only save i+1
            vec0 = np.round(qcoords[i] - org, 8)
            vec1 = np.round(qcoords[i+1] - org, 8)

            dvec = np.linalg.norm(vec0-vec1)
            dcos = np.round(1-np.dot(vec0, vec1)/np.linalg.norm(vec0)/np.linalg.norm(vec1), 8)
            if dvec < 1e-8 or dcos < 1e-8: # overlapped q point or different direction
                path.append(a_path)
                a_path = [i+1]
                org = qcoords[i+1]
            else:
                a_path.append(i+1)
        path.append(a_path)

        nqpoint = 0; npath = 0.; pathpt = []
        for p in path:
            nqpoint += len(p)
            npath += 1
            pathpt.append(len(p))

        # header
        header = dict(calculator=self._phonopy.calculator,
                      length_unit='angstrom',
                      nqpoint=int(nqpoint),
                      npath=int(npath),
                      segment_nqpoint=pathpt,
                      reciprocal_lattice=rlattmx.tolist(),
                      natom=int(self.natom),
                      lattice=lattmx.tolist())
        points = []
        for crd, ele, mas in zip(self._phonopy.primitive.scaled_positions,
                                 self._phonopy.primitive.symbols,
                                 self._phonopy.primitive.masses):
            points.append(dict(symbol=str(ele), coordinates=crd.tolist(), mass=float(mas)))

        file = open(filename, 'w')
        yaml.add_representer(float, _header_rep)
        yaml.dump(header, file, sort_keys=False, default_flow_style=None)
        yaml.dump({'points' : points}, file, sort_keys=False, default_flow_style=None)
        file.write('\n')
        del header, points

        # phonons
        phonon = self._write_phonon(rlattmx, 'band', False)

        yaml.add_representer(float, _phonon_rep)
        yaml.dump({'phonon' : phonon}, file, sort_keys=False, default_flow_style=None)
        file.write('\n')
        file.close()
        return

    def _write_phonon(self, rlattmx, method, eigvec):
        """Internal method for dumping 'mesh', 'qpoints' and 'band' files.

        Args:
            rlattmx (array): Reciprocal lattice matrix.
            method (str): 'mesh', 'qpoints' or 'band'.
            eigvec (bool): Whether to dump eigenvectors.
        Returns:
            phonon (list[dict]): Phonon information.
        """
        if not hasattr(self, 'qpoint') or not hasattr(self, 'frequency'):
            raise Exception("The object must have phonon information.")
        if not hasattr(self, 'eigenvector'):
            if eigvec == True:
                warn('The object does not have eigenvector attribute.')
            eigvec = False

        method = method.lower()
        if method not in ['qpoints', 'band', 'mesh']:
            raise Exception("Unknown method: '{}'.".format(method))

        # Frequency
        phonon = []
        if method == 'band':
            count_dist = 0.
        if write_eigenvector == False:
            for crd, wei, freq in zip(qcoords, qweight, self.freqency):
                dic = {'q-position' : crd.tolist(),
                       'band' : [dict(frequency=float(i)) for i in freq]}

                if method == 'mesh':
                    dic['distance_from_gamma'] = float(np.linalg.norm(crd@rlattmx))
                    dic['weight'] = int(wei)
                elif method == 'qpoints':
                    pass
                elif method == 'band':
                    count_dist += np.linalg.norm((crd-self.qpoint[0, 0:3]) @ rlattmx)
                    dic['distance'] = float(count_dist)
                phonon.append(dic)
        else:
            eigvt = np.stack([self.eigenvector.real, self.eigenvector.imag], axis=4)
            for crd, wei, freq, eigv in zip(qcoords, qweight, self.freqency, eigvt):
                dic = {'q-position' : crd.tolist(),
                       'band' : []}
                bd = []
                for im in range(self.nmode):
                    bd.append({'frequency' : float(freq[iq, im]),
                               'eigenvector' : eigv[im, :, :, :].tolist()})
                dic['band'] = bd

                if method == 'mesh':
                    dic['distance_from_gamma'] = float(np.linalg.norm(crd@rlattmx))
                    dic['weight'] = int(wei)
                elif method == 'qpoints':
                    pass
                elif method == 'band':
                    count_dist += np.linalg.norm((crd-self.qpoint[0, 0:3]) @ rlattmx)
                    dic['distance'] = float(count_dist)
                phonon.append(dic)

        return phonon


def read_FC(input='FORCE_CONSTANTS'):
    """Read force constant (Hessian) matrix in Phonopy/VASP FORCE_CONSTANTS format.
    Units: 'eV' and ':math:`\\AA`'.

    Args:
        input (str): The input file name
    Returns:
        hess (array): nMode\*nMode array. Mass-unweighted Hessian matrix.
    """
    from pandas import DataFrame

    file = open(input, 'r')
    df = DataFrame(file)
    file.close()

    na = df[0].loc[0].strip().split()
    if len(na) == 1:
        na1 = int(na[0]); na2 = int(na[0])
    else:
        na1 = int(na[0]); na2 = int(na[1])

    headers = df[df[0].str.contains(r'^\s*[0-9]+\s+[0-9]+\s*$')].index.tolist()[1:]
    nmode = np.max([na1, na2])*3
    hess = np.zeros([nmode, nmode])
    for ih in range(len(headers)):
        mx = np.array(
            df[0].loc[headers[ih]+1:headers[ih]+3].map(lambda x: x.strip().split()).tolist(),
            dtype=float
        )
        at1, at2 = df[0].loc[headers[ih]].strip().split()
        at1 = int(at1); at2 = int(at2)
        at1 = int(3*at1-3); at2 = int(3*at2-3)
        hess[at1:at1+3, at2:at2+3] = mx

    if na1 < na2:
        repeat = na2 // na1
        for i in range(2, repeat):
            hess[int((na1-1)*i*3):int(na1*i*3), :] = hess[0:int(na1*3), :]
    elif na1 > na2:
        repeat = na1 // na2
        for i in range(2, repeat):
            hess[:, int((na2-1)*i*3):int(na2*i*3)] = hess[:, 0:int(na2*3)]
    return hess


def write_FC(hess, output='FORCE_CONSTANTS'):
    """Write force constant (Hessian) matrix into Phonopy/VASP FORCE_CONSTANTS format.
    Input units: 'eV' and ':math:`\\AA`'.

    Args:
        hess (array): nMode\*nMode array. Mass-unweighted Hessian matrix.
        output (str): The output name
    Returns:
        None
    """
    hess = np.array(hess, dtype=float, ndmin=2)
    natom = int(hess.shape[0] / 3)

    file = open(output, 'w')
    file.write('%4i%4i\n' % (natom, natom))
    for i in range(natom):
        for j in range(natom):
            file.write('%4i%4i\n' % (i + 1, j + 1))
            submx = hess[int(3*i):int(3*i+3), int(3*j):int(3*j+3)]
            for d in submx:
                file.write('%22.15f%22.15f%22.15f\n' % (d[0], d[1], d[2]))
    file.close()
    return






