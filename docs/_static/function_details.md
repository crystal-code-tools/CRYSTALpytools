# Details of the CRYSTALpytools

<h3>adsorb</h3>
<h6>sub_ads_indices(structure)</h6>

<b>Parameters</b>
<ul>
<li> <b>structure</b>: pymatgen Slab object. The function needs the 'adsorbate' or 'surface' site properties to be defined
</ul>
<b>Returns</b> dict (keys='adsorbate','substrate')(values=indices of the adsorbate/substrate atoms)

<h3>calculate</h3>

<h6>cry_ads_energy(e_full_system, e_substrate, e_adsorbate)</h6>

<b>Parameters</b>
<ul>
<li> <b>e_full_system</b>: energy of the surface+adsorbate
<li> <b>e_substrate</b>: energy of the substrate
<li> <b>e_adsorbate</b>: energy of the adsorbate
</ul>
<b>Returns</b> adsorption energy (float)

<h6>cry_shrink(structure, spacing=0.2)</h6>

<b>Parameters</b>
<ul>
<li> <b>structure</b>: pymatgen Structure object
<li> <b>spacing</b>: grid spacing in Å^(-1)
</ul>
<b>Returns</b>: the SHRINK value to be used in the CRYSTAL input.

<h6></h6>
<b>Parameters</b>
<ul>
<li> <b></b>:
</ul>
<b>Returns</b>

<h3>convert</h3>

<h6>cry_bands2pmg(output, bands, labels=None)</h6>
<b>Parameters</b>
<ul>
<li> <b>output</b>: CRYSTAL output object
<li> <b>bands</b>:CRYSTAL bands object
<li> <b>labels</b>: <b>k</b> point labels
</ul>
<b>Returns</b>

<h6>cry_gui2pmg(gui_file)</h6> Add lower dimensionality
<b>Parameters</b>
<ul>
<li> <b>gui_file</b>: CRYSTAL gui file
</ul>
<b>Returns</b>: pymatgen structure object

<h6>cry_out2pmg(output, initial=False, dimensionality=3, vacuum=10)</h6>
<b>Parameters</b>
<ul>
<li> <b>output</b>: CRYSTAL output file
<li> <b>initial</b>: read the first geometry in a geometry optimisation
<li> <b>dimensionality</b>: dimansionality of the system
<li> <b>vacuum</b>: vacuum to be used in the generation of the pymatgen object for slabs
</ul>
<b>Returns</b>: pymatgen Structure object

<h3>execute</h3>

<h6>runcry(file_name, guessp=None)</h6>
<b>Parameters</b>
<ul>
<li> <b>file_name</b>: name of the CRYSTAL input
<li> <b>guessp</b>: name of the .f20 file (if GUESSP used in the input)
</ul>
<b>Returns</b>: Confirmation of successful calculation

<h6>runprop(prop_name,wf_file)</h6>
<b>Parameters</b>
<ul>
<li> <b>prop_name</b>: name of the CRYSTAL properties input
<li> <b>wf_file</b>: name of the wavefunction file (.f9 or .f98)
</ul>
<b>Returns</b>: Confirmation of successful calculation

<h6>set_runcry_path(path)</h6>
<b>Parameters</b>
<ul>
<li> <b>path</b>: path to the runcry scrips (str)
</ul>
<b>Returns</b>: None

<h6>set_runprop_path(path)</h6>
<b>Parameters</b>
<ul>
<li> <b>path</b>: path to the runprop scrips (str)
</ul>
<b>Returns</b>: None

<h3>crystal_io</h3>

<h6>class Crystal_input(input_name)</h6>
<b>Parameters</b>
<ul>
<li> <b>input_name</b>: name of the input file (str)
</ul>
<b>Properties</b>:
<ul>
<li> is_basisset: whether the basis set is defined by using the BASISSET keyword (bool)
<li> geom_block: input geometry block (list)
<li> bs_block: input basis set block (list)
<li> func_block: input functional block (list)
<li> scf_block: input scf parameters block (list)
</ul>

<h6>add_ghost(ghost_atoms)</h6>
<b>Parameters</b>
<ul>
<li> <b>ghost_atoms</b>: indices of the atoms to transform into ghosts
</ul>

<h6>opt_to_sp</h6>
Removes the OPTGEOM part of the input.

<h6>opt_to_sp</h6>
Adds the OPTGEOM part to the input.


<h6>class Crystal_output(output_name)</h6>
<b>Parameters</b>
<ul>
<li> <b>output_name</b>: name of the output file (str)
</ul>
<b>Properties</b>:
<ul>
<li> name: output name (str)
<li> converged: whether the calculation ended correctly (bool)
</ul>

<h6>final_energy()</h6>
<b>Properties</b>
<ul>
<li> <b>energy</b>: final energy in eV (float)
</ul>

<h6>scf_convergence(all_cycles=False)</h6>
<b>Parameters</b>
<ul>
<li> <b>all_cycles</b>: whether to return the scf convergence energy for all cycles. This is relevant in a geometry optimisation calculation. If set to False only the last scf will be analysed.
</ul>
<b>Properties</b>
<ul>
<li> <b>scf_energy</b>: energy of each scf cycle (list)
<li> <b>scf_deltae</b>: energy difference wrt the previous cycle
</ul>

<h6>num_cycles()</h6>
<b>Properties</b>
<ul>
<li> <b>num_cycles_scf</b>: number of cycles of the last scf.
</ul>


<h6>fermi_energy()</h6>
<b>Properties</b>
<ul>
<li> <b>efermi</b>: Fermi energy in eV (float)
<li> <b>efermi_from</b>: which value was used for the Fermi energy (possible values: 'scf' or 'scf_top_valence')
</ul>


<h6>primitive_lattice(initial=True)</h6>
<b>Parameters</b>
<ul>
<li> <b>initial</b>: initial (True) or optimised primitive lattice (False)
</ul>
<b>Properties</b>
<ul>
<li> <b>primitive_vectors</b>: primitive cell lattice vectors in Å (3x3 numpy array)
</ul>

<h6>reciprocal_lattice(initial=True)</h6>
<b>Parameters</b>
<ul>
<li> <b>initial</b>: initial (True) or optimised reciprocal lattice (False)
</ul>
<b>Properties</b>
<ul>
<li> <b>reciprocal_vectors</b>: primitive cell reciprocal vectors in Å^(-1)(3x3 numpy array)
</ul>

<h6>get_band_gap</h6>
<b>Properties</b>
<ul>
<li> <b>spin_pol</b>: set to True if the calculation is spin polarised. False otherwise.
<li> <b>band_gap</b>: band gap in eV (float)
</ul>

<h6>extract_last_geom(write_gui_file=True,print_cart=False)</h6>
<b>Parameters</b>
<ul>
<li> <b>write_gui_file</b>: write the geometry into a .gui file (same name as the output file)
<li> <b>print_cart</b>: print the geometry to output
</ul>
<b>Properties</b>
<ul>
<li> <b>opt_converged</b>: whether the geometry optimisation has converged (True) or not (False)
<li> <b>trans_matrix</b>: primitive to conventional transformation matrix (3x3 numpy array)
<li> <b>n_atoms</b>: number of atoms
<li> <b>atom_symbols</b>: atom symbols (list)
<li> <b>atom_numbers</b>: atoms numbers (list)
<li> <b>atom_positions_cart</b>: cartesian coordinates of the atoms (numpy array)
<li> <b>cart_coords</b>: numpy array containing for each line the atomic number and its coordinates.
</ul>


<h6>symm_ops()</h6>
<b>Properties</b>
<ul>
<li> <b>symmops</b>: symmetry operators (N_atoms x 12 array)
</ul>


<h6>forces(initial=False,grad=False)</h6>
<b>Parameters</b>
<ul>
<li> <b>initial</b>: whether to read the last (False) or first (True) geometry optimisation
<li> <b><grad/b>: whether to read the gradient (True) or not (False)
</ul>
<b>Properties</b>
<li> <b>forces_atoms</b>: forces on the atoms
<li> <b>forces_cell</b>: forces on the cell
<li> <b>grad</b>: gradient
<li> <b>rms_grad</b>: rms on the gradient
<li> <b>disp</b>: displacement
<li> <b>rms_disp</b>: rms on the displacement


<h6>config_analysis</h6>
To use after a SCELCONF calculation.
<b>Properties</b>
<ul>
<li> <b>atom_type1</b>: indices of the first atomic species (list of lists)
<li> <b>atom_type2</b>: indices of the second atomic species (list of lists)
</ul>

<h6></h6>
<b>Properties</b>
<ul>
<li> <b></b>:
</ul>

<h3>plot</h3>
<h6>plot_cry_bands(bands, k_labels=None, energy_range=None, title=False, not_scaled=False, mode='single', linestl='-', linewidth=1, color='blue', fermi='forestgreen', k_range=None, labels=None, figsize=None, scheme=None, sharex=True, sharey=True)<\h6>
<ul>
<li>
<li>
<li>
<li>
<li>
</ul>
<h3>properties</h3>
<ul>
<li>
<li>
<li>
<li>
<li>
</ul>
