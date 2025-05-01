Example 02: FreeLigand
======================

This example demonstrates how to use the FreeLigand class to perform a multi-state RESP calculation.

The key difference between this and :doc:`01_LazyLigand` is that this
example uses the gaussian rotation calculations to calculate many resp calcualtions over a wide-range of rotations. These types of rotations
ared useful for ensuring that the choice of grid points in the RESP calculation is not biased by the initial orientation of the molecule.

Learning Outcomes:
------------------

1) Learn how to use the FreeLigand class to perform a gaussian rotational RESP calculation to parametrize a ligand using multi-state resp fitting.
2) Learn how to list and execute stages in the pipeline.

Files 
-----
The files for this example can be found in the `LigandParameterization/examples/02_FreeLigand` directory of the source code.


Tutorial 
--------

To start with any of the premade recipes, you need to import the recipe from the module. Using the wild-card character will import
all the recipes available in the module, whereas, importing a specific recipe will only import that recipe. Here, we import the FreeLigand recipe. 


.. literalinclude :: ../../../../examples/02_FreeLigand/param_thio.py
    :language: python
    :end-at: gaussian_scratch


The next step is to load the pdb file into the recipe of your choosing, and set various machine parameters. Full parameters that can be selected
are available in the documentation for the class (for instance, for FreeLigand, see :class:`ligandparam.module.FreeLigand`).

.. literalinclude :: ../../../../examples/02_FreeLigand/param_thio.py
    :language: python
    :start-at: parametrize_ligand =
    :end-before: Select the pre-initialized stages


Note that here, we have chosen the pdb thiophenol.pdb and set the net charge to 0 using the *net_charge* parameter. If you have a charged ligand, you need to select the
proper net charge, as later when Gaussian runs it will need to know the net charge of the molecule.  Note that this class is called nearly identically to the LazyLigand class, 
but with the FreeLigand class used instead.

The next step is to select the pre-initialized stages for the FreeLigand class. This can be done using the *setup* method. The *nproc* and *mem* parameters are used to set the 
number of processors and the memory to be used by Gaussian, respectively.


.. literalinclude :: ../../../../examples/02_FreeLigand/param_thio.py
    :language: python
    :start-at: Select the pre-initialized stages
    :end-at: dry_run=False 

The FreeLigand class has a number of stages that are pre-initialized. 

In brief, these are:

1) Convert the pdb file to a mol2 and assign initial charges
2) Run Gaussian to assign resp charges
3) Do a Gaussian rotation calculation where the ligand is rotated in 3D space and RESP calculations are repeated.
4) Write out lib/frcmod files for the ligand


The output files will be generated in the same directory as the input pdb file, and will have the same name as the pdb file, but with different extensions.

These files are:

- thiophenol.resp.mol2 - The final mol2 file with the RESP charges.

- thiophenol.frcmod - The frcmod file for the ligand.

- thiophenol.off - The off(lib) parameter file for the ligand.

Note that the charges ion the mol2 file should be similar, but not exactly the same as the charges that you obtained from the LazyLigand class.

Full code
---------
.. literalinclude :: ../../../../examples/02_FreeLigand/param_thio.py
    :language: python
