# Efficient Bayesian phase estimation using mixed priors

This code is provided to allow reproducibility of the results in the paper:
E. van den Berg, "Efficient Bayesian phase estimation using mixed priors"
[arXiv:2007.11629](https://arxiv.org/abs/2007.11629).

The 'cache' directory contains all pre-computed results. These can be
plotted using the relevant Python scripts:

*  generate_figure_2.py
*  generate_figure_3ab.py
*  generate_figure_3c.py
*  generate_figure_3de.py
*  generate_figure_3f.py
*  generate_figure_4ab.py
*  generate_figure_4c.py
*  generate_figure_4d.py
*  generate_figure_4e.py
*  generate_figure_4f.py
*  generate_figure_5a.py
*  generate_figure_5b.py
*  generate_figure_5d.py
*  generate_figure_5e.py
*  generate_figure_6.py

The tables in the paper can be generated using:

*  generate_tables_1_2.py

The output of the generate figure scripts is one or more figures in PDF format,
stored in the 'fig' directory. Automatic cropping of the figures using 'pdfcrop'
has been disabled because it is highly system dependent. However, it can be
enabled by uncommenting the relevant section of the 'exportFigure' function in
the 'generic.py' file.

The generate-figure scripts load precomputed results from the 'cache' directory
and generate the plots. To reproduce the results themselves, delete the relevant
experiment files 'cache/experiment_hybrid_<index>.dat' (likewise for the files
'cache/experiment_transition_<index>.dat') and regenerate the data as follows:

1. Run `make run_experiment_hybrid` or manually compile the files
   (run_experiment_hybrid.c, bpe.c, and bpe_plugins.c)

2. Create the 'cache' directory if it does not already exist

3. On Mac and Linux run the experiment and redirect the output to file:
     `./run_experiment_hybrid <index>  > ./cache/experiment_hybrid_<index>.dat`
   where `<index>` is the experiment index needed by the figure generation script.

4. Once all data files have been generated, run the figure generation script.

For missing `cache/experiment_transition_<index>.dat` files, compile the
executable using `make run_experiment_transition`, and run it with
output redirected to `cache/experiment_transition_<index>.dat`. This part of
the code, as well as the generation of the test problems could be improved.
We decided to leave the code as-is to avoid introducing errors or inconsistencies
with the data that was already generated.

Some verification of the code or equations is provided by
the verify* scripts. The `verify_clib_fourier.py` script first calls make
to compile `verify_clib_fourier.c` and then runs this code for comparison
with the python implementation of the Fourier representation in `bpe.py`.

The `compress.py` script truncates data files for problems in the 200 and 300
ranges, since only data from the last several iterations is actually used
(aside from two specific problem instances used for plotting). We added this
script to reduce the file size for these problems and make it easier to check
out the repository. When locally regenerating the data it can optionally be
truncated by running the `compress.py` script.

