# Structure of the codes

* **Generic classes.** These are self-contained classes or sets of functions used throughout the code base.

    * **filters.hpp**: contains the ```filter_database``` class that allows you to load broad-band filters from ASCII or FITS files and automatically re-normalize them and (if asked) convert them to the right filter type (integrate photons or energy, and flux density per unit wavelength or frequency).

    * **metrics.hpp**: contains the ```metrics``` class that holds a set of PSF metrics including the quadrupole moments, the size-squared, and the ellipticities. This class has well-defined mathematical operators that act on all metrics at once, and secondary metrics (size-squared and ellipticities) can be re-computed from the quadrople moments on request.

    * **psf_moments.hpp**: contains the ```psf_moments``` class that holds an array of monochromatric PSF quadrupole moments, and functions to integrate them into a broadband PSF in a way that reproduces perfectly the output of the Euclid PSF toolkit.

    * **rebin.hpp**: contains functions to re-bin an SED into a new grid using various methods.

    * **tree.hpp**: contains a tree-like data structure that has fast lookup and insertions but rather slow sequential traversal.


* **Generator classes.** These are classes used to *generate* data to be fed to a *consumer* class (fitter, averager, etc.).

    * **egg-analytic.hpp**: generate noise-less photometry of galaxies using the EGG analytical generator.

    * **psf-averager-photoz.hpp**: extends ```egg-analytic``` to generate noisy photometry, compute true PSF moments, interface with photo-z codes, and provide a fine control of what data gets stored on disk (everything, or averages).

    * **psf-averager-stars.hpp**: generates photometry of stars using a simple model, generates noisy photometry, computes true PSF moments, and interfaces for SED fitting codes.

    * **fitter-wrapper.hpp**: generate photometry by reading it from a catalog on disk.

    * **generator-single.hpp** [deprecated]: contains a generator with a 1D familly of SEDs and no redshift.


* **Fitter classes.** These classes use photometry to fit a model and produce outputs, including PSF moments.

    * **fitter-base.hpp**: root class for all fitters which defines the interface, also contains a "null" fitter (which does nothing, and is just a placeholder for interfaces that require one).
    * **fitter-eazy.hpp**: implements the EAzY galaxy model.
    * **fitter-bpz.hpp**: implements the BPZ galaxy model.
    * **fitter-star.hpp**: implements a simple star fitter.


* **Generate libraries of SEDs.**

    * **make_color_cube.cpp**: uses a catalog of photometry (with or without noise) and spectro-photometry (medium bands) to build an empirical, observer-frame galaxy SED library. The galaxies in the input catalog are gridded in color space and their SEDs are averaged together.


* **Fit fluxes in a catalog.**

    * **eazy.cpp**: uses ```fitter-eazy``` and ```fitter-wrapper``` to run the EAzY model on fluxes provided in a catalog.

    * **bpz.cpp**: uses ```fitter-bpz``` and ```fitter-wrapper``` to run the BPZ model on fluxes provided in a catalog.

    * **get_psf_color_cube.cpp**: use a galaxy observer-frame SED library (no redshifting) to fit models on fluxes provided in a catalog.


* **Fit fluxes produced on-the-fly by a model.**

    * **psf-averager-eazy.cpp**: uses ```fitter-eazy``` and ```psf-averager-photoz``` to run the EAzY model on fluxes generated on-the-fly by the EGG analytical generator.

    * **psf-averager-bpz.cpp**: uses ```fitter-bpz``` and ```psf-averager-photoz``` to run the BPZ model on fluxes generated on-the-fly by the EGG analytical generator.

    * **psf-averager-egg.cpp**: uses the "null" fitter and ```psf-averager-photoz``` to produce only the true PSF metrics and noise-less photometry (no noisy estimates) from SEDs produced on-the-fly by the EGG analytical generator.

    * **psf-averager-binning.cpp**: uses the "null" fitter and ```psf-averager-photoz``` to produce true PSF metrics and estimates obtained from coarsely binned SEDs produced on-the-fly by the EGG analytical generator.

    * **psf-averager-stars.cpp**: uses the ```fitter-star``` and ```psf-averager-star``` to run the simple star fitting model on fluxes and SEDs generated on-the-fly from a stellar library.
