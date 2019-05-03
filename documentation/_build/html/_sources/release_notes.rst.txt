Release Notes
=============

1.1
---

Small changes for the plotting and some performance optimization for the montecarlo simulation.

Plotting
~~~~~~~~
Changed from the plt interface to the fig interface for all plots.

Montecarlo simulation
~~~~~~~~~~~~~~~~~~~~~

Added multiprocessing for the montecarlo simulation. :code:`__init__()` now takes a new keyword argument :code:`threads`.
:code:`threads` specifies how many threads should be used for the montecalo simulation, the default value is 1.
The function :code:`run_montecarlo_simulation` also has the keyword :code:`threads`, to specify the number us used threads.
If no value is set the default value, set in the :code:`__init__()` will be used.

Examples
~~~~~~~~

The river wheaton example as been updated to show the usage of the multiprocessing for the montecarlo simulation.