matTFA
=======

Requirements
------------

You will need to have `Git-LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/pytfa.git /path/to/pytfa
    cd /path/to/pytfa
    git lfs install
    git lfs pull

The scripts have been developed with Matlab 2012a, and CPLEX 12.5 (freely downloadable with the `IBM Academic initiative <https://developer.ibm.com/academic/>`_), and successfully ran on several other versions of both softwares. However, it is important to respect the IBM compatibility specs sheets between Matlab, CPLEX, and the computer OS - available `on IBM's website <https://www.ibm.com/software/reports/compatibility/clarity/index.html>`_

We recommend the stable combination of MATLAB 2016a and CPLEX 12.7 (also freely downloadable from the IBM Academic initiative).


License
=======
The software in this repository is put under a GPLv3.0 licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/matTFA/blob/master/LICENSE>`_ file for more details.

This software uses open source components. These are included in the `ext/ <https://github.com/EPFL-LCSB/matTFA/blob/master/ext>`_ folder. See the `NOTICE <https://github.com/EPFL-LCSB/matTFA/blob/master/ext/NOTICE.rst>`_ for more information on their specific licensing schemes
