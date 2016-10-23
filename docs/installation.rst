Instalation
***********

Dependencies
============

There is number of packages that the pipeline requires.

CONTRAfold
----------

Download and install `CONTRAfold <http://contra.stanford.edu/contrafold/download.html>`_.
It might be that you experience an error when compiling CONTRAfold. Something like this:

.. code-block:: bash

    In file included from LBFGS.hpp:52:0,
                     from InnerOptimizationWrapper.hpp:12,
                     from OptimizationWrapper.hpp:12,
                     from Contrafold.cpp:16:
    LBFGS.ipp: In instantiation of ‘Real LBFGS<Real>::Minimize(std::vector<T>&) [with Real = double]’:
    OptimizationWrapper.ipp:260:9:   required from ‘void OptimizationWrapper<RealT>::LearnHyperparameters(std::vector<int>, std::vector<T>&) [with RealT = double]’
    Contrafold.cpp:451:9:   required from ‘void RunTrainingMode(const Options&, const std::vector<FileDescription>&) [with RealT = double]’
    Contrafold.cpp:68:54:   required from here
    LBFGS.ipp:112:105: error: ‘DoLineSearch’ was not declared in this scope, and no declarations were found by argument-dependent lookup at the point of instantiation [-fpermissive]
    LBFGS.ipp:112:105: note: declarations in dependent base ‘LineSearch<double>’ are not found by unqualified lookup
    LBFGS.ipp:112:105: note: use ‘this->DoLineSearch’ instead
    make: *** [Contrafold.o] Error 1



To fix it:

* add -fpermissive flag to CSXXFLAGS in Makefile:

.. code-block:: c

    CXXFLAGS = -O3 -DNDEBUG -W -pipe -Wundef -Winline --param large-function-growth=100000 -Wall -fpermissive

    instead of
    CXXFLAGS = -O3 -DNDEBUG -W -pipe -Wundef -Winline --param large-function-growth=100000 -Wall

* add in Utilities.hpp:

.. code-block:: c

  #include <limits.h>

Jobber
------

Download and setup Jobber python library for workflow managment.

.. code-block:: bash

    pip install Jobber

After installation start the Jobber daemon:

.. code-block:: bash

    $ nohup jobber_server > jobber.log 2>&1 &


.. note::

    If you installed Jobber as user you might not have an access to the jobber_server. By
    default the binary location is $HOME/.local/bin and you have to export it in bash:

    .. code-block:: python

        $ export PATH="$HOME/.local/bin:$PATH"


    or add this statement to .bashrc file.

    jobber_server produces ~/.jobber/jobber.pid file that indicates whether the Jobber is already
    running. If the file exists one cannot start new instance of the jobber_server. This file is
    not clean when jobber_server is killed - only when it was stopped with stop command. Thus,
    after some crash one have to remove this file in order to start jobber_server again.


This will automatically create a ~/.jobber and ~/jobber/log directories and
it will put there config.py and executers.py files. Look at them and adjust
according to your needs.

This should create a jobber.sqlite file next to config.py where jobs will be stored (all in ~/.jobber).
Now you can create pipelines that will be managed with a python script.


To stop the jobber daemon, run following command:

.. code-block:: bash

    $ jobber_server -stop

You can watch and control your jobs and pipelines present in the database using simple we interface.
To launch it type:

.. code-block:: bash

    $ jobber_web

or

.. code-block:: bash

    $ jobber_web --ip Your.IP.addres --port YourPort

.. note::
    If you would like to run snoRNAHybridSearch pipeline locally without DRMAA change executer
    in config.py file from "drmaa" to "local"


BEDTools
--------

Please refer to `BEDTools website <http://bedtools.readthedocs.io/en/latest/>`_ for detailed
installation instructions.

ViennaRNA package
-----------------

Please refer to `ViennaRNA website <http://www.tbi.univie.ac.at/RNA/>`_ for detailed
installation instructions.


SAM Tools
---------

Please refer to `SAM Tools website <http://samtools.sourceforge.net/>`_ for detailed
installation instructions.

Bowtie 2
--------

Please refer to `Bowtie 2 website <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ for detailed
installation instructions.

Python
------

Install required python modules:
 * Jobber (see upper paragraph)
 * drmaa (if you are going to submit it to the cluster)
 * statsmodels
 * pandas
 * BioPython
 * numpy
 * scipy
 * swalign
 * configobj
 * HTSeq
 * MetaProfile
 * ushuffle (one can use this `repo <https://github.com/guma44/ushuffle.git>`_)

For documentation build and not necessary for run:
 * sphinx
 * sphinx-argparse
 * sphinx_rtd_theme


Download
========

The pipeline code is available as a git repository on GitHub or on our website:

.. code-block:: bash

    git clone https://github.com/guma44/snoRNAHybridSearchPipeline.git

    OR

    wget http://www.clipz.unibas.ch/snoRNAchimeras/snoRNAHybridSearchPipeline.tar.gz

In order to run the example and to run pipeline it is neccessary to provide number
of additional files including genome, annotations and snoRNA sequences. Preprepared
files for GRCh37 can be downloaded from our website. If you would like to prepare
your own data it is recomended to look at these files, too:

.. code-block:: bash

    wget http://www.clipz.unibas.ch/snoRNAchimeras/snoRNAHybridSearchData.tar.gz


You can also download whole package including additional data from our website:

.. code-block:: bash

    wget http://www.clipz.unibas.ch/snoRNAchimeras/snoRNAHybridSearch.tar.gz
