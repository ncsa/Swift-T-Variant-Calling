Installation
------------

Dependencies
~~~~~~~~~~~~

First, you need Swift/T installed in your system. Depending on your system, the instructions below will guide you through the process:

 http://swift-lang.github.io/swift-t/guide.html#_installation

Next, depending on the analysis step you like, you also need the installation path of the following tools in your system:

+-------------------+-----------------------------------------------------------------------------------+
|     **Step**      |   **Tool options**                                                                |
+===================+===================================================================================+
| Alignment         | `Bwa mem <https://github.com/lh3/bwa>`__ or `Novoalign <http://novocraft.com/>`__ |
+-------------------+-----------------------------------------------------------------------------------+
| Sorting           | `Novosort <http://novocraft.com/>`__                                              |
+-------------------+-----------------------------------------------------------------------------------+
|                   | `Samblaster <https://github.com/GregoryFaust/samblaster>`__,                      | 
+                   +-----------------------------------------------------------------------------------+
| Marking Duplicates| `Novosort <http://novocraft.com/>`__, or                                          | 
+                   +-----------------------------------------------------------------------------------+
|                   | `Picard <https://broadinstitute.github.io/picard/>`__                             | 
+-------------------+-----------------------------------------------------------------------------------+
| IndelRealignment  |                                                                                   |
+-------------------+                                                                                   +
| BaseRecalibration |                                                                                   |
+-------------------+                                                                                   +
| Variant Calling   | `GATK <https://software.broadinstitute.org/gatk/download/>`__                     |
+-------------------+                                                                                   +
| Joint Genotyping  |                                                                                   |
+-------------------+-----------------------------------------------------------------------------------+
| Miscellaneous     | `Samtools <http://samtools.github.io/>`__, and                                    |
+                   +-----------------------------------------------------------------------------------+
|                   | `Novosort <http://novocraft.com/>`__                                              |
+-------------------+-----------------------------------------------------------------------------------+


Workflow Installation
~~~~~~~~~~~~~~~~~~~~~

Simply, clone the repository::

 git clone https://github.com/ncsa/Swift-T-Variant-Calling/

Additionally, you may need ``R`` installed along with the following packages ``shiny``, ``lubridate``, ``tidyverse`` and ``forcats``. Detailed instructions are on the Logging functionality section here :doc:`UserGuide`



