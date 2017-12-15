Troubleshooting
---------------

**General Troubleshooting Tips**

Regardless of the platform, one can use the following environmental
variables to better debug the workflow:

``ADLB_DEBUG_RANKS=1`` One can see if the processes are spread across
the nodes correctly

``TURBINE_LOG=1`` Makes the Swift-T log output very verbose
``TURBINE_LOG_FILE=<filePath>`` Changes the Swift-T log output from
StdOut to the file of choice

More debug info can be found
`here <http://swift-lang.github.io/swift-t/guide.html>`__

FAQs
~~~~~

-  The pipeline seems to be running, but then prematurely stops at one
   of the tools?
- * Solution: make sure that all tools are specified in your runfile up
   to the executable itself (or the jar file if applicable)

-  The realignment/recalibration stage produces a lot of errors or
   strange results?
-  * Solution: make sure you are preparing your reference and extra files
   (dbsnp, 1000G,...etc) according to the guidelines in the `Data
   Preparation <#data-preparation>`__ section 
 
-  Things that should be running in parallel appear to be running
   sequencially
-  * Solution: make sure you are setting the ``-n`` flag to a value at
   least one more than ``PROGRAMS_PER_NODE`` \* ``NODES``, as this
   allocates processes for Swift/T itself to run on 
 
-  The job is killed as soon as BWA is called?
-  * Solution: make sure there is no space in front of ``BWAMEMPARAMS``
   -  DO-THIS: ``BWAMEMPARAMS=-k 32 -I 300,30``
   -  NOT-THIS: ``BWAMEMPARAMS= -k 32 -I 300,30`` 
 
-  I'm not sure how to run on a cluster that uses torque as a resource
   manager?
-  * Solution: Clusters are typically configured to kill head node jobs that run
   longer than a few minutes, to prevent users from hogging the head
   node. Therefore, you may qsub the initial job, the swift-t command
   with its set variables, and it will qsub everybody else from its
   compute node. 
