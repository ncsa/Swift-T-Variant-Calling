These are the test details of this pipeline on the stand alone server of CBSB lab at the University of Khartoum. The server has no resource manager in place, so most of the settings are similar to how the pipeline is run on a PC.

The data used for testing here corresponds to a Whole Chromosome Sequencing sample of 30X provided as a training data for the H3ABioNet node [accreditation exercise](http://h3data.cbio.uct.ac.za/assessments/NextGenVariantCalling/practice/)

The options for running the pipeline are provided in the `runfile` (for pipeline options), `sampleinfo` (for what/where the samples are), and `settings.sh` (for swift/t options). Here we are testing with 6 MPI ranks, which can also be seen in the final swift-t logging report: `log.swift_t_run`. 

`Timing.log` is the timing trace from running the swift-t pipeline as per these options, and is typically put in the `OUTPUTDIR/deliverables/docs`. The figure below is a visualization of this log using our [Plotting app](https://github.com/ncsa/Swift-T-Variant-Calling#logging-functionality)

<img src=./media/Plotting_App_server.png width="700">

Finally, the `profile.complete_pipeline.log` file  contains the corresponding resource usage details as extracted via dstat. The figure below visualizes the CPU usage per app. 

<img src=./media/CPU_usage.png width="400"> 
