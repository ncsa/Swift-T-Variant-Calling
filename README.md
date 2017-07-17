# Logging functionality

The provided scripts allow you to check out the trace of a successful run of the pipeline. To invoke it, and for the time being, you need R installed in your environment along with the `shiny` package. 

To do so, proceed as follows:

1. Go to the [R-project webpage](http://ftp.heanet.ie/mirrors/cran.r-project.org/), and follow the instructions based on your system
2. Once the step above is completed and R is installed, open a terminal window, type `R`, then proceed as follows:


```
install.packages('shiny')
runGitHub(repo = "jacobrh91/Swift-T-Variant-Calling", ref = "dev-logging",
          subdir = "src/plotting_app" )
```

The first time you run these commands in your system it will also install some libraries for you in case you don't have them already, namely: `lubridate, tidyverse and forcats`.

Once all is done, a webpage should open up for you to actually take a look at your trace files. For a taste of how things look, you may take a look at the sample `Timing.log` file provided [in the repo](https://github.com/jacobrh91/Swift-T-Variant-Calling/tree/dev-logging/src/plotting_app)

To take a look at your own analysis trace, you need to have a copy of this branch first, Run it on you samples, and then find your own `Timing.log` file within `Results_folder_path/delivery/docs`. Simply upload this file, and start using the app.

## Important Notes:

One problem spotted from using the app with 2 samples is that the analysis is done for only one of them (the realignment/recalibration stages are problemetic, where sampleNames get swapped haphazardly, and only one sample gets fully analyzed, which is what the supplied example `Timing.log` file shows - **this needs a closer look**)

It should also be noted that running this pipeline in its current form is expected to be more expensive than normal, due to the manual logging involved. The alternative is to use the native `MPE` library (or equivalent), which requires re-compiling the Swift/T source. This approach is **currently limited at the moment**, but some discussions with the Swift/T team on this is found on their [repo ](https://github.com/swift-lang/swift-t/issues/118)

