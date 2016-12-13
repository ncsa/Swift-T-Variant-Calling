#!/bin/bash

export PPN=2
swift-t -l -n 2 -r $PWD/pipelinefunctions  VCcallingPipeline.swift --runfile=HgG0.lowcoverage.chr20.parameters-azza
