#!/bin/sh

julia --project=docs -e 'using Pkg; Pkg.resolve(); using LiveServer; using CRTBPNaturalMotion; servedocs()'