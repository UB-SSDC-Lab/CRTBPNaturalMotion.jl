#!/bin/sh

julia --project=docs -e 'using LiveServer; using CRTBPNaturalMotion; servedocs()'