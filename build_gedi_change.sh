#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )

####install requirements packages
conda env update -f ${basedir}/gedi_change_env.yml
#pushd ${HOME}
#source activate gedi_pa_env

conda run --no-capture-output --name gedi_change_env Rscript ${basedir}/install.R
