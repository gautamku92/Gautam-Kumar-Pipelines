#!/bin/bash
#
#$ -N dlc_himem
#$ -l h_vmem=100G
#$ -o /home/tn2283/logs/
#$ -e /home/tn2283/logs/

source /mnt/mfs/cluster/bin/CONDA/etc/profile.d/conda.sh
conda activate /home/tn2283/conda/my-envs/deeplabcut
export DLClight=True
export http_proxy=http://bcp3.cumc.columbia.edu:8080
export https_proxy=http://bcp3.cumc.columbia.edu:8080
export TMPDIR=/mnt/mfs/hgrcgrid/homes/tn2283/tmp

python - <<END
import wx,os,tensorflow
import deeplabcut
import pandas as pd
import numpy as np
from pathlib import Path

main_dir = os.getcwd()

path_config_file = "/mnt/mfs/hgrcgrid/homes/tn2283/$1/config.yaml"
video_file = "'/mnt/mfs/hgrcgrid/homes/tn2283/$1/videos/$2' + '.mp4'"

deeplabcut.create_training_dataset(path_config_file, windows2linux=True)
deeplabcut.train_network(path_config_file, displayiters=50000, saveiters=50000, maxiters=200000)
deeplabcut.evaluate_network(path_config_file, plotting=True)
deeplabcut.analyze_videos(path_config_file,["'/mnt/mfs/hgrcgrid/homes/tn2283/$1/videos/$2.mp4'"
],save_as_csv=True)
deeplabcut.create_labeled_video(path_config_file,["'/mnt/mfs/hgrcgrid/homes/tn2283/$1/videos/$2.mp4'"
])
END
