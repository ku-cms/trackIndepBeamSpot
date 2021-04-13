# Tracking Independent Beamspot in CMS

This repository contains python scripts used for reading layer 1 pixel clusters and fitting distriubtions based on the clusters.

## Reading Clusters

The script readClusters.py should be run in python 2 and requires these packages:
- numpy
- pyROOT

You can run readClusters.py like this:
```
cd python
python readClusters.py
```

However, the script takes a long time to finish, especially when running over a large number of events.
The command "nohup" can be used so that the script can continue running even if you logout of your session.
Here we chose the log file name "readClusters_SingleMuon_v1.log."
```
cd python
mkdir logs
nohup python readClusters.py > logs/readClusters_SingleMuon_v1.log 2>&1 &
```

You can watch the log with tail (-f for follow as the log file grows, and -n for the max number of lines).
```
tail -f -n 100 logs/readClusters_SingleMuon_v1.log
```

You can see your jobs with this command (only before logging out):
```
jobs -l
```

You can find your process with this command (even after logging out and logging back in): 
```
ps aux | grep readClusters.py
```

When the script is finished, you should see the run time at the end of the log file, and the process should no longer be running. You can check that you don't find the process anymore when you run
```
ps aux | grep readClusters.py
```

The script will create a .npy file when it finished running.

## Fitting

The script plotOccupancy.py should be run in python 3 and requires these packages:
- numpy
- matplotlib
- iminuit


Creating virtual environment for python (here we picked the name python3_env):
```
python3 -m venv python3_env
```

To activate:
```
source python3_env/bin/activate
```

To check python path and version:
```
which python
python --version
```

Install packages after activating environment:
```
pip install numpy
pip install matplotlib
pip install iminuit
```

To deactivate:
```
deactivate
```

The script plotOccupancy.py expects these directories to exists within the trackIndepBeamSpot directory (not in the python directory):
```
cd <full_path>/trackIndepBeamSpot
mkdir data
mkdir plots
```

Then you can copy .npy files to the data directory and add the file names to plotOccupancy.py. 

You will need to activate the virtual environment before running plotOccupancy.py.
The activation can be done from any directory, and it will be maintained when changing directories. 
You can also use a full path to activate.

```
source <full_path>/python3_env/bin/activate
python python/plotOccupancy.py
```

This script will save plots as pngs in the plots directory.

