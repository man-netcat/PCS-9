# PCS-9


Make sure that you have the following pip3 packages installed: `numpy, matplotlib, scikit-image, pillow`
Also check if you installed `ffmpeg`


To manually run the simulation, first run `gen_bifurcation.py` with parameters of your choice (run `gen_bifurcation.py` with the flag `--help` for parameter explanations) and then run `blood_flow.py` with parameters of your choice (run `blood_flow.py` with the flag `--help` for parameter explanations)


An example of running the code would be the following:
```
python3 gen_bifurcation.py
python3 blood_flow.py out/bifurcation_0.5000.png video
```
You will get results equal to example.mp4 and example.png in the source directory.


To run all tests, run the file `run_tests.sh` in the source directory with your favourite shell (preferably bash)
The results will appear in some newly made directories for velocity and density each. This script creates for all geometries the density and velocity plots and videos. This  takes however a very long time so we recommend to just test the above example.
