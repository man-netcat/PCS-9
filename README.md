# PCS-9

Make sure that you have the following pip3 packages installed: numpy, matplotlib, skimage, pillow

To run tests, run the file `run_tests.sh` in the source directory with your favourite shell (preferably bash)

The results will appear in some newly made directories for velocity and density each. To manually run the simulation, first run `gen_bifurcation.py` with parameters of your choice (run `gen_bifurcation.py` with the flag `--help` for parameter explanations) and then run `blood_flow.py` with parameters of your choice (run `blood_flow.py` with the flag `--help` for parameter explanations)

An example of running the code would be the following:

```
python3 gen_bifurcation.py
python3 blood_flow.py out/bifurcation_0.5000.png video
```
