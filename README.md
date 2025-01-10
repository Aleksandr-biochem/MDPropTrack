# MDPropTrack tools

MDPropTrack is a Python3 mini-library designed to facilitate quick property extraction from Molecular Dynamics trajectories, time series plotting and assessment of convergence.  

**Contents:**
1. [What is it convenient for?](#sec1) </br> 
2. [Installation](#sec2) </br>
3. [Import and use](#sec3) </br>
4. [Command line tool usage](#sec4) </br>
5. [Questions and feedback](#sec5) </br>

<a name="sec1"></a>
## 1.  What is it convenient for?

- Quickly read and plot data from GROMACS edr files.

- Compute and plot properties from trajectories: radius of gyration, RMSD, area per lipid, bilayer thickness, lipid tail order parameter. Configure calculation of custom properties.

- Assess convergence using time series autocorrelation.

- Plot sequential simulation steps together.

- Save and load DataFrames.

<a name="sec1"></a>
## 2.  Installation

Clone this repository to your local system:

```
git clone git@github.com:Aleksandr-biochem/MDPropTrack.git

# move to repository folder
cd MDPropTrack
```

To install or update MDPropTrack into your environment use pip:

```
pip install .
```

Now you should be able to call the command line tool:

```
mdpt -h
```

File `requirements.txt` lists the exact versions of the dependencies used during development.

```
# download libraries using pip
pip install -r requirements.txt

# after that run
pip install .
```

<a name="sec3"></a>
## 3. Import and use

The full potential of MDPropTrack is unlocked, when it's used as a python library. See `tutorial.ipynb` for all the details. Here is a quick start snippet:

```
from MDPropTrack.analysis import PropertyAnalyser

# create analyser for your system and get properties from an edr file
pa = PropertyAnalyser(edr='path/to/simulation.edr')

# call the following method to extract properties into DataFrame
pa.extract_properties()

# look at the DataFrame
pa.data

# plot extracted properties as time series
pa.plot()

# or plot autocorrelation to assess convergence
pa.plot(plot_convergence=True)
```

<a name="sec4"></a>
## 4. Command line tool usage

MDPropTrack is also equipped with a command line tool `mdpt` that will read your input files and generate a plot for the requested properties. `tutorial.ipynb` will help you to better understand the `mdpt` arguments. Examples:


```
# see help
mdpt -h

# read a edr file and plot a default list of properties 
mdpt -e file.edr -o plot.jpg 

# read 2 trajectories as inconsequent steps 
# calculate backbone RMSD for chain A and Rg for protein
# use time in ps when plotting
mdpt -t file1.xtc,file2.xtc -p file1.tpr -i -tu ps --rg 'protein' --rmsd 'backbone,backbone and chain A' -o plot.jpg 

# read edrs and trajectories
# center the DPPC bilayer
# calculate some MARTINI lipid properties with stride 10
# plot convergence of the properties
mdpt -e file1.edr -t file1.xtc,file2.xtc -p file1.tpr --center_group 'DPPC' --step 10 --apl 'name PO4' --thickness 'name PO4' --scc 'name ??A,name ??B' -o plot.jpg --plot_convergence 
```

<a name="sec5"></a>
## 5. Questions and feedback

If you have a question, which does not seem to be answered in this manual or tutorial, or if you want to report an issue, please, do so in the `Issues` tab at the GitHub page of this repository. Your feedback is very much appreciated!


