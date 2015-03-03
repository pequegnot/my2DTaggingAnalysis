# Presentation

This is the analysis tools for the 2DTagging (jet flavour) study.

# Setup my2DTaggingAnalysis

## Get the code

Get the code by executing this command:

Using HTTPS:

```bash
myWorkingDir> git clone https://github.com/pequegnot/my2DTaggingAnalysis.git
```

## Setup the environment

Every time you want to work with MttTools, you have to setup the environment: 

```bash
myWorkingDir> cd my2DTaggingANalysis
my2DTaggingAnalysis> source setup_lyoserv_sl6_env.sh
```

# How-to use/run my2DTaggingAnalysis

## Step 1 : generate the tree containing the informations you want to analyse and compare, both for data and MC.

**Code needed :** 

* main code for data: 
`data2012_2DTagging_data.cpp`

* main code for MC :
`data2012_2DTagging_MC.cpp`

* common codes :
`common.h` (code with the different functions common to several main codes)
`binning.h` (contains the flavour binning and number of flavours, i. e. the tagger flavour number associated with the name of the flavour, and the 2D tagging zones binning and number of zones)

* classes used in data and MC main codes (called in common.h) :
`ptBinning`
`alphaBinning`

* QG-Smearing :
`QGSyst.cpp`
`QGSyst.h`
`SystDatabase_doubleMin.txt`

* Makefile to run the code :
`Makefile` 

**Files needed in inputs of the code :**

The rootfiles you need in input must be in `input_rootfile` directory: you have to create it:
```bash
my2DTaggingAnalysis> mkdir input_rootfile 
```

* examples of MC ROOTfiles : (you can only copy the one you need)
input_rootfile/PhotonJet_G_PFlowAK5chs.root
input_rootfile/PhotonJet_QCD_PFlowAK5chs.root (this is QCD EM-enriched)
input_rootfile/PhotonJet_MC_TOT_PFlowAK5chs.root (this is QCD EM-enriched)

* examples of data ROOTfiles :
input_rootfile/PhotonJet_Photon_Run2012_residuals_PFlowAK5chs.root 
(or data without residuals : input_rootfile/PhotonJet_Photon_Run2012_PFlowAK5chs.root)


**Run the code:**

Create an images result directory in the main directory :
```bash
my2DTaggingAnalysis> mkdir images2DTagging
my2DTaggingAnalysis> cd images2DTagging/
images2DTagging> mkdir 2DTaggingZones (will contain plots of 2DTaggingPlans for each flavour and each pt zones)
images2DTagging> mkdir FlavourFractions (will contain plots of flavour jet percentages for each pt zones)
images2DTagging> cd ../
```


Create an « output_rootfile » directory.
```bash
my2DTaggingAnalysis> mkdir output_rootfile
```


Make sure that you are using the right input rootfile for the main codes. You can change the name of your output rootfile (which will be stored in the directory output_rootfile you created).
For MC, choose the sample you want to run on : select the right input rootfile and choose a corresponding output rootfiles (both general and matrix rootfiles. 
For example, if you want to run on MC gamma+jet only, do:

```json
TFile *f=TFile::Open("input_rootfile/PhotonJet_G_PFlowAK5chs.root");
TFile *out_matrix = new TFile("output_rootfile/outputMatrix2DTagging_MC_G.root", "recreate");
TFile *out = new TFile("output_rootfile/output2DTagging_MC_G.root", "recreate");
```


gamma+jet
input_rootfile/PhotonJet_G_PFlowAK5chs.root
output_rootfile/outputMatrix2DTagging_MC_G.root
output_rootfile/output2DTagging_MC_G.root

QCD
input_rootfile/PhotonJet_QCD_PFlowAK5chs.root
output_rootfile/outputMatrix2DTagging_MC_QCD.root
output_rootfile/output2DTagging_MC_QCD.root

Total (G+QCD)
input_rootfile/PhotonJet_MC_TOT_PFlowAK5chs.root
output_rootfile/outputMatrix2DTagging_MC_TOT.root
output_rootfile/output2DTagging_MC_TOT.root


In the Makefile code, enter the name of the code you want to run :
```json
NAME = data2012_2DTagging_data  (to run on data) ;
NAME = data2012_2DTagging_MC  (to run on MC file) ;
```

Compile the script and run the executable, by excuting these commands:
```bash
my2DTaggingAnalysis> make
my2DTaggingAnalysis> ./data2012_2DTagging_MC > output_rootfile/MC_G.txt
```
The last command will also save the output in a `txt` file, usefull to retrieve table with MC flavour fractions for example.




## Step 2 : 2DTagging analysis

### General plot results

```bash
my2DTaggingAnalysis> cd my2DTaggingPlots
```

**Code needed:**
`my2DTaggingPlots.C`

**Rootfiles needed:**
the 2DTagging rootfiles you create in step 1 (data and MC), localised in my2DTaggingAnalysis/output_rootfile

**Run the code:**
In my2DTaggingAnalysis/my2DTaggingPlot directory, you will need an output directory containing your plot results:
```bash
my2DTaggingPlots> mkdir plotsResult
my2DTaggingPlots> cd plotsResult
plotsResult> mkdir fractionEvolution (will contain the evolution of the flavour fractions as a function of gammapt plots)
plotsResult> mkdir RmpfPerZonePerPt (will contain the different MPF responses per 2DTagging zones, filled with data and flavours MC. Each plots are done for each gammapt bin)
plotsResult> mkdir tagger (will contain the data/MC btag CSV and QGLikelihood taggers plots for each gammapt bin)
plotsResult> cd ..
```

Modify the name and localisation of your rootfiles to be analysed. Then run the code with ROOT (you are in my2DTaggingAnalysis /my2DTaggingPlot):
```bash
my2DTaggingPlots> root -l
root> .L ../ptBinning.cpp++
root> .L my2DTaggingPlots.C++
root> my2DTaggingPlots ()
```


### Apply the 2DTaggingMethod with 4x4 matrix (4 2DTaggingZones)

```bash
my2DTaggingPlots> cd ../
my2DTaggingAnalysis> cd my2DTaggingMethod
my2DTaggingMethod> ls
matrix4x4  matrix6x4 (method not implemented yet)
my2DTaggingMethod> cd matrix4x4
```

**Code needed:**
`my2DTaggingMethod_4x4.C`


**Rootfiles needed:**
the 2DTagging rootfiles you create in step 1 (data and MC), localised in my2DTaggingAnalysis/output_rootfile

**Run the code:**
In my2DTaggingAnalysis /my2DTaggingMethod/matrix4x4/ directory, you will need an output directory containing your plot results :
```bash
matrix4x4> mkdir plotsResult
```

Modify the name and localisation of your rootfiles to be analysed. Then run the code with ROOT (you are in my2DTaggingAnalysis /my2DTaggingPlot):
```bash
matrix4x4> root -l
root> .L ../../ptBinning.cpp++
root> .L my2DTaggingMethod.C++
root> my2DTaggingMethod ()
```
