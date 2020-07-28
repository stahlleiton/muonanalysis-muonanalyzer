# MuonAnalysis-MuonAnalyzer

Package to run tag/probe ntuples for Muon POG on both AOD and miniAOD format.


## Setup
```
cmsrel CMSSW_10_5_0 
cd CMSSW_10_5_0/src
cmsenv
git cms-init
git clone https://gitlab.cern.ch/cms-muonPOG/muonanalysis-muonanalyzer.git MuonAnalysis/MuonAnalyzer
scram b -j 8
```

## Usage
```
cmsRun MuonAnalysis/MuonAnalyzer/test/run_muonAnalizer_cfg.py
```
