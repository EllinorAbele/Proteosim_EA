# Proteosim_EA

A Python package for simulating proteomics experiments — from protein digestion to MS1/MS2 spectrum generation.  

---

# What Does This Project Do?

Proteosim_EA performs a complete **in-silico simulation of proteomic mass spectrometry data**.  
Given protein FASTA sequences, the package generates peptides, precursor spectra (MS1), and fragment spectra (MS2).

The workflow consists of three core stages:

---

## 1. Digestion

Proteins are enzymatically digested into peptides according to a user-selected enzyme (e.g., trypsin).

**Output:**
- List of peptide sequences  

---

## 2. HPLC / LC Simulation*
After digestion, relative LC retention times are predicted for all peptides.
Optionally, the retention time distribution can be visualized as a histogram, and peptides can be filtered based on a selected retention time window.

**Output:**  
- predicted retention times (per peptide)  
- optional: retention time histogram  
- optional: filtered peptide list based on LC time window

---

## 3. MS1 Simulation – Precursor Spectrum

All peptides are converted into precursor ions (typically doubly charged).  
The mass spectrometer simulation then produces an MS1 spectrum containing all peptide precursors.

**Output:**
- Precursor m/z table  
- Simulated MS1 spectrum covering all peptides

---

## 4. MS2 Simulation – Fragmentation

A single peptide is selected for fragmentation (CID-like b/y series).  
The resulting fragment ions are assembled into an MS2 spectrum.

**Output:**
- List of fragment ions  
- Simulated MS2 spectrum for this peptide

---

# Installation

To install this package and its dependencies:

```bash
pip install -r requirements.txt
pip install .

