# MELTS Streamlit GUI — Installation Guide

This package contains a Streamlit-based graphical interface for:
- MELTS liquidus calculations
- Sub-liquidus equilibrium
- GRD08 viscosity
- Mueller et al. (2011) relative viscosity

## 1. Requirements

- Windows 10/11 (64-bit)
- Conda (Anaconda or Miniconda)
- Python 3.10–3.12 recommended
- Streamlit
- MELTSdynamic (included in /app)

## 2. Install the environment

From the top-level directory:

    conda env create -f environment.yml
    conda activate melts_gui

## 3. Run the app

Windows:

    run_app.bat

Mac/Linux:

    ./run_app.sh

Or manually:

    streamlit run app/MELTS_streamlit_GUI.py

## 4. MELTSdynamic

The MELTS engine DLLs are included:
- libalphamelts.dll
- libgsl*.dll
- libxml2.dll
- zlib1.dll
- etc.

These must remain in the same folder as:
- meltsdynamic.py
- meltsengine.py
- meltsstatus.py

No compilation is required.

## 5. Troubleshooting

If MELTS fails to load:
- Ensure DLLs remain next to the Python modules
- Ensure you are using the conda environment
- Ensure Python is 64-bit
