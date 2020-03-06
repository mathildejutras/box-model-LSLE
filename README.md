# box-model-LSLE
Box Model of the Lower St Lawrence Estuary

This repository contains the python code and GUI software for a box model of nutrient cycling in the Lower St. Lawrence Estuary (LSLE) developed at McGill University by Mathilde Jutras, Alfonso Mucci, Bjorn Sundby, Yves Gratton and Sergei Katsev. The model is presented in the following paper: Jutras, M., et al., (Under review), Nutrient cycling in the Lower St. Lawrence Estuary: response to environmental perturbations, Estuarine, Coastal and Shelf Science.

Here are the steps to install the GUI program on Linux or Windows:

1) Make sure all the required python packages are installed:
- pyinstaller
- Tkinter
- tkMessageBox
- numpy

2) Download the box_model_gui_science.py file.

3) In the directory of the file, run pyinstaller --onefile box_model_gui_science.py

This should create an executable file in the dist directory. To launch the software, double click on the box_model_gui_science icon.

Contact: Mathilde Jutras, mathilde.jutras@mail.mcgill.ca
