This directory includes code to build a custom supernova model in SNEWPY from energy and luminosity files. It currently can create a tarball, but is having issues interfacing 
with SNOwGLoBES generate_time_series() function. If used, it will have to be modified in some way to make it compatible with the base code (something about the way
python and the snowglobes.py file interact with the custom tarball creates value and/or unit errors), but if anyone wants a place to start with, this code is good for that! 
--> only thing that is not compatible is turning tarball into a timeseries, have not tested it in regards to more typical SNEWPY and SNOwGLoBES energy binned functions
REQUIRES miniconda3, globes-3.2.18, snowglobes, SNEWPY 
