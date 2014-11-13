run command's - 
	python final_model windowentropy
	python final_model plainentropy
	
NOTE:
	- Runs on python >= 2.7.2 and requires packages numpy and matplotlib 
	- Because of large data processing, running time may take few minutes
	- All input files are in data/
	- Every run outputs a numpy plot and appends the points in the follwing format to output_windowentropy/output_plainentropy file
			real x-points
			real y-points

			model x-points
			model y-points
	