# Plant cell division and growth analysis

Code developed for the data analysis for (Ref manuscript)

- morphographx_to_tissue.jl: it takes in input a triangular mesh from MorphographX and it returns an input file for Tissue (.init) (computational software for plant simultion).
- cell_wall_growth.jl: it takes in input two .init files (two different time points) and a file containing the parent lables between the two time point segmentations from MorphographX. It returns the growth map of the cell walls.
- shortestPath_beforeDivision_with PGD.jl: it takes in input two .init files and the labels of the dividing cells. It compute the shortest path for the dividing cells and compare it with the actual plane of division.
- ImageJ analysis macro for flourescence analysis (Dr Matthieu Bourdon): macro code for signal quantification in immunostaining images (highly dividing zone vs slower dividing zone). 



