"""
Description:
    Extract data from modem output files (.dat and .rho) and produce csv files for further visualization
    The output look like: StationName, Lat, Long, X, Y, Z, Log(Resistivity)
    where (X,Y,Z) are relative distances in meters from the mesh's origin.
    Projection/Coordinate system must be known in order to associate (Lat, Long) to (X, Y)
    
References: 
    https://gajira.atlassian.net/browse/ALAMP-49

CreationDate:   8/09/2017
Developer:      fei.zhang@ga.gov.au

LastUpdate:     8/09/2017   FZ
LastUpdate:     dd/mm/yyyy  Who
"""



