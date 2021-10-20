'''
Created on 21 feb. 2018
Updated on 28 jan 2021
Last update on 12 Feb 2021

@author: thomasgumbricht
'''

# Package application imports

from geoimagine.setup_processes import SetupProcessesRegions, ModisTileCoords, Ease2PolarTileCoords, Ease2GlobalTileCoords

def SetupDefaultRegions(prodDB):
    ''' Setup of global default regions and projection system tiling etc
    '''
    
    ''' DegreeTiles creates a global system of 1deg by 1deg tiles and adds these to the db
         You can only run this command once, then comment it out with a "#"'''
    
    #DegreeTiles(prodDB)

    # DefaultRegions runs a set of scripts installing region categories and global default regions
    DefaultRegions = True
    
    MODIS = False
    
    EASE2 = False
    
    Landsat = False
    
    Sentinel = False
    
    Ancillary = False
    
    Climate = False
    
    '''Link to project file that sets up default regions, arbitrary regions and special regions. 
    '''
    if DefaultRegions:
        
        projFN = 'regions_karttur_setup_20210117.txt'
        
        SetupProcessesRegions('regiondoc', projFN, prodDB)
        
    if MODIS:
        '''Stand alone script that defines the MODIS tile coordinates'''
        
        ModisTileCoords(prodDB)
        
        '''
        exitstr = 'The script ModisTileCoords() produced a shape with all MODIS SIN tiles projected to Geographic coordiantes.\n \
            Copy the shape data sourse: %(fpn)s,\n and edit the xml file for importing this layer.\n \
            Then comment out the "exit" command and re-run the module.' %{'fpn':FPN}
        exit(exitstr)
        '''
        
        #projFN = 'modis_karttur_setup_20210127.txt'
        
        #SetupProcessesRegions('modisdoc', projFN, prodDB)
        
    if EASE2:
        '''Stand alone script that defines EASE grid tile coordinates'''
        
        Ease2PolarTileCoords(prodDB)
        
        Ease2GlobalTileCoords(prodDB)
        
    if Landsat:
        
        '''Link to project file that imports the Landsat WRS system'''
        
        projFN = 'landsat_karttur_setup_20210128.txt'
        
        SetupProcessesRegions('landsatdoc',projFN, prodDB)
    
    if Sentinel:
        
        '''Link to project file that sets up the Sentinel tiling system'''
        
        projFN = 'sentinel_karttur_setup_20210129.txt'
        
        SetupProcessesRegions('sentineldoc',projFN, prodDB)
        
    if MODIS and Sentinel:
        
        '''Stand alone script that links sentinel and modis, requires that all sentinel tiles are in the db'''
       
        pass
    
        #LinkSentineModisTiles()

    if MODIS and Landsat:
        
        pass
    
        #LinkLandsatModisTiles()

    if Sentinel and Landsat:
        
        pass
    
        #LinkSentinelLandsatTiles()

    if Ancillary:
        
        ''' link to project file that imports default ancillary data'''
        
        projFN = 'ancillary_karttur_setup_20180221_0.txt'
        
        SetupProcessesRegions('ancildoc',projFN,prodDB)
    
    if Climate:
        ''' Climate data'''
        
        projFN = 'climate_karttur_setup_20181116.txt'
        
        #Setup('climatedoc',projFN,verbose)
    
def SetupProcesses(prodDB):
    ''' Setupprocesses links to a text file defining Framework processes to define
    '''
    
    relativepath = 'dbdoc'
    
    txtfilename = 'process_karttur_setup_20210104.txt'
    
    SetupProcessesRegions(relativepath, txtfilename, prodDB)
       
if __name__ == "__main__":
    
    '''
    This module should only be run at the very startup of building the Karttur Geo Imagine framework.
    To run, remove the comment "#prodDB" and set the name of your production DB ("YourProdDB")
    '''
    
    # Set the name of the productions db cluster
    # prodDB = 'YourProdDB' #'e.g. postgres or geoimagine
    prodDB = 'geoimagine'
    
    # Setupprocesses links to a text file defining Framework processes to define
    SetupProcesses(prodDB)
    
    # SetupDefaultRegions starts a subroutine with different region processing
    #SetupDefaultRegions(prodDB)
    
    exit(' REACHED THE END! ')