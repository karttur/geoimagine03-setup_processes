"""
setup_processes
==========================================

Package belonging to Karttur´s GeoImagine Framework.

Author
------
Thomas Gumbricht (thomas.gumbricht@karttur.com)

"""
from .version import __version__, VERSION, metadataD

#from .setup_process_class import PGsession, ProcessDefaultRegions

from .setup_process_process import SetupProcessesRegions, ModisTileCoords, Ease2PolarTileCoords, Ease2GlobalTileCoords

__all__ = ['PGsession', 'ProcessDefaultRegions','SetupProcessesRegions','ModisTileCoords',
           'Ease2NTileCoords','Ease2NTileCoords6391']