'''
Created on 20 Jan 2021

@author: thomasgumbricht
'''

# Standard library imports

from os import path, makedirs

from sys import exit

# Third party imports

import psycopg2

from base64 import b64encode, b64decode

import netrc

from pyproj import Proj, transform

# Package application imports

from geoimagine.ancillary import ProcessAncillary

from geoimagine.postgresdb import ManageAncillary

from geoimagine.modis import ProcessModis

from geoimagine.postgresdb import ManageProcess, ManageRegion, ManageMODIS

from geoimagine.gis import kt_gis as ktgis

from geoimagine.params import JsonParams

from geoimagine.params.layers import VectorLayer

from geoimagine.postgresdb.easegrid import ManageEASEgrid

from geoimagine.region import IntersectRegions


def DbConnect(db):
    '''
    '''
    # the HOST must exist in the .netrc file in the users home directory
    HOST = 'karttur'

    # Retrieve login and password from the .netrc file
    secrets = netrc.netrc()

    # Authenticate username, account and password
    username, account, password = secrets.authenticators( HOST )

    # Encode the password before sending it
    password = b64encode(password.encode())

    # Create a query dictionary for connecting to the Postgres server
    query = {'db':db, 'user':username, 'pswd':password}

    return query

class PGsession:
    """Connect to postgres server"""

    def __init__(self, query):
        """Connect to selected database"""

        query['pswd'] = b64decode(query['pswd']).decode('ascii')

        conn_string = "host='localhost' dbname='%(db)s' user='%(user)s' password='%(pswd)s'" %query

        self.conn = psycopg2.connect(conn_string)

        self.cursor = self.conn.cursor()

        self.name = 'Setup_process'

    def _DictToSelect(self, queryD):
        '''
        Converts a dictionary to Select statement
        '''
        selectL = []
        for key in queryD:
            #statement = key operator value
            statement = ' %(col)s %(op)s \'%(val)s\'' %{'col':key.replace('#',''), 'op':queryD[key]['op'], 'val':queryD[key]['val']}
            selectL.append(statement)
        self.select_query = "WHERE %(where)s" %{'where':' AND '.join(selectL)}
        return self.select_query

    def _SelectRootProcess(self,queryD):
        '''
        '''
        self.cursor.execute("SELECT rootprocid, minuserstratum FROM process.subprocesses WHERE subprocid = '%(subprocid)s';" %queryD)

        record = self.cursor.fetchone()

        return record

    def _SelectUserCred(self, queryD):
        '''
        '''

        sql = "SELECT userid, usercat, stratum FROM userlocale.users WHERE userid = '%(user)s';" %queryD

        self.cursor.execute(sql)

        self.userid, self.usercat, self.stratum = self.cursor.fetchone()

    def _SelectTractDefRegion(self, queryD):
        '''
        '''
        #First check if this region is itself a defregion

        sql = "SELECT regionid FROM regions.defregions WHERE regionid = '%(tract)s';" %queryD

        self.cursor.execute(sql)

        rec = self.cursor.fetchone()

        if rec != None:

            return (rec[0], 'D')

        sql = "SELECT parentid FROM regions.tracts WHERE tractid = '%(tract)s';" %queryD

        self.cursor.execute(sql)

        rec = self.cursor.fetchone()

        if rec == None:

            return rec

        return (rec[0], 'T')

    def _SelectProcessSystem(self, queryD, paramL):
        ''' Select system for this process
        '''

        queryD['cols'] = " ,".join(paramL)

        sql = "SELECT %(cols)s FROM process.procsys WHERE subprocid = '%(subprocid)s' and system = '%(system)s';" %queryD

        self.cursor.execute(sql)

        record = self.cursor.fetchone()

        if record == None:

            self.cursor.execute("SELECT srcsystem, dstsystem, srcdivision, dstdivisio FROM process.procsys WHERE subprocid = '%(subprocid)s' and system = '*';" %queryD)

            record = self.cursor.fetchone()

        if record == None:

            exitstr = 'No records in _setup_process_class.PGsession.SelectProcessSystem'

            exit(exitstr)

        return dict(zip(paramL,record))

    def _SingleSearch(self,queryD, paramL,  table, schema, pq = False):
        #self._GetTableKeys(schema, table)
        selectQuery = {}
        for item in queryD:

            if isinstance(queryD[item],dict):
                #preset operator and value
                selectQuery[item] = queryD[item]
            else:
                selectQuery[item] = {'op':'=', 'val':queryD[item]}
        wherestatement = self._DictToSelect(selectQuery)
        cols =  ','.join(paramL)
        selectQuery = {'schema':schema, 'table':table, 'select': wherestatement, 'cols':cols}

        query = "SELECT %(cols)s FROM %(schema)s.%(table)s %(select)s" %selectQuery
        if pq:
            print ('SingleSearch query',query)
        self.cursor.execute(query)
        self.records = self.cursor.fetchone()
        return self.records

    def _MultiSearch(self,queryD, paramL, schema, table):
        ''' Select multiple records from any schema.table
        '''

        selectQuery = {}
        for item in queryD:

            if isinstance(queryD[item],dict):
                #preset operator and value
                selectQuery[item] = queryD[item]
            else:
                selectQuery[item] = {'op':'=', 'val':queryD[item]}
        wherestatement = self._DictToSelect(selectQuery)

        if len(paramL) == 1:
            cols = paramL[0]
        else:
            cols =  ','.join(paramL)
        selectQuery = {'schema':schema, 'table':table, 'select': wherestatement, 'cols':cols}
        query = "SELECT %(cols)s FROM %(schema)s.%(table)s %(select)s" %selectQuery

        self.cursor.execute(query)

        self.records = self.cursor.fetchall()

        return self.records

    def _SetSystem(self,system):
        '''
        '''
        self.system = system

    def _Close(self):

        self.cursor.close()

        self.conn.close()

class ProcessProcess:
    """"class for processes defining other processes"""

    def __init__(self, pp):
        """"The constructor requires an instance of the main process, and the json object defining the process to setup"""

        
        
        
        self.pp = pp

        self.verbose = pp.process.verbose

        if self.verbose:

            infostr = '        Processing %s' %(self.pp.process.processid)

            print (infostr)

        self.session = ManageProcess(self.pp.postgresdb.db)
        
        

        if self.pp.process.processid == 'addrootproc':

            if self.verbose:

                print ('            %s' %(self.pp.process.parameters.rootprocid))

            queryD = {'rootprocid':self.pp.process.parameters.rootprocid,
                    'title':self.pp.process.parameters.title,
                  'label':self.pp.process.parameters.label,
                  'creator':self.pp.userproject.userid}

            self.session._ManageRootProcess(self.pp.process, queryD)

        elif self.pp.process.processid == 'addsubproc':

            if self.verbose:

                print ('            %s' %(self.pp.process.parameters.subprocid))

            queryD = {'rootprocid':self.pp.process.parameters.rootprocid,
                      'subprocid':self.pp.process.parameters.subprocid,
                      'title':self.pp.process.parameters.title,
                      'label':self.pp.process.parameters.label,
                      'version':self.pp.process.parameters.version,
                      'minuserstratum':self.pp.process.parameters.minuserstratum,
                      'creator':self.pp.userproject.userid}

            self.session._ManageSubProcess(self.pp.process, queryD)

        else:

            exitstr = 'subprocess %s not defined in manageprocess' %(self.process.processid)

            exit( exitstr )

class ProcessDefaultRegions(IntersectRegions):
    '''
    '''
    def __init__(self, pp ):

        self.pp = pp
        
        self.verbose = self.pp.process.verbose

        self.session = ManageRegion(self.pp.postgresdb.db, self.verbose)
        
        IntersectRegions.__init__(self)

        #direct to subprocess
        if self.pp.process.processid == 'RegionCategories':

            self.session._InsertRegionCat(self.pp.process)

        elif self.pp.process.processid == 'DefaultRegionFromCoords':

            self._DefaultRegionFromCoords()

        elif self.pp.process.processid == 'DefaultRegionFromVector':

            self._DefaultRegFromVec()

                
        elif self.pp.process.processid.lower() == 'linkdefaultregiontiles':

            self._LinkDefaultRegionTiles()
                                       
        else:
            
            
            exitstr = 'No process %s under setup_process.Processregion' %(self.pp.process.processid)
            exit(exitstr)

    def _DefaultRegionRegister(self,layer):
        '''
        '''

        # Get the projection
        projection = ktgis.GetVectorProjection(layer.FPN)

        #Set lonlat projection
        lonlatproj = ktgis.MjProj()

        lonlatproj.SetFromEPSG(4326)

        # Get the boundary
        boundsD = ktgis.GetFeatureBounds(layer.FPN,'REGIONID')

        if len(boundsD) != 1:

            exitstr = 'Default regions must consist on only one (1) feature (polygon or multipolygon): %s' %(layer.FPN)

            exit(exitstr)

        k = list(boundsD)[0]

        layer._SetBounds(projection.epsg,boundsD[k][0], boundsD[k][1], boundsD[k][2], boundsD[k][3] )

        #Get the corners in lonlat
        llD = ktgis.ReprojectBounds(layer.BoundsPtL,projection.proj_cs,lonlatproj.proj_cs)

        queryD = {'regionid': self.pp.process.parameters.regionid,
                  'regionname': self.pp.process.parameters.regionname,
                  'parentid': self.pp.process.parameters.parentid,
                  'regioncat': self.pp.process.parameters.regioncat,
                  'parentcat': self.pp.process.parameters.parentcat,
                  'title': self.pp.process.parameters.title,
                  'label':self.pp.process.parameters.label,
                  'epsg':self.pp.procsys.dstepsg}

        if ' ' in queryD['regionid'] or ' ' in queryD['parentid']:

            exit('regionid or parentid with whuite space in setup_process_process')

        self.session._InsertDefRegion(layer, queryD, boundsD[self.pp.process.parameters.regionid], llD, self.pp.process.overwrite, self.pp.process.delete )

    def _DefaultRegionFromCoords(self):
        '''
        '''

        for locus in self.pp.dstLayerD:

            for datum in self.pp.dstLayerD[locus]:

                for comp in self.pp.dstLayerD[locus][datum]:

                    layer = self.pp.dstLayerD[locus][datum][comp]

                    #The destination region must be forced,this is because the locus to be created did not exists when checking for the default locus

                    layer.locus.locus = self.pp.process.parameters.regionid.lower()

                    layer.locus.path = self.pp.process.parameters.regionid.lower()

                    layer._SetPath()

                    fieldDD = self._SetfieldD()

                    layer.CreateVectorAttributeDef(fieldDD)

                    layer._SetBounds(self.pp.procsys.dstepsg,
                                     self.pp.process.parameters.minx,
                                     self.pp.process.parameters.miny,
                                     self.pp.process.parameters.maxx,
                                     self.pp.process.parameters.maxy)

                    projection = ktgis.MjProj()

                    projection.SetFromEPSG(self.pp.procsys.dstepsg)

                    if not layer._Exists() or self.pp.process.overwrite:

                        ktgis.CreateESRIPolygonPtL(layer.FPN, layer.fieldDefL, layer.BoundsPtL, projection.proj_cs, self.pp.process.parameters.regionid)

                    self._DefaultRegionRegister(layer)


    def _DefaultRegFromVec(self):
        '''
        '''

        # dstLayerD and srcLayerD are almost identical
        for locus in self.pp.dstLayerD:

            for datum in self.pp.dstLayerD[locus]:

                for comp in self.pp.dstLayerD[locus][datum]:

                    srcLayer = self.pp.srcLayerD[locus][datum][comp]

                    if not path.isfile(srcLayer.FPN):

                        exitstr = 'No source layer in _DefaultRegFromVec', srcLayer.FPN

                        exit(exitstr)

                    p = self.pp.process.parameters

                    fieldL = [p.vector_db_id, p.vector_db_name,
                               p.vector_db_category, p.vector_db_parentid,
                               p.vector_db_parentcat, p.vector_db_stratum,
                               p.vector_db_title, p.vector_db_label]

                    fieldD = ktgis.GetFeatureAttributeList(srcLayer.FPN, fieldL, p.vector_db_id)

                    if not fieldD:
                        exit('setup_process_class: fieldD failed in _DefaultRegFromVec')

                    for key in fieldD:

                        # Convert the field data to a dict
                        params = ['regionid', 'regionname', 'regioncat', 'stratum', 'parentid', 'parentcat', 'title', 'label']

                        values = [ str(fieldD[key][p.vector_db_id]).lower().replace(' ', '-'),
                                    str(fieldD[key][p.vector_db_name]),
                                    str(fieldD[key][p.vector_db_category].lower()),
                                    int(fieldD[key][p.vector_db_stratum]),
                                    str(fieldD[key][p.vector_db_parentid]).lower().replace(' ', '-'),
                                    str(fieldD[key][p.vector_db_parentcat].lower()),
                                    str(fieldD[key][p.vector_db_title]),
                                    str(fieldD[key][p.vector_db_label]) ]

                        d = dict(zip(params, values))

                        # Replace the process class parameter with the dict
                        self.pp.process.parameters = lambda:None

                        for k,v in d.items():

                            setattr(self.pp.process.parameters, k, v)

                        fieldDD = self._SetfieldD()

                        regionid = self.pp.process.parameters.regionid

                        #Construct the locus for this region
                        locusD = {'locus':regionid,'path':regionid}

                        # Recreate the composition
                        compDstCopy = self.pp.dstLayerD[locus][datum][comp].comp

                        # Set layerid and prefix to "defreg"
                        compDstCopy.layerid = compDstCopy.prefix = 'defreg'

                        # Set content ot roi (region of interest)
                        compDstCopy.content = 'roi'

                        # Reset the compid
                        compDstCopy._SetCompid()

                        # Recreate the vector Layer
                        dstLayer = VectorLayer(compDstCopy, locusD, self.pp.dstPeriod.datumD[datum])

                        dstLayer.CreateVectorAttributeDef(fieldDD)

                        fieldname = p.vector_db_id

                        valueLL = [[fieldD[key][p.vector_db_id]]]

                        if not dstLayer._Exists() or self.pp.process.overwrite: #or overwrite

                            ktgis.ExtractFeaturesToNewDS(srcLayer.FPN, dstLayer.FPN, fieldname,valueLL, dstLayer.fieldDefL)

                        self._DefaultRegionRegister(dstLayer)

                        '''
                            fieldname = 'REGIONID'

                            #Get the epsg and bounds
                            boundsD = ktgis.GetFeatureBounds(dstLayer.FPN,fieldname)

                            if len(boundsD) != 1:

                                exitstr = 'Default regions must consist on only one (1) feature (polygon or multipolygon): %s' %(dstLayer.FPN)

                                exit(exitstr)

                            projection = ktgis.GetVectorProjection(dstLayer.FPN)

                            k = list(boundsD)[0]

                            bounds = boundsD[k]

                            dstLayer._SetBounds(projection.epsg,boundsD[k][0], boundsD[k][1], boundsD[k][2], boundsD[k][3] )

                            _DefaultRegionRegister(self,dstLayer, projection)

                            #Set lonlat projection

                            lonlatproj = ktgis.MjProj()

                            lonlatproj.SetFromEPSG(4326)

                            #Get the corners in lonlat

                            llD = ktgis.ReprojectBounds(dstLayer.BoundsPtL,projection.proj_cs,lonlatproj.proj_cs)


                            title = label = 'default region %s' %(regionid)

                            query = {'regionname':regionname,'regioncat':regioncat, 'parentid':parentid, 'parentcat':parentcat,'regionid':regionid, 'title':title,'label':label,'epsg':projection.epsg}

                            session._InsertDefRegion(self.process, dstLayer, query, bounds, llD )
                        '''


    def _SetfieldD(self):
        ''' Set the fields for default region layers
        '''

        fieldDD = {}

        fieldDD['REGIONID'] = {'name':'REGIONID', 'type':'string','width':32,
                               'precision':0,'transfer':'constant','source':self.pp.process.parameters.regionid }

        fieldDD['NAME'] = {'name':'NAME', 'type':'string','width':64,
                           'precision':0,'transfer':'constant','source':self.pp.process.parameters.regionname }

        fieldDD['CATEGORY'] = {'name':'CATEGORY', 'type':'string','width':32,
                               'precision':0,'transfer':'constant','source':self.pp.process.parameters.regioncat }

        fieldDD['STRATUM'] = {'name':'STRATUM', 'type':'integer','width':4,
                              'precision':0,'transfer':'constant','source':self.pp.process.parameters.stratum }

        fieldDD['PARENTID'] = {'name':'PARENTID', 'type':'string','width':32,
                               'precision':0,'transfer':'constant','source':self.pp.process.parameters.parentid }

        fieldDD['PARENTCAT'] = {'name':'PARENTCAT', 'type':'string','width':32,
                                'precision':0,'transfer':'constant','source':self.pp.process.parameters.parentcat }

        return fieldDD

def SetupProcessesRegions(docpath, projFN, db):
    '''
    Setup processes
    '''

    srcFP = path.join(path.dirname(__file__),docpath)

    projFPN = path.join(srcFP,projFN)

    # Get the full path to the project text file
    dirPath = path.split(projFPN)[0]

    if not path.exists(projFPN):

        exitstr = 'EXITING, project file missing: %s' %(projFPN)

        exit( exitstr )

    infostr = 'Processing %s' %(projFPN)

    print (infostr)

    # Open and read the text file linking to all json files defining the project
    with open(projFPN) as f:

        jsonL = f.readlines()

    # Clean the list of json objects from comments and whithespace etc
    jsonL = [path.join(dirPath,x.strip())  for x in jsonL if len(x) > 10 and x[0] != '#']

    # Get the user and password for connecting to the db
    query = DbConnect(db)

    # Connect to the Postgres Server
    session = PGsession(query)

    ProcPar = JsonParams(session)

    processL = []

    #Loop over all json files and create Schemas and Tables
    for jsonObj in jsonL:

        infostr = '        reading json file: %s' %(path.split(jsonObj)[1])

        print (infostr)

        processL.append( ProcPar._JsonObj(jsonObj) )

    # Close the db connection for getting processes and user
    session._Close()

    for processD in processL:

        for k in range(len(processD)):

            print('        ', path.split(processD[k]['PP'].jsonFPN)[1] )

            print ('        ',k, processD[k])

            if processD[k]['PP'].rootprocid == 'manageprocess':

                ProcessProcess(processD[k]['PP'])

            elif processD[k]['PP'].rootprocid == 'ManageRegion':

                #ProcessDefaultRegions(db, process, self.procsys, self.userproject, self.userid, self.usercat, self.stratum)
                ProcessDefaultRegions(processD[k]['PP'])

            elif processD[k]['PP'].rootprocid == 'Ancillary':

                session = ManageAncillary(db)

                ProcessAncillary(processD[k]['PP'], session)

                session._Close()

            elif processD[k]['PP'].rootprocid == 'MODISProc':

                session = ManageMODIS(db)

                ProcessModis(processD[k]['PP'], session)

                session._Close()

            else:

                print (processD[k]['PP'].rootprocid)
                print (processD[k]['PP'].subprocid)

def ModisTileCoords(db, verbose = 1):
    ''' Create the MODIS defaut tiling system
    '''

    #Open the db session for MODIS
    session = ManageMODIS(db)

    SINproj = ktgis.MjProj()

    SINproj.SetFromProj4('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs')

    LatLonproj = ktgis.MjProj()

    LatLonproj.SetFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')


    ptL = []

    for lon in range(360):

        ptL.append((lon-180,90))

    for lat in range(180):

        ptL.append((180,-1*(lat-90)))

    for lon in range(360):

        ptL.append((-1*(lon-180),-90))

    for lat in range(180):

        ptL.append((-180,lat-90))


    worldgeom = ktgis.ShapelyPolyGeom(ptL)

    worldgeom.ShapelyToOgrGeom()

    worldgeom.GeoTransform(LatLonproj,SINproj)

    worldgeom.OgrGeomToShapely()

    home = path.expanduser("~")

    tarShpFP =  path.join(path.dirname(__file__),'data')

    if not path.exists(tarShpFP):

        makedirs(tarShpFP)

    FN = 'modtiles-multi_karttur_global_epsg6842.shp'

    tarShpFPN = path.join(tarShpFP,FN)

    fieldDefD = {'type':'string','transfer':'constant','source':'globe','width':8}

    fieldDefL = [ktgis.FieldDef('name',fieldDefD)]

    ktgis.CreateESRIPolygonGeom(tarShpFPN, fieldDefL, worldgeom, SINproj.proj_cs, 'globe')

    # Create a shape file for all individual tiles in SIN proj

    FN = 'modtiles-single_karttur_global_epsg6842.shp'

    tarShpFPN = path.join(tarShpFP,FN)

    tarDS,tarLayer = ktgis.ESRICreateDSLayer(tarShpFPN, SINproj.proj_cs, 'polygon', 'tiles', fieldDefL)

    # Create a shape file for all individual tiles in Geographic coordinates
    FN = 'modtiles_karttur_global_0.shp'

    tarShpFPN = path.join(tarShpFP,FN)

    tarDSLonLat,tarLayerLonLat = ktgis.ESRICreateDSLayer(tarShpFPN, LatLonproj.proj_cs, 'polygon', 'tiles', fieldDefL)

    #create a region with all tiles
    tlen = 20015109.3539999984204769

    tlen /= 18

    for h in range(36):

        minx = tlen*(18-36)+h*tlen

        maxx = minx+tlen

        for v in range(18):

            maxy = tlen*(9-18)+(18-v)*tlen

            miny = maxy-tlen

            ptL = [(minx,maxy),(maxx,maxy),(maxx,miny),(minx,miny)]

            tilegeom = ktgis.ShapelyMultiPointGeom(ptL)

            #convert to ogr
            tilegeom.ShapelyToOgrGeom()

            #write target feature
            tilegeom.GeoTransform(SINproj,LatLonproj)

            tilegeom.OgrGeomToShapely()

            coordL = []

            for point in [ptgeom for ptgeom in tilegeom.shapelyGeom]:

                coordL.extend([list(point.coords)[0][0],list(point.coords)[0][1]])

            ullon, ullat, urlon, urlat, lrlon, lrlat, lllon, lllat = coordL

            tilepoly = ktgis.ShapelyPolyGeom([(minx, maxy), (maxx, maxy), (maxx, miny), (minx,miny)])

            #Test if this tile is inside the globe
            if tilepoly.shapelyGeom.intersects(worldgeom.shapelyGeom):

                if h < 10:

                    htile = 'h0%s' %(h)

                else:

                    htile = 'h%s' %(h)

                if v < 10:

                    vtile = 'v0%s' %(v)

                else:

                    vtile = 'v%s' %(v)

                hvtile = '%s%s' %(htile,vtile)

                polytilegeom = ktgis.ShapelyPolyGeom(ptL)

                polytilegeom.ShapelyToOgrGeom()

                fieldDefD = {'type':'string','transfer':'constant','source':hvtile,'width':8}

                fieldDefL = [ktgis.FieldDef('name',fieldDefD)]

                #create target feature
                tarFeat = ktgis.ogrFeature(tarLayer)

                tarFeat.CreateOgrFeature(polytilegeom.ogrGeom, fieldDefL)

                if h == 17:

                    pass

                else:

                    #to be correct 5 points are needed and also the lat must be fitted
                    if h < 18 and ullon > 0:

                        ullon = -180

                    if h < 18 and lllon > 0:

                        lllon = -180

                    if h < 18 and urlon > 0:

                        urlon = -180

                    if h < 18 and lrlon > 0:

                        lrlon = -180

                    if h > 18 and urlon < 0:

                        urlon = 180

                    if h > 18 and lrlon < 0:

                        lrlon = 180

                    if h > 18 and ullon < 0:

                        ullon = 180

                    if h > 18 and lllon < 0:

                        lllon = 180

                    if hvtile == 'h24v01':

                        urlon = 180

                    if hvtile == 'h24v16':

                        lrlon = 180

                    if hvtile == 'h11v01':

                        ullon = -180

                    if hvtile == 'h11v16':

                        lllon = -180

                if ullon > urlon:

                    print ('ERROR','ullon > urlon',hvtile,ullon,urlon)

                if lllon > lrlon:

                    print ('ERROR','lllon > lrlon',hvtile, lllon, lrlon)

                #
                polytilegeom = ktgis.ShapelyPolyGeom([(ullon, ullat), (urlon, urlat), (lrlon, lrlat), (lllon,lllat)])

                polytilegeom.ShapelyToOgrGeom()

                #polytilegeom.GeoTransform(SINproj,LatLonproj)

                #create target feature
                tarLonLatFeat = ktgis.ogrFeature(tarLayerLonLat)

                tarLonLatFeat.CreateOgrFeature(polytilegeom.ogrGeom, fieldDefL)

                west,south,east,north = polytilegeom.shapelyGeom.bounds

                session._InsertModisTileCoord(hvtile,h,v,
                        round(minx,8), round(maxy,8), round(maxx,8), round(miny,8),
                        round(west,8), round(south,8), round(east,8), round(north,8),
                        round(ullat,8), round(ullon,8), round(lrlon,8), round(lrlat,8),
                        round(urlon,8), round(urlat,8), round(lllon,8), round(lllat,8))

                query = {'system':'system','table':'regions','h':h,'v':v,'hvtile':hvtile,'regionid':'global','regioncat':'global','regiontype':'default','delete':False}

                session._InsertModisRegionTile(query)

    tarDS.CloseDS()

    tarDSLonLat.CloseDS()

    session._Close()

    if verbose > 1:

        print ('Check the shape file',tarShpFPN)

    return (tarShpFPN)


def Ease2PolarTileCoords(db, verbose = 1):
    ''' Create the Ease2 polar tiling system
    '''

    eD = {'ease2n':6931, 'ease2s':6932}

    for easegrid in ['ease2n','ease2s']:

        latlonProj = Proj('epsg:4326') # 4326 represents geographic coordinates

        projstr = 'epsg:%(e)d' %{'e':eD[easegrid]}
        # Set the target projection (EASE-grid)
        easeProj = Proj(projstr) # 6933 represents the global/tropial EASE grid

        session = ManageEASEgrid(db)

        home = path.expanduser("~")

        Ease2proj = ktgis.MjProj()

        Ease2proj.SetFromEPSG(eD[easegrid])

        # Create a shape file for all individual tiles in Geographic coordinates
        FN = '%(e)stiles-multi_karttur_epsg%(e)s.shp' %{'e': eD[easegrid]}

        tarShpFPN = path.join(home,FN)

        fieldDefD = {'type':'string','transfer':'constant','source':'globe','width':8}

        fieldDefL = [ktgis.FieldDef('name',fieldDefD)]

        # Create a shape file for all individual tiles in SIN proj
        tarDS,tarLayer = ktgis.ESRICreateDSLayer(tarShpFPN, Ease2proj.proj_cs, 'polygon', 'tiles', fieldDefL)

        # Define the side of a tile
        tileside = 900000

        # Set initial maxx
        maxx = -9000000

        for x in range(20):

            maxx += tileside

            minx = maxx-tileside

            maxy = -9000000

            for y in range(20):

                maxy += tileside

                miny = maxy-tileside

                ptL = [(minx,maxy),(maxx,maxy),(maxx,miny),(minx,miny)]

                tilegeom = ktgis.ShapelyPolyGeom(ptL)

                #convert to ogr
                tilegeom.ShapelyToOgrGeom()

                if x < 10:

                    xtile = 'x0%s' %(x)

                else:

                    xtile = 'x%s' %(x)

                if y < 10:

                    ytile = 'y0%s' %(y)

                else:

                    ytile = 'y%s' %(y)

                xytile = '%s%s' %(xtile,ytile)

                fieldDefD = {'type':'string','transfer':'constant','source':xytile,'width':8}

                fieldDefL = [ktgis.FieldDef('name',fieldDefD)]

                #create target feature
                tarFeat = ktgis.ogrFeature(tarLayer)

                tarFeat.CreateOgrFeature(tilegeom.ogrGeom, fieldDefL)


                west,south,east,north = tilegeom.shapelyGeom.bounds

                corners = ['ul','ur','lr','ll']

                llD = {}
                for z, pt in enumerate(ptL):

                    lat,lon = transform(easeProj, latlonProj, pt[0], pt[1])
                    key = '%(c)slat' %{'c':corners[z]}
                    llD[key] = round(lat,5)

                    key = '%(c)slon' %{'c':corners[z]}
                    llD[key] = round(lon,5)

                # Write tile to db
                session._InsertTileCoord(easegrid,xytile,x,y,round(minx,2),round(maxy,2),round(maxx,2),
                                    round(miny,2),round(west,2),round(south,2),round(east,2),round(north,2),
                                    llD['ullat'],llD['ullon'],llD['lrlat'],llD['lrlon'],llD['urlat'],llD['urlon'],llD['lllat'],llD['lllon'])

                query = {'system':easegrid,'easegrid': easegrid, 'table':'regions','xtile':x,'ytile':y,'xytile':xytile,'regionid':'global','regioncat':'global','regiontype':'default','delete':False}

                session._InsertRegionTile(query)

        tarDS.CloseDS()

        if verbose > 0:

            print ('Check the shape file',tarShpFPN)

def Ease2GlobalTileCoords(db, verbose = 1):
    ''' Create the Ease2 global tiling system
    '''

    easegrid = 'ease2t'

    latlonProj = Proj('epsg:4326') # 4326 represents geographic coordinates

    # Set the target projection (EASE-grid)
    easeProj = Proj('epsg:6933') # 6933 represents the global/tropial EASE grid

    session = ManageEASEgrid(db)

    home = path.expanduser("~")

    Ease2proj = ktgis.MjProj()

    Ease2proj.SetFromEPSG(6933)

    # Create a shape file for all individual tiles in Geographic coordinates
    FN = '%(e)stiles-multi_karttur_epsg%(e)s.shp' %{'e': '6933'}

    tarShpFPN = path.join(home,FN)

    fieldDefD = {'type':'string','transfer':'constant','source':'globe','width':8}

    fieldDefL = [ktgis.FieldDef('name',fieldDefD)]

    # Create a shape file for all individual tiles in SIN proj
    tarDS,tarLayer = ktgis.ESRICreateDSLayer(tarShpFPN, Ease2proj.proj_cs, 'polygon', 'tiles', fieldDefL)

    # Define the side of a tile

    tileside = 936837.98

    xtileside = 1441284.79

    ytileside = 936837.98-180163.01

    # Set initial maxx
    maxx = -17367530.45

    for x in range(36):

        if x == 0 or x == 35:

            maxx += xtileside

        else:

            maxx += tileside

        if x == 0 or x == 35:

            minx = maxx-xtileside

        else:

            minx = maxx-tileside

        # Maxy is not the real maxy, the last rows of data will instead be omitted to keep the dimensions
        maxy = -7314540.83

        for y in range(1,17):

            if y == 1 or y == 16:

                maxy += ytileside

                miny = maxy-ytileside

            else:

                maxy += tileside

                miny = maxy-tileside


            print (x,y,minx,miny,maxx,maxy)

            ptL = [(minx,maxy),(maxx,maxy),(maxx,miny),(minx,miny)]

            tilegeom = ktgis.ShapelyPolyGeom(ptL)

            #convert to ogr
            tilegeom.ShapelyToOgrGeom()

            if x < 10:

                xtile = 'x0%s' %(x)

            else:

                xtile = 'x%s' %(x)

            if y < 10:

                ytile = 'y0%s' %(y)

            else:

                ytile = 'y%s' %(y)


            xytile = '%s%s' %(xtile,ytile)

            fieldDefD = {'type':'string','transfer':'constant','source':xytile,'width':8}

            fieldDefL = [ktgis.FieldDef('name',fieldDefD)]

            #create target feature
            tarFeat = ktgis.ogrFeature(tarLayer)

            tarFeat.CreateOgrFeature(tilegeom.ogrGeom, fieldDefL)

            # Extract bounds
            west,south,east,north = tilegeom.shapelyGeom.bounds

            # Reproject to lat/lon

            # Transform the coordinate point
            corners = ['ul','ur','lr','ll']
            llD = {}
            for z, pt in enumerate(ptL):

                lat,lon = transform(easeProj, latlonProj, pt[0], pt[1])
                key = '%(c)slat' %{'c':corners[z]}
                llD[key] = round(lat,5)

                key = '%(c)slon' %{'c':corners[z]}
                llD[key] = round(lon,5)

            # Write tile to db
            session._InsertTileCoord(easegrid,xytile,x,y,round(minx,2),round(maxy,2),round(maxx,2),
                                    round(miny,2),round(west,2),round(south,2),round(east,2),round(north,2),
                                    llD['ullat'],llD['ullon'],llD['lrlat'],llD['lrlon'],llD['urlat'],llD['urlon'],llD['lllat'],llD['lllon'])

            query = {'system':easegrid,'easegrid': easegrid, 'table':'regions','xtile':x,'ytile':y,'xytile':xytile,'regionid':'global','regioncat':'global','regiontype':'default','delete':False}

            session._InsertRegionTile(query)

    tarDS.CloseDS()

    if verbose > 0:

        print ('Check the shape file',tarShpFPN)

def Ease2GlobalTileCoordsOld(db, verbose = 1):
    ''' Create the Ease2 polar tiling system
    '''

    session = ManageEASEgrid(db)

    home = path.expanduser("~")

    Ease2proj = ktgis.MjProj()

    Ease2proj.SetFromEPSG(6933)

    # Create a shape file for all individual tiles in Geographic coordinates
    FN = '%(e)stiles-multi_karttur_epsg%(e)s.shp' %{'e': '6933'}

    tarShpFPN = path.join(home,FN)

    fieldDefD = {'type':'string','transfer':'constant','source':'globe','width':8}

    fieldDefL = [ktgis.FieldDef('name',fieldDefD)]

    # Create a shape file for all individual tiles in SIN proj
    tarDS,tarLayer = ktgis.ESRICreateDSLayer(tarShpFPN, Ease2proj.proj_cs, 'polygon', 'tiles', fieldDefL)

    # Define the side of a tile
    #tileside = 720644.42
    tileside = 1801611.043

    # Set initial maxx
    #maxx = -17367530.45
    maxx = -17115304.904

    for x in range(19):

        maxx += tileside

        minx = maxx-tileside

        # Maxy is not the real maxy, the last rows of data will instead be omitted to keep the dimensions
        # maxy = -7090347.50497537
        maxy = -7206444.16748768

        for y in range(8):

            maxy += tileside

            miny = maxy-tileside

            print (x,y,minx,miny,maxx,maxy)

            ptL = [(minx,maxy),(maxx,maxy),(maxx,miny),(minx,miny)]

            tilegeom = ktgis.ShapelyPolyGeom(ptL)

            #convert to ogr
            tilegeom.ShapelyToOgrGeom()

            if x < 10:

                xtile = 'x0%s' %(x)

            else:

                xtile = 'x%s' %(x)

            if y < 10:

                ytile = 'y0%s' %(y)

            else:

                ytile = 'y%s' %(y)

            xytile = '%s%s' %(xtile,ytile)

            fieldDefD = {'type':'string','transfer':'constant','source':xytile,'width':8}

            fieldDefL = [ktgis.FieldDef('name',fieldDefD)]

            #create target feature
            tarFeat = ktgis.ogrFeature(tarLayer)

            tarFeat.CreateOgrFeature(tilegeom.ogrGeom, fieldDefL)

            # Write to db
            west,south,east,north = tilegeom.shapelyGeom.bounds

            session._InsertTileCoord('ease2t',xytile,x,y,minx,maxy,maxx,miny,west,south,east,north)

            query = {'system':'system','easegrid':'ease2t', 'table':'regions','x':x,'y':y,'xytile':xytile,'regionid':'global','regioncat':'global','regiontype':'default','delete':False}

            session._InsertRegionTile(query)

    tarDS.CloseDS()

    if verbose > 0:

        print ('Check the shape file',tarShpFPN)
