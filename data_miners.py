from astropy.time import Time
import numpy as np
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy import units as u
from astropy.coordinates import SkyCoord

import urllib
from bs4 import BeautifulSoup


"""
Name: kws_grab()
Description: Download data for a given object in the
KWS (Kamogata/Kiso/Kyoto Wide-field Survey)
Inputs:
object: name of object. It is resolved using Simbad.
coords: list [RA, Dec] in decimal degrees.Only used if object not specified
name: Name you would like the file saved under. If not specified, Name
will default to object name. This is either the given name or the MAIN_ID
found in Simbad.
@todo add bands
"""


def kws_grab(object=None, coords=None, name=None):
    # check for object and retrieve if not given
    if not object:
        RA = coords[0]
        Dec = coords[1]
        result_table = Simbad.query_region(SkyCoord(RA, Dec,\
                       unit=(u.deg, u.deg), frame='icrs'), radius='0d0m2s')
        if result_table:
            object = result_table[0]['MAIN_ID']#.decode("ascii")
        else:
            print("No object found at RA = {} Dec = {}".format(RA, Dec))
            return
    #build and send url
    url = 'http://kws.cetus-net.org/~maehara/VSdata.py'

    values = {'object': object,
              'resolver': 'simbad',
              'p_band': 'All',
              'plot': '0' }


    data = urllib.parse.urlencode(values)
    data = data.encode('ascii') # data should be bytes

    #create request
    req = urllib.request.Request(url, data)

    with urllib.request.urlopen(req) as response:
        the_page = response.read()

    #grab table data
    soup = BeautifulSoup(the_page, 'html.parser')
    tables = soup.find_all('table')
    if len(tables) > 0:
        table = tables[0]

        table_data = []

        for row in table.find_all('tr'):
            row_data = []
            columns = row.find_all('td')
            for column in columns:
                row_data.append(column.get_text())

            table_data.append(row_data)

        header = np.array(table_data[0])

        table_data.pop(0)

        table_data = np.vstack(table_data)

        #convert times to JD

        times = []
        for x in range(len(table_data)):
            row = table_data[:,1][x].split(' ')
            time = row[0]+'T'+row[1]
            times.append(time)

        time_jd = Time(times).jd
        #create table data as structured array
        table_data_rev = np.column_stack((time_jd.astype(float),table_data[:,2].astype(float), table_data[:,3].astype(float), table_data[:,4]))
        table_data_rev = list(zip(*table_data_rev.T))

        dt = [("time","float"),("mag","float"),("magerr","float"),("filter", "U1")]

        table_data_rev = np.array(table_data_rev, dtype=dt)
        table_data_V = table_data_rev[table_data_rev['filter']=='V']
        table_data_B = table_data_rev[table_data_rev['filter']=='B']
        table_data_Ic = table_data_rev[table_data_rev['filter']=='Ic']

        #remove spaces in output before creating file
        header = 'Time (JD) '+header[2]+' '+header[3]+' '+header[4]
        if not name:
            name = object

        name = name.replace(" ", "_")
        fmt = '%f %.3f %.3f %s'

        if len(table_data_V) > 0:
            print('V data was retrieved')
            np.savetxt(name+'_kws_V.dat', table_data_V, header=header, fmt=fmt)
        if len(table_data_B) > 0:
            print('B data was retrieved')
            np.savetxt(name+'_kws_B.dat', table_data_B, header=header, fmt=fmt)
        if len(table_data_V) > 0:
            print('Ic data was retrieved')
            np.savetxt(name+'_kws_Ic.dat', table_data_Ic, header=header, fmt=fmt)

    else:
        print("Sorry no data was found for {} in the KWS survey".format(name))

    return

"""
Name: asas_grab
Description: Download data for a given object from ASAS (All sky automated survey)
Inputs:
object: name of object. It is resolved using Simbad.
coords: list [RA, Dec] in decimal degrees
name: Object name specified in the saved filename. If not specified, name
will default to object name. If neither is specified, the ASAS ID will be used.

"""

def asas_grab(object=None, coords=None, name=None, asas_id=None):

    if not asas_id:
        #find coordinates
        if object:

            try:
                result_table = Simbad.query_object(object)
                RA = result_table['RA'][0]
                Dec = result_table['DEC'][0]
                ra_split = RA.split(' ')
                dec_split = Dec.split(' ')
                ra_h = ra_split[0]; ra_m = ra_split[1]; ra_s = ra_split[2]
                dec_d = dec_split[0]; dec_m = dec_split[1]; dec_s = dec_split[2]

            except:
                """ searching VSX for target """
                print("Sorry target {} not found in Simbad, searching vsx".format(object))
                #query vizier and find vsx catalog
                result = Vizier.query_object(object)
                indices = [i for i, s in enumerate(result.keys()) if 'vsx' in s]
                if len(indices) == 1:
                    vsx_table = result[indices[0]]
                    vsx_table = vsx_table[vsx_table['Name'] == object]
                #find RA and dec and convert to hms/dms
                    RA = vsx_table['RAJ2000'][0]
                    Dec = vsx_table['DEJ2000'][0]
                    c = SkyCoord(ra=RA*u.degree, dec=Dec*u.degree, frame='icrs')
                    ra_h = int(c.ra.hms[0]); ra_m = c.ra.hms[1]; ra_s = c.ra.hms[2]
                    dec_d = int(c.dec.dms[0]); dec_m = abs(c.dec.dms[1]); dec_s = abs(c.dec.dms[2])
                    if int(dec_d) > 0:
                        dec_d = '+'+str(dec_d)
                    else:
                        dec_d = str(dec_d)

                else:
                    print("Found too many tables in Vizier that match query.")
        elif coords:
            RA = coords[0] ; Dec = coords[1]
            c = SkyCoord(ra=RA*u.degree, dec=Dec*u.degree, frame='icrs')
            ra_h = int(c.ra.hms[0]); ra_m = c.ra.hms[1]; ra_s = c.ra.hms[2]
            dec_d = int(c.dec.dms[0]); dec_m = abs(c.dec.dms[1]); dec_s = abs(c.dec.dms[2])
            if int(dec_d) > 0:
                dec_d = '+'+str(dec_d)
            else:
                dec_d = str(dec_d)
        else:
            print("You must specify either the object or it's coordinates")
            return

        """
        #Configuring ASAS ID
        A = 'RA hours'
        B = 'RA minutes'
        C = 'RA seconds roundest to the nearest second'
        D = 'declination Degrees only'
        E = 'Degree minutes including seconds rounded to 1 decimal place'
        ASAS ID = A+B+C+D+E
        """

        #zfill pads with zeros if number is less than 10
        #two in most cases but 4 for dms because the period and the
        #number after the decimal count
        print(ra_h, ra_m, ra_s)
        asas_ra_id = str(int(ra_h)).zfill(2)+str(int(ra_m)).zfill(2)\
        +str(round(float(ra_s))).zfill(2)
    #    print(dec_d, dec_m, dec_s)
        dms = round(float(dec_m)+float(dec_s)/60., 1)


    #    print(dms, str(dms).zfill(4))
        asas_dec_id = dec_d.zfill(2)+str(dms).zfill(4)
        asas_id = asas_ra_id+asas_dec_id

    print('asas id is', asas_id)
    # Build URL

    """
    url consists base+asas_id+end where:
    base = 'www.astrouw.edu.pl/cgi-asas/asas_cgi_get_data?'
    asas_id = asas_id
    tail = ,asas3

    """
    url = 'http://www.astrouw.edu.pl/cgi-asas/asas_cgi_get_data'
    tail = ',asas3'
    end = '?'+asas_id+tail
    #end = end.encode('ascii')
    full_url = url+end

    #Get Data

    with urllib.request.urlopen(full_url) as response:
        the_page = response.read()

    the_page = the_page.decode('utf-8') #decode bytes
    lines = the_page.count('\n')

    #Check if data is available

    if lines < 54:

        # ASAS is often down so check
        # against survey parameters to see if data
        # should be available
        dec_max = 28.0
        dec_star = float(dec_d)
        if not name:
            if object:
                name = object
            else:
                name = "these coordinates "+str(RA)+", "+str(Dec)

        no_data_message = """No Data was found for {}, but it is within
        the surveys coverage area. This often means that the ASAS server
        is having problems. Please try again at a later date. """.format(name)

        if dec_star <= dec_max:
            print(' '.join(no_data_message.split()))

        else:

            print("No data was found for {}".format(name))

    else:
        if not name:
            if object:
                name = object
            else:
                name = asas_id
        print("Retrieved data for star {}".format(name))
        name = name.replace(" ", "_")
        f = open(name+'_asas.dat', 'w')
        f.write(the_page)
        f.close()

    return
