#!/usr/bin/env python



import os
import sys
import os.path as op
import mtpy.utils.filehandling as MTfh
import mtpy.core.edi as EDI
import fnmatch

def main():

    if len(sys.argv) < 2:
        print 'usage:  kml_from_edis.py <edi_folder> [<output filename>]'
        return

    edifolder = sys.argv[1]
    edilist = []
    try:
        if not op.isdir(edifolder):
            raise
        edilist = fnmatch.filter(os.listdir(edifolder),'*.[Ee][Dd][Ii]')
        edilist = [op.abspath(op.join(edifolder,i)) for i in edilist]
        if len(edilist) == 0:
            raise
    except:
        print 'No EDI files in folder {0}'.format(edifolder)
        return

    out_fn = None

    if len(sys.argv) > 2:
        out_fn = sys.argv[2]
        try:
            out_fn = op.abspath(op.join('.',out_fn))
            if not out_fn.lower().endswith('.kml'):
                out_fn += '.kml'
        except:
            print 'Cannot write to file {0}, using generic filename'.format(out_fn)
            out_fn = None


    outfilename = convert_edi_coordinates_to_kml_file(edilist, out_fn)
   
    print 'written to file {0}'.format(outfilename) 


def convert_edi_coordinates_to_kml_file(edi_filelist, outfilename = None):


    kml_string = ''
    kml_string += '<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n'
    kml_string += '<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n'
    kml_string += '<Document>\n'

    for edi_idx, edi in enumerate(edi_filelist):
        
        e = EDI.Edi()
        e.readfile(edi)
        lat = e.lat
        lon = e.lon
        ele = e.elev
        station = e.head['dataid']

        kml = []
       
        description = 'File: {0}'.format(e.filename)
    
        kml.append('  <Placemark>')
        kml.append('    <name>%s</name>' % station)
        kml.append('    <description>')
        kml.append('        <p>%s</p>' % description)
        kml.append('      </description>')
        kml.append('    <Point>')
        kml.append('      <coordinates>%f,%f,%f</coordinates>' % (lon, lat,ele))
        kml.append('    </Point>')
        kml.append('  </Placemark>')
    
        kml_string += '\n'.join(kml)

    kml_string += '\n</Document>\n'
    kml_string += '</kml>\n'




    if outfilename is None:
        outfilename = op.abspath('edi_coordinates.kml')
    else:
        try: 
            outfilename = op.abspath(op.join('.',outfilename))
        except:
            outfilename = op.abspath('edi_coordinates.kml')


    outfilename =   MTfh.make_unique_filename(outfilename)

    fileObj = open(outfilename, 'w' )
    fileObj.write(kml_string)
    fileObj.close()


    return outfilename




if __name__=='__main__':
    main()
