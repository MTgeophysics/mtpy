# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 11:54:38 2011

@author: a1185872
"""
#==============================================================================
import numpy as np
import scipy as sp
import os
import subprocess
import shutil
import fnmatch
from operator import itemgetter
import time
import matplotlib.colorbar as mcb
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
import mtpy.core.edi as mtedi
import mtpy.modeling.winglinktools as wlt
import matplotlib.pyplot as plt
import mtpy.utils.latlongutmconversion as utm2ll
#==============================================================================

occamdict = {'1':'resxy','2':'phasexy','3':'realtip','4':'imagtip','5':'resyx',
             '6':'phaseyx'}

def getdatetime():

    return time.asctime(time.gmtime())


def makestartfiles(parameter_dict):
    """
    create a startup file from a parameter dictionary with keys:
    
    ======================= ===================================================
    keys                      description
    ======================= ===================================================
    n_sideblockelements     number of padding cells in the horizontal direction
    n_bottomlayerelements   number of padding cells in the vertical direction
    itform                  iteration form
    description             desctription of inversion
    datetime                date and time of inversion run
    iruf                    roughness type
    idebug                  debug level
    nit                     first iteration number
    pmu                     first mu value
    rlast                   initial roughness value
    tobt                    starting misfit value
    ifftol                  misfit reached
    ================ ==========================================================
    """

    read_datafile(parameter_dict)

    parameter_dict['n_sideblockelements'] = 7
    parameter_dict['n_bottomlayerelements'] = 4

    parameter_dict['itform'] = 'not specified'
    parameter_dict['description'] = 'N/A'

    parameter_dict['datetime'] = getdatetime()
    

    parameter_dict['iruf'] = 1
    parameter_dict['idebug'] = 1
    parameter_dict['nit'] = 0
    parameter_dict['pmu'] = 5.0
    parameter_dict['rlast'] = 1.0E+07
    parameter_dict['tobt'] = 100.
    parameter_dict['ifftol'] = 0
    
    blocks_elements_setup(parameter_dict)
    
    get_model_setup(parameter_dict)
    
    writemeshfile(parameter_dict)
    writemodelfile(parameter_dict)
    writestartupfile(parameter_dict)

    MeshF = parameter_dict['meshfn']
    ModF  = parameter_dict['inmodelfn']
    SF    = parameter_dict['startupfn']
    
    return (MeshF,ModF,SF)

def writemeshfile(parameter_dict):
    """
    create a startup file from a parameter dictionary with keys:
    
    ======================= ===================================================
    keys                      description
    ======================= ===================================================
    mesh_positions_vert     
    ================ ==========================================================
    """
    mesh_positions_vert = parameter_dict['mesh_positions_vert']
    mesh_positions_hor  = parameter_dict['mesh_positions_hor']
    n_nodes_hor         = parameter_dict['n_nodes_hor']
    n_nodes_vert        = parameter_dict['n_nodes_vert']
    
    fh_mesh = file(parameter_dict['meshfn'],'w')
    mesh_outstring =''

    temptext = "MESH FILE FROM MTpy\n"
    mesh_outstring += temptext

    temptext = "%i %i %i %i %i %i\n"%(0,n_nodes_hor,n_nodes_vert,0,0,2)
    mesh_outstring += temptext

    temptext = ""
    for i in range(n_nodes_hor-1):
        temptext += "%.1f "%(mesh_positions_hor[i])
    temptext +="\n"
    mesh_outstring += temptext

    temptext = ""
    for i in range(n_nodes_vert-1):
        temptext += "%.1f "%(mesh_positions_vert[i])
    temptext +="\n"
    mesh_outstring += temptext

    mesh_outstring +="%i\n"%(0)

    for j in range(4*(n_nodes_vert-1)):
        tempstring=''
        tempstring += (n_nodes_hor-1)*"?"
        tempstring += '\n'
        mesh_outstring += tempstring
    

    fh_mesh.write(mesh_outstring)
    fh_mesh.close()



def writemodelfile(parameter_dict):
    "needed : filename,binding_offset,startcolumn, n_layers,layer_thickness,block_width"

    modelblockstrings = parameter_dict['modelblockstrings']
    nfev              = parameter_dict['nfev']
    lo_colnumbers     = parameter_dict['lo_colnumbers']
    boffset           = float(parameter_dict['binding_offset'])
    n_layers          = int(float(parameter_dict['n_layers']))

    
    fh_model = file(parameter_dict['inmodelfn'],'w')
    model_outstring =''

    temptext = "Format:           %s\n"%("OCCAM2MTMOD_1.0")
    model_outstring += temptext
    temptext = "Model Name:       %s\n"%(parameter_dict['modelname'])
    model_outstring += temptext
    temptext = "Description:      %s\n"%("Random Text")
    model_outstring += temptext
    temptext = "Mesh File:        %s\n"%(os.path.basename(parameter_dict['meshfn']))
    model_outstring += temptext
    temptext = "Mesh Type:        %s\n"%("PW2D")
    model_outstring += temptext
    temptext = "Statics File:     %s\n"%("none")
    model_outstring += temptext
    temptext = "Prejudice File:   %s\n"%("none")
    model_outstring += temptext
    temptext = "Binding Offset:   %.1f\n"%(boffset)
    model_outstring += temptext
    temptext = "Num Layers:       %i\n"%(n_layers)
    model_outstring += temptext

    for k in range(n_layers):
        n_meshlayers  = nfev[k]
        n_meshcolumns = lo_colnumbers[k]
        temptext="%i %i\n"%(n_meshlayers, n_meshcolumns)
        model_outstring += temptext

        temptext = modelblockstrings[k]
        model_outstring += temptext
        #model_outstring += "\n"
        

    temptext = "Number Exceptions:%i\n"%(0)
    model_outstring += temptext
    

    fh_model.write(model_outstring)
    fh_model.close()



def writestartupfile(parameter_dict):


    fh_startup = file(parameter_dict['startupfn'],'w')
    startup_outstring =''

    temptext = "Format:           %s\n"%(parameter_dict['itform'])
    startup_outstring += temptext
    temptext = "Description:      %s\n"%(parameter_dict['description'])
    startup_outstring += temptext
    temptext = "Model File:       %s\n"%(os.path.basename(parameter_dict['inmodelfn']))
    startup_outstring += temptext
    temptext = "Data File:        %s\n"%(os.path.basename(parameter_dict['datafile']))
    startup_outstring += temptext
    temptext = "Date/Time:        %s\n"%(parameter_dict['datetime'])
    startup_outstring += temptext
    temptext = "Max Iter:         %i\n"%(int(float(parameter_dict['n_max_iterations'])))
    startup_outstring += temptext
    temptext = "Req Tol:          %.1g\n"%(float(parameter_dict['targetrms']))
    startup_outstring += temptext
    temptext = "IRUF:             %s\n"%(parameter_dict['iruf'])
    startup_outstring += temptext
    temptext = "Debug Level:      %s\n"%(parameter_dict['idebug'])
    startup_outstring += temptext
    temptext = "Iteration:        %i\n"%(int(float(parameter_dict['n_max_iterations'])))
    startup_outstring += temptext
    temptext = "PMU:              %s\n"%(parameter_dict['pmu'])
    startup_outstring += temptext
    temptext = "Rlast:            %s\n"%(parameter_dict['rlast'])
    startup_outstring += temptext
    temptext = "Tlast:            %s\n"%(parameter_dict['tobt'])
    startup_outstring += temptext
    temptext = "IffTol:           %s\n"%(parameter_dict['ifftol'])
    startup_outstring += temptext
    temptext = "No. Parms:        %i\n"%(int(float(parameter_dict['n_parameters'])))
    startup_outstring += temptext
    temptext = ""
    for l in range(int(float(parameter_dict['n_parameters']))):
        temptext += "%.1g  "%(2.0)
    temptext += "\n"
    startup_outstring += temptext
     

    fh_startup.write(startup_outstring)
    fh_startup.close()

def read_datafile(parameter_dict):


    df = parameter_dict['datafile']
    F  = file(df,'r')
    datafile_content = F.readlines()
    F.close()

    #RELYING ON A CONSTANT FORMAT, ACCESSING THE PARTS BY COUNTING OF LINES!!!:
    
    n_sites = int(datafile_content[2].strip().split()[1])
    sitenames = []
    for i in range(n_sites):
        sitenames.append(datafile_content[3+i].strip())

    sitelocations=[]
    for i in range(n_sites):
        idx = 4+n_sites+i
        sitelocations.append(float(datafile_content[idx].strip()))

    n_freqs = int(datafile_content[2*n_sites+4].strip().split()[1])
    freqs=[]
    for i in range(n_freqs):
        idx = 2*n_sites+5+i
        freqs.append(float(datafile_content[idx].strip()))


    n_data = int(datafile_content[2*n_sites+5+n_freqs].strip().split()[2])
    

    parameter_dict['lo_site_names']     = sitenames
    parameter_dict['lo_site_locations'] = sitelocations
    parameter_dict['n_sites']           = n_sites
    parameter_dict['n_datapoints']      = n_data
    parameter_dict['n_freqs']           = n_freqs
    parameter_dict['lo_freqs']          = freqs
    
    

def get_model_setup(parameter_dict):
    """
    get model setup from a parameter dictionary:
    
    ======================= ===================================================
    keys                      description
    ======================= ===================================================
    ncol0                   number of columns
    n_layer                 number of desired layers
    nfe                     
    ================ ==========================================================
    """

    ncol0      = int(float(parameter_dict['ncol0']))
    n_layer    = int(float(parameter_dict['n_layers']))
    nfe        = parameter_dict['nfe']
    thickness  = parameter_dict['thickness']
    width      = parameter_dict['width']
    trigger    =  float(parameter_dict['trigger'])
    dlz = parameter_dict['dlz']

    modelblockstrings = []
    lo_colnumbers     = []

    ncol     = ncol0
    num_params = 0


    
    for layer_idx in range(n_layer):
        block_idx = 1

        #print layer_idx,len(thickness),len(width),ncol,len(dlz)
    
        while block_idx+2 < ncol-1 :

            #PROBLEM : 'thickness' has only "n_layer'-1 entries!!
            if not dlz[layer_idx] > (trigger*(width[block_idx]+width[block_idx+1])):
                block_idx += 1
                continue

            else:
                width[block_idx] += width[block_idx+1]
                nfe[block_idx]   += nfe[block_idx+1]

                for m in range(block_idx+2,ncol):
                    width[m-1] = width[m]
                    nfe[m-1]   = nfe[m]

                ncol -=1

        lo_colnumbers.append(ncol)

        tempstring = ""
        for j in range(ncol):
            tempstring += "%i "%(nfe[j])
        tempstring += "\n"
        modelblockstrings.append(tempstring)

        num_params += ncol

        #completely unnecessary!!! :
        if layer_idx == 0:
            mcol = ncol
        

    parameter_dict['modelblockstrings'] = modelblockstrings
    parameter_dict['lo_colnumbers']     = lo_colnumbers
    parameter_dict['n_parameters']      = np
    parameter_dict['n_cols_max']        = mcol
    

def blocks_elements_setup(parameter_dict):
    

    lo_sites = parameter_dict['lo_site_locations']
    n_sites  = len(lo_sites)
    maxwidth = float(parameter_dict['max_blockwidth'])

    nbot     = int(float(parameter_dict['n_bottomlayerelements']))
    nside    = int(float(parameter_dict['n_sideblockelements']))

    # j: index for finite elements
    # k: index for regularisation bricks
    # Python style: start with 0 instead of 1

    sitlok      = []
    sides       = []
    width       = []
    dly         = []
    nfe         = []
    thickness   = []
    nfev        = []
    dlz         = []
    bot         = []

    
    j = 0
    sitlok.append(lo_sites[0])

    for idx in range(1,n_sites-1):
        
        spacing           = lo_sites[idx] - lo_sites[idx-1]
        n_localextrasites = int(spacing/maxwidth) + 1

        for idx2 in range(n_localextrasites):
            sitlok.append(lo_sites[idx-1] + (idx2+1.)/float(n_localextrasites)*spacing )
            j += 1

    # nrk: number of total dummy stations
    nrk = j
    print "%i dummy stations defined"%(nrk)

    
    spacing1 = (sitlok[1]-sitlok[0])/2.
    sides.append(3*spacing1)


    for idx in range(1,nside):
        curr_side = 3*sides[idx-1]
        if curr_side > 1000000.:
            curr_side = 1000000.
        sides.append(curr_side)
        
    #-------------------------------------------

    j = 0
    k = 0

    firstblockwidth = 0.
    
    for idx in range(nside-1,-1,-1):
        firstblockwidth += sides[idx]
        dly.append(sides[idx])
        j += 1
        
    width.append(firstblockwidth)
    nfe.append(nside)
    
    dly.append(spacing1)
    dly.append(spacing1)
    j += 2
    nfe.append(2)
    width.append(2*spacing1)

    block_offset = width[1]

    k += 1

    dly.append(spacing1)
    dly.append(spacing1)
    j += 2
    nfe.append(2)
    width.append(2*spacing1)
    
    block_offset += spacing1

    k += 1

    #------------------------
    
    for idx in range(1,nrk-1):
        spacing2 = (sitlok[idx+1]-sitlok[idx])/2.
        dly.append(spacing1)
        dly.append(spacing2)
        j += 2
        nfe.append(2)
        width.append(spacing1+spacing2)
        k += 1
        spacing1 = spacing2

    dly.append(spacing2)
    dly.append(spacing2)

    j += 2
    nfe.append(2)
    width.append(2*spacing2)
    k += 1

    dly.append(spacing2)
    dly.append(spacing2)
    
    j += 2
    nfe.append(2)
    width.append(2*spacing2)
    k += 1

    width[-1] = 0.
    sides[0] = 3*spacing2

    #------------------------------
    
    for idx in range(1,nside):
        curr_side = 3*sides[idx-1]
        if curr_side > 1000000.:
            curr_side = 1000000.
        sides[idx] = curr_side


    lastblockwidth= 0.
    for idx in range(nside):
        j += 1
        lastblockwidth += sides[idx]
        dly.append(sides[idx])

    width[-1] = lastblockwidth

    #---------------------------------

    k+= 1
    nfe.append(nside)

    nodey = j+1
    ncol0 = k

    block_offset = sitlok[0] - block_offset

    #----------------------------------

    layers_per_decade     = float(parameter_dict['n_layersperdecade'])
    first_layer_thickness = float(parameter_dict['firstlayer_thickness'])

    t          = 10.**(1./layers_per_decade)
    t1         = first_layer_thickness
    thickness.append(t1)
    
    d1 = t1
    
    n_layers = int(float(parameter_dict['n_layers']))

    for idx in range(1,n_layers-1):
        d2 = d1*t
        curr_thickness = d2 - d1
        if curr_thickness < t1:
            curr_thickness = t1
        thickness.append(curr_thickness)
        d1 += curr_thickness
    
    
    bot.append(3*thickness[n_layers-2])

    for idx in range(1,nbot):
        bot.append(bot[idx-1]*3)

    #--------------------------------------------------

    k = 0
    
    dlz.append(thickness[0]/2.)
    dlz.append(thickness[0]/2.)
    nfev.append(2)

    k += 2

    dlz.append(thickness[1]/2.)
    dlz.append(thickness[1]/2.)
    nfev.append(2)

    k += 2

    for idx in range(2,n_layers-1):
        k += 1
        nfev.append(1.)
        dlz.append(thickness[idx])

    for idx in range(nbot):
        k += 1
        dlz.append(bot[idx])

    nfev.append(nbot)

    nodez = k+1
    

    parameter_dict['ncol0']             = ncol0
    parameter_dict['nfe']               = nfe
    parameter_dict['nfev']              = nfev
    parameter_dict['thickness']         = thickness
    parameter_dict['width']             = width
    parameter_dict['binding_offset']    = block_offset
    #parameter_dict['y_nodes']           = nodey 
    #parameter_dict['z_nodes']           = nodez
    parameter_dict['dlz']               = dlz
    #parameter_dict['dly']               = dly
    parameter_dict['mesh_positions_vert'] = dlz
    parameter_dict['mesh_positions_hor']  = dly
    parameter_dict['n_nodes_hor']         = nodey
    parameter_dict['n_nodes_vert']        = nodez

class OccamPointPicker(object):
    """
    This class helps the user interactively pick points to mask and add 
    error bars. 
    
    Useage:
    -------
    To mask just a single point right click over the point and a gray point 
    will appear indicating it has been masked
    
    To mask both the apparent resistivity and phase left click over the point.
    Gray points will appear over both the apparent resistivity and phase.  
    Sometimes the points don't exactly matchup, haven't quite worked that bug
    out yet, but not to worry it picks out the correct points
    
    To add error bars to a point click the middle or scroll bar button.  This
    only adds error bars to the point and does not reduce them so start out
    with reasonable errorbars.  You can change the increment that the error
    bars are increased with reserrinc and phaseerrinc.
    
    Arguments:
    ----------
        **axlst** : list of the resistivity and phase axis that have been 
                    plotted as [axr_te,axr_tm,axp_te,axp_tm]
        
        **linelst** : list of lines used to plot the responses, not the error 
                      bars as [res_te,res_tm,phase_te,phase_tm]
        
        **errlst** : list of the errorcaps and errorbar lines as 
                   [[cap1,cap2,bar],...]
                 
        **reserrinc** : increment to increase the errorbars for resistivity.
                        put .20 for 20 percent change. *Default* is .05
        
        **phaseerrinc** : increment to increase the errorbars for the phase
                          put .10 for 10 percent change. *Defualt* is .02 
                    
        **marker** : marker type for masked points.  See matplotlib.pyplot.plot
                    for options of markers.  *Default* is h for hexagon.
                   
    Attributes:
    -----------
    
        **axlst** : axes list used to plot the data
        
        **linelst** : line list used to plot the data
        
        **errlst** : error list used to plot the data
        
        **data** : list of data points that were not masked for each plot.
        
        **fdict** : dictionary of frequency arrays for each plot and data set.
        
        **fndict** : dictionary of figure numbers to corresponed with data.
        
        **cidlst** : list of event ids.
        
        **reserrinc** : increment to increase resistivity error bars
        
        **phaseinc** : increment to increase phase error bars
        
        **marker** : marker of masked points
        
        **fignum** : figure numbers
        
        **occamlines** : list of lines to write into the occam data file.
        
    :Example: ::
        
        >>> ocd = occam.Occam2DData()
        >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/Data.dat"
        >>> ocd.plotMaskPoints()
    """    
    
    def __init__(self,axlst,linelst,errlst,reserrinc=.05,phaseerrinc=.02,
                 marker='h'):
          
        #give the class some attributes
        self.axlst=axlst
        self.linelst=linelst
        self.errlst=errlst
        self.data=[]
        self.error=[]
        self.fdict=[]
        self.fndict={}
        #see if just one figure is plotted or multiple figures are plotted
        self.ax=axlst[0][0]
        self.line=linelst[0][0]
        self.cidlst=[]
        for nn in range(len(axlst)):
            self.data.append([])
            self.error.append([])
            self.fdict.append([])
        
            #get data from lines and make a dictionary of frequency points for 
            #easy indexing
            for ii,line in enumerate(linelst[nn]):
                self.data[nn].append(line.get_data()[1])
                self.fdict[nn].append(dict([('{0:.5g}'.format(kk),ff) for ff,kk in 
                                        enumerate(line.get_data()[0])]))
                self.fndict['{0}'.format(line.figure.number)]=nn
                
                #set some events
                if ii==0:
                    cid1=line.figure.canvas.mpl_connect('pick_event',self)
                    cid2=line.figure.canvas.mpl_connect('axes_enter_event',
                                                        self.inAxes)
                    cid3=line.figure.canvas.mpl_connect('key_press_event',
                                                        self.on_close)
                    cid4=line.figure.canvas.mpl_connect('figure_enter_event',
                                                        self.inFigure)
                    self.cidlst.append([cid1,cid2,cid3,cid4])
        
            #read in the error in a useful way so that it can be translated to 
            #the data file.  Make the error into an array
            for ee,err in enumerate(errlst[nn]):
                errpath=err[2].get_paths()
                errarr=np.zeros(len(self.fdict[nn][ee].keys()))
                for ff,epath in enumerate(errpath):
                    errv=epath.vertices
                    errarr[ff]=abs(errv[0,1]-self.data[nn][ee][ff])
                self.error[nn].append(errarr)
        
        #set the error bar increment values
        self.reserrinc=reserrinc
        self.phaseerrinc=phaseerrinc
        
        #set the marker
        self.marker=marker
        
        #set the figure number
        self.fignum=self.line.figure.number
        
        #make a list of occam lines to write later
        self.occamlines=[]
    
        
        
    
    def __call__(self,event):
        """
        When the function is called the mouse events will be recorder for 
        picking points to mask or change error bars.  The axes is redrawn with
        a gray marker to indicate a masked point and/or increased size in 
        errorbars.
        
        Arguments:
        ----------
            **event** : type mouse_click_event
                
        Useage:
        -------
        
            **Left mouse button** will mask both resistivity and phase point
        
            **Right mouse button** will mask just the point selected
        
            **Middle mouse button** will increase the error bars
        
            **q** will close the figure.
        """
        self.event=event
        #make a new point that is an PickEvent type
        npoint=event.artist
        #if the right button is clicked mask the point
        if event.mouseevent.button==3:
            #get the point that was clicked on
            ii=event.ind
            xd=npoint.get_xdata()[ii]
            yd=npoint.get_ydata()[ii]
            
            #set the x index from the frequency dictionary
            ll=self.fdict[self.fignum][self.jj]['{0:.5g}'.format(xd[0])]
            
            #change the data to be a zero
            self.data[self.fignum][self.jj][ll]=0
            
            #reset the point to be a gray x
            self.ax.plot(xd,yd,ls='None',color=(.7,.7,.7),marker=self.marker,
                         ms=4)
        
        #if the left button is clicked change both resistivity and phase points
        elif event.mouseevent.button==1:
            #get the point that was clicked on
            ii=event.ind
            xd=npoint.get_xdata()[ii]
            yd=npoint.get_ydata()[ii]
            
            #set the x index from the frequency dictionary
            ll=self.fdict[self.fignum][self.jj]['{0:.5g}'.format(xd[0])]
            
            #set the data point to zero
            self.data[self.fignum][self.jj][ll]=0
            
            #reset the point to be a gray x
            self.ax.plot(xd,yd,ls='None',color=(.7,.7,.7),marker=self.marker,
                         ms=4)
            
            #check to make sure there is a corresponding res/phase point
            try:
                #get the corresponding y-value 
                yd2=self.data[self.fignum][self.kk][ll]
                
                #set that data point to 0 as well
                self.data[self.fignum][self.kk][ll]=0
                
                #make that data point a gray x
                self.axlst[self.fignum][self.kk].plot(xd,yd2,ls='None',
                                        color=(.7,.7,.7),marker=self.marker,
                                        ms=4)
            except KeyError:
                print 'Axis does not contain res/phase point'
                
        #if click the scroll button or middle button change increase the 
        #errorbars by the given amount
        elif event.mouseevent.button==2:
            ii=event.ind
            xd=npoint.get_xdata()[ii]
            yd=npoint.get_ydata()[ii]
            
            #get x index
            ll=self.fdict[self.fignum][self.jj]['{0:.5g}'.format(xd[0])]
            
            #make error bar array
            eb=self.errlst[self.fignum][self.jj][2].get_paths()[ll].vertices
            
            #make ecap array
            ecapl=self.errlst[self.fignum][self.jj][0].get_data()[1][ll]
            ecapu=self.errlst[self.fignum][self.jj][1].get_data()[1][ll]
            
            #change apparent resistivity error
            if self.jj==0 or self.jj==1:
                nebu=eb[0,1]-self.reserrinc*eb[0,1]
                nebl=eb[1,1]+self.reserrinc*eb[1,1]
                ecapl=ecapl-self.reserrinc*ecapl
                ecapu=ecapu+self.reserrinc*ecapu
                
            #change phase error
            elif self.jj==2 or self.jj==3:
                nebu=eb[0,1]-eb[0,1]*self.phaseerrinc
                nebl=eb[1,1]+eb[1,1]*self.phaseerrinc
                ecapl=ecapl-ecapl*self.phaseerrinc
                ecapu=ecapu+ecapu*self.phaseerrinc
                
            #put the new error into the error array    
            self.error[self.fignum][self.jj][ll]=abs(nebu-\
                                        self.data[self.fignum][self.jj][ll])
            
            #set the new error bar values
            eb[0,1]=nebu
            eb[1,1]=nebl
            
            #reset the error bars and caps
            ncapl=self.errlst[self.fignum][self.jj][0].get_data()
            ncapu=self.errlst[self.fignum][self.jj][1].get_data()
            ncapl[1][ll]=ecapl
            ncapu[1][ll]=ecapu
            
            #set the values 
            self.errlst[self.fignum][self.jj][0].set_data(ncapl)
            self.errlst[self.fignum][self.jj][1].set_data(ncapu)
            self.errlst[self.fignum][self.jj][2].get_paths()[ll].vertices=eb
            
        #redraw the canvas
        self.ax.figure.canvas.draw()

    #get the axis number that the mouse is in and change to that axis
    def inAxes(self,event):
        """
        gets the axes that the mouse is currently in.
        
        Arguments:
        ---------
            **event**: is a type axes_enter_event
                
        Returns:
        --------
        
            **OccamPointPicker.jj** : index of resistivity axes for axlst
            
            **OccamPointPicker.kk** : index of phase axes for axlst
        
        """
        
        self.event2=event
        self.ax=event.inaxes
        for jj,axj in enumerate(self.axlst):
            for ll,axl in enumerate(axj):
                if self.ax==axl:
                    self.jj=ll
                
        #set complimentary resistivity and phase plots together
        if self.jj==0:
            self.kk=2
        if self.jj==1:
            self.kk=3
        if self.jj==2:
            self.kk=0
        if self.jj==3:
            self.kk=1
        
    #get the figure number that the mouse is in
    def inFigure(self,event):
        """
        gets the figure number that the mouse is in
        
        Arguments:
        ----------
            **event** : figure_enter_event
            
        Returns:
        --------
            **OccamPointPicker.fignum** : figure number that corresponds to the
                                          index in the axlst, datalst, errorlst
                                          and linelst.
                        
        """
        self.event3=event
        self.fignum=self.fndict['{0}'.format(event.canvas.figure.number)]
        self.line=self.linelst[self.fignum][0]
    
    #type the q key to quit the figure and disconnect event handling            
    def on_close(self,event):
        """
        close the figure with a 'q' key event and disconnect the event ids
        
        Arguments:
        ----------
            **event** : key_press_event
               
        Returns:
        --------
            print statement saying the figure is closed
        """
        self.event3=event
        if self.event3.key=='q':
            for cid in self.cidlst[self.fignum]:
               event.canvas.mpl_disconnect(cid)
            plt.close(event.canvas.figure)
            print 'Closed figure ',self.fignum  

            
class Occam2DData(object):
    """
    Occam2DData covers all aspects of dealing with data for an Occam 2D
    inversion using the code of Constable et al. [1987] and deGroot-Hedlin and 
    Constable [1990] from Scripps avaliable at 
    http://marineemlab.ucsd.edu/Projects/Occam/2DMT/index.html.
    
    
    """    
    
    def __init__(self,data_fn=None):
        self.data_fn = data_fn
        self.survey_lst = None
        self.freq_lst = None
        self.thetar = 0
        self.edipath = None
        self.station_lst = None
        self.proj_angle = None
        self.line_orientation = 'ew'
        
        self.resxy_err = 10
        self.resyx_err = 10
        self.phasexy_err = 5
        self.phaseyx_err = 5
        self.tipper_err = None
        
        self.freq_dict = None
        self.freq_step = 1
        self.ftol = 0.05
        
        self.model_mode = 'both'
        self.save_path = None
        self.title = None
        self.thetar = 0
        self._ss = ' '*3
        self._string_fmt = ' 2.6f'


    def _get_station_data(self, edipath, station_lst=None, thetar=0):
        """
        get the important information from the edifiles to invert for
        
        """
        
        self.thetar = thetar
        if abs(self.thetar) > 2*np.pi:
            self.thetar = self.thetar*(np.pi/180)
            
        #create a list to put all the station dictionaries into
        self.survey_lst = []
        self.pstation_lst = []
        self.freq_lst = []
        
        self.edipath = edipath
        self.station_lst = station_lst
        
        #get edi files for all stations in edipath if station_lst is None
        if station_lst == None:
            self.station_lst = [edifile[:-4] 
                for edifile in os.listdir(edipath) if edifile.find('.edi')]
        
        for kk, station in enumerate(self.station_lst):
            #search for filenames in the given directory and match to station 
            #name
            for filename in os.listdir(edipath):
                if fnmatch.fnmatch(filename,station+'*.edi'):
                    print 'Found station edifile: ', filename
                   
                    #create a dictionary for the station data and info 
                    surveydict = {} 
                    edifile = os.path.join(edipath,filename) 
                    
                    #read in edi file
                    z1 = mtedi.Edi()
                    z1.readfile(edifile)
                    freq = z1.freq
                    
                    #rotate data
                    if self.thetar != 0:
                        z1.Z.rotate(self.thetar)
                        if z1.Tipper.tipper is not None:
                            z1.Tipper.rotate(self.thetar)
                    
                    #get resistivity and phase
                    res, phase, res_err, phase_err = z1.Z.res_phase
                    
                    #check to see if the frequency is in descending order
                    if freq[0]<freq[-1]:
                        freq = freq[::-1]
                        res = res[::-1]
                        res_err = res_err[::-1]
                        phase = phase[::-1]
                        phase_err = phase_err[::-1]
                        if z1.Tipper.tipper is not None:
                            tip = z1.Tipper.tipper[::-1,:,:]
                            tipvar = z1.Tipper.tipper_err[::-1,:]
                        
                        print ('Flipped frequency to descending for station: '
                               '{0}'.format(station))
                    else:
                        if z1.Tipper.tipper is not None:
                            tip = z1.Tipper.tipper
                            tipvar = z1.Tipper.tipper_err
                            
                    #get eastings and northings so everything is in meters
                    zone, east, north = utm2ll.LLtoUTM(23, z1.lat, z1.lon)
                    
                    #put things into a dictionary to sort out order of stations
                    surveydict['station'] = station
                    surveydict['east'] = east
                    surveydict['north'] = north
                    surveydict['zone'] = zone
                    
                    surveydict['resxy'] = res[:, 0, 1]
                    surveydict['resxy_err'] = res_err[:, 0, 1]
                    surveydict['resyx'] = res[:, 1, 0]
                    surveydict['resyx_err'] = res_err[:, 1, 0]
                    
                    surveydict['phasexy'] = phase[:, 0, 1]
                    surveydict['phasexy_err'] = phase_err[:, 0, 1]
                    surveydict['phaseyx'] = phase[:, 1, 0]
                    surveydict['phaseyx_err'] = phase_err[:, 1, 0]
                    surveydict['freq'] = freq
                    
                    if z1.Tipper.tipper is not None:
                        surveydict['tipper'] = tip
                        surveydict['tippervar'] = tipvar
                        
                    surveydict['lat'] = z1.lat
                    surveydict['lon'] = z1.lon
                    
                    self.freq_lst.append(freq)
                    self.pstation_lst.append(station)
                    self.survey_lst.append(surveydict)
        
    def _project_stations(self, proj_angle=None, plot_yn='y'):
        """
        project stations onto a line
        
        """
        
        #get information from suvey_lst
        east_lst = np.array([sdict['east'] for sdict in self.survey_lst])
        north_lst = np.array([sdict['north'] for sdict in self.survey_lst])
        
        #get bestfitting line
        p = sp.polyfit(east_lst, north_lst, 1)
        
        #proj_slope needs to be orthogonal to the projection angle and inverse
        #need to subtract 90 because measured from 0, which is assumed to be
        #north and positive clockwise, where the projection assumes angle
        #is positive counter-clockwise
        if proj_angle is not None:
            print 'projecting stations onto {0} deg line'.format(proj_angle)
            if proj_angle == 0:
                proj_angle = 0.0001
            self.proj_angle = proj_angle
            self.proj_slope = -1./np.arctan(np.deg2rad(proj_angle))
            p[0] = np.arctan(np.deg2rad(proj_angle))
        else:
            self.proj_angle = np.rad2deg(np.arctan(p[0]))
            self.proj_slope = -1/p[0]
        
        new_east_lst = np.zeros_like(east_lst)
        new_north_lst = np.zeros_like(north_lst)
        
        #point is projected by calculating the orthogonal line and matching
        #the two at a common point
        # y1 = m1 x +b1; y2 = m2 x + b2 --> m2 = -1/m1 
        #calculate b2 and match y1 == y2
        for east, north, ii in zip(east_lst, north_lst, 
                                   range(new_east_lst.shape[0])):
                
            new_b = north-self.proj_slope*east
            new_east = (new_b-p[1])/(p[0]-self.proj_slope)
            new_north = p[0]*new_east+p[1]
            new_east_lst[ii] = new_east
            new_north_lst[ii] = new_north
            self.survey_lst[ii]['east_proj'] = new_east
            self.survey_lst[ii]['north_proj'] = new_north
        
        #plot stations on profile line
        if plot_yn == 'y':
            lfig = plt.figure(4, dpi=200)
            plt.clf()
            lax = lself.fig.add_subplot(1, 1, 1,aspect='equal')
            
            #plot the line that stations have been projected onto
            ploty = sp.polyval(p, new_east_lst)
            lax.plot(new_east_lst, ploty, '-b', lw=2)

            #plot stations
            lax.scatter(new_east_lst, new_north_lst, marker='v', s=50, 
                        color='k')
            for sdict in self.survey_lst:
                lax.text(sdict['east_proj'], sdict['north_proj']+100,
                         sdict['station'], verticalalignment='baseline',
                        horizontalalignment='center', 
                        fontdict={'size':12, 'weight':'bold'})
                        
            lax.set_title('Projected Stations')
            lax.set_ylim(new_north_lst.min()-1000., 
                         new_north_lst.max()+1000.)
            lax.set_xlim(new_east_lst.min()-1000., 
                         new_east_lst.max()+1000.)
            lax.set_xlabel('Easting (m)', 
                           fontdict={'size':12, 'weight':'bold'})
            lax.set_ylabel('Northing (m)',
                           fontdict={'size':12, 'weight':'bold'})
            plt.show()

        
    def _get_station_offsets(self, line_orientation='ew'):
        """
        get relative offsets between stations, assuming stations have been 
        projected onto a line
        
        returns a sorted list in ascending offsets
        """
        
        self.line_orientation = line_orientation
        for ii in range(len(self.survey_lst)):
            if self.survey_lst[ii]['zone'] != self.survey_lst[0]['zone']:
                print 'Different zone for {0}'.format(
                                            self.survey_lst[ii]['station'])
            
            #need to figure out a way to account for zone changes
            
            if self.line_orientation == 'ew': 
                if self.survey_lst[0]['east_proj'] < \
                   self.survey_lst[ii]['east_proj']:
                    self.survey_lst[ii]['offset'] = \
                                    np.sqrt((self.survey_lst[0]['east_proj']-
                                    self.survey_lst[ii]['east_proj'])**2+
                                    (self.survey_lst[0]['north_proj']-
                                    self.survey_lst[ii]['north_proj'])**2)
                elif self.survey_lst[0]['east_proj'] > \
                     self.survey_lst[ii]['east_proj']:
                    self.survey_lst[ii]['offset'] = \
                                  -1*np.sqrt((self.survey_lst[0]['east_proj']-
                                  self.survey_lst[ii]['east_proj'])**2+
                                  (self.survey_lst[0]['north_proj']-
                                  self.survey_lst[ii]['north_proj'])**2)
                else:
                    self.survey_lst[ii]['offset'] = 0
                    
            elif self.line_orientation == 'ns': 
                if self.survey_lst[0]['north_proj'] < \
                   self.survey_lst[ii]['north_proj']:
                    self.survey_lst[ii]['offset'] = \
                                    np.sqrt((self.survey_lst[0]['east_proj']-
                                    self.survey_lst[ii]['east_proj'])**2+
                                    (self.survey_lst[0]['north_proj']-
                                    self.survey_lst[ii]['north_proj'])**2)
                elif self.survey_lst[0]['north_proj'] > \
                     self.survey_lst[ii]['north_proj']:
                    self.survey_lst[ii]['offset'] = \
                                -1*np.sqrt((self.survey_lst[0]['east_proj']-
                                self.survey_lst[ii]['east_proj'])**2+
                                (self.survey_lst[0]['north_proj']-
                                self.survey_lst[ii]['north_proj'])**2)
                else:
                    self.survey_lst[ii]['offset'] = 0
                    
        #sort by ascending order of distance from first station
        self.survey_lst = sorted(self.survey_lst,key=itemgetter('offset'))
             
    def make2DdataFile(self, edipath, **kwargs):
        """
        Make a data file that Occam can read.  At the moment the inversion line
        is the best fit line through all the stations used for the inversion.
        
        Arguments:
        ----------
            **edipath** : path to edifiles
            
            **model_mode** : modes to invert for.  Can be: 
                            * 'both' -> will model both TE and TM modes
                            * 'TM'   -> will model just TM mode
                            * 'TE'   -> will model just TE mode
                        
            **save_path** : path to save the data file to, this can include the 
                           name of the data file, if not the file will be 
                           named: save_path/Data.dat or edipath/Data.dat if 
                           save_path=None
                           
            **station_lst** : list of stations to put in the data file, doesn't 
                             need to be in order, the relative distance will be
                             calculated internally.  If station_lst:None, it
                             will be assumed all the files in edipath will be 
                             input into the data file
                             
            **title** : title input into the data file. *Default* is None
            
            **thetar** : rotation angle (deg) of the MT response to align 
                         TE mode along strike which is assumed to be 
                         perpendicular to the inversion line and TM
                         parallel to the inversion line.  Angle is on 
                         the unit circle with an orientation that north is 0 
                         degree, east 90. Rotations are postive clockwise.
                         *Default* is 0 (North).  
                         
                         If proj_strike = 'yes' then the stations will be 
                         projected onto a line that is perpendicular to 
                         the thetar, which is assumed to be the strike
                         direction. 
            
            **resxy_err** : percent error in the res_xy component (TE), 
                          can be entered as 'data' where the errors from the 
                          data are used otherwise enter as a percentage.
                          enter 10 for 10 percent. *Default* is 10
            
            **resyx_err** : percent error in the res_yx component (TM), 
                          can be entered as 'data' where the errors from the 
                          data are used otherwise enter as a percentage.
                          enter 10 for 10 percent. *Default* is 10
            
            **phasexy_err** : percent error in the phase_xy component (TE), 
                             can be entered as 'data' where the errors from the 
                             data are used otherwise enter as a percentage.
                             enter 10 for 10 percent. *Default* is 5  
            
            **phaseyx_err** : percent error in the phase_yx component (TM), 
                             can be entered as 'data' where the errors from the 
                             data are used otherwise enter as a percentage.
                             enter 10 for 10 percent. *Default* is 5  
            
            **ss** : string
                    is the spacing parameter for the data file. 
                     *Default* is ' '*3 (3 spaces) 
            
            **string_fmt** : format of the numbers for the data file, see string 
                      formats for a full description. *Default* is '%+2.6f
            
            **freq_step** : take frequencies at this step, so if you want to 
                           take every third frequency enter 3.  
                           Can input as a list of specific frequencies.  
                           Note that the frequencies must match the frequencies
                           in the EDI files within ftol, otherwise they will 
                           not be input.  
            
            **ftol** : tolerance level (decimal %) to match frequencies to 
                       freq_step if input as a list.  *Default* is .05
            
            **plotyn** : y or n to plot the stations that are pojected on to 
                         the best fitting line through the stations.
            
            **line_orientation** : predominant line orientation with respect to 
                          geographic north
                          ew for east-west line-> will orientate so first 
                                                  station is farthest to the 
                                                  west
                          ns for north-south line-> will orientate so first 
                                                    station is farthest to the
                                                    south
                                                    
            
            **proj_angle** : angle in degrees
                             stations will be projected onto a line with this
                             angle assuming North is 0 and east is 90.
                             If None, then the stations will be projected onto
                             a best fitting line through all the stations.
                             projection is orthogonal.
            
            **tipper_err** : error for tipper in percent.  If this value is 
                            entered than the tipper will be included in the 
                            inversion, if the value is None than the tipper 
                            will not be included. 
                            Can be entered as 'data' where the errors from the 
                            data are used otherwise enter as a percentage.
                            enter 10 for 10 percent.
                            
        Returns:
        --------
            **Occam2DData.data_fn** : full path of data file
            
        :Example: ::
            
            >>> import mtpy.core.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> #define a path to where the edifiles are
            >>> edipath = r"/home/EDIfiles"
            >>>
            >>> #create a save path to put the data file
            >>> svpath = r"/home/Occam2D/Line1/Inv1"
            >>>
            >>> #create a station list of stations to use in inversion
            >>> slst = ['MT0{0}'.format(ii) for ii in range(1,10)]
            >>> 
            >>> #write data file that is rotated 50 degrees east of north and
            >>> #projected on to the strike direction 50 degrees east of north
            >>> ocd.make2DdataFile(edipath,station_lst=slst,save_path=svpath,
            >>> ...                thetar=50,proj_strike='yes',
            >>> ...                line_orientation='ew')
            >>>                       
            >>> 'Wrote Occam2D data file to: /home/Occam2D/Line1/Inv1/Data.dat' 

                     
        """
        # --> get information from key word arguments
        self.model_mode = kwargs.pop('model_mode', self.model_mode)
        self.station_lst = kwargs.pop('station_lst', self.station_lst)
        self.title = kwargs.pop('title', self.title)
        self.thetar = kwargs.pop('thetar', self.thetar)
        self.resxy_err = kwargs.pop('resxy_err', self.resxy_err)
        self.resyx_err = kwargs.pop('resyx_err', self.resyx_err)
        self.phasexy_err = kwargs.pop('phasexy_err', self.phasexy_err)
        self.phasexy_err = kwargs.pop('phasexy_err', self.phaseyx_err)
        self.tipper_err = kwargs.pop('tipper_err', self.tipper_err)
        self.freq_step = kwargs.pop('freq_step', self.freq_step)
        plot_yn = kwargs.pop('plot_yn', 'y')
        self.proj_angle = kwargs.pop('proj_angle', self.proj_angle)
        self.line_orientation = kwargs.pop('line_orientation', 
                                       self.line_orientation)
        self.save_path = kwargs.pop('save_path', self.save_path)
        self._ss = kwargs.pop('ss', self._ss)
        self._string_fmt = kwargs.pop('string_fmt', self._string_fmt)
        
        #make local variables with shorter names for writing to a file
        ss = self._ss
        sfmt = self._string_fmt

        
        #--> get important information from the edi files and rotate data
        self._get_station_data(edipath, station_lst=self.station_lst, 
                               thetar=self.thetar)
        
        #-----------------------------------------------------------------            
        #project stations onto a best fitting line taking into account the 
        #strike direction to get relative MT distances correct
        #-----------------------------------------------------------------
        self._project_stations(proj_angle=self.proj_angle, plot_yn=plot_yn)
        
        #--> get relative offsets between stations in meters
        self._get_station_offsets(line_orientation=self.line_orientation)

        #number of stations read    
        nstat=len(self.survey_lst)    
        
        #--------------------------Match Frequencies---------------------------
        #a dictionary is created with the frequency as the key and the value is
        #the frequency number in the list. Each edi file is iterated over 
        #extracting only the matched frequencies.  This makes it necessary to 
        #have the same frequency content in each edifile.  If the frequencies
        #do not match then you can specify a tolerance to look around for 
        #each frequency.
        
        #make a list to iterate over frequencies
        if type(self.freq_step) is list or type(self.freq_step) is not int:
            if type(self.freq_step[0]) is int:
                #find the median frequency list
                maxflen = max([len(ff) for ff in self.freq_lst])
                farray = np.zeros((nstat,maxflen))
                for ii in range(nstat):
                    farray[ii,0:len(self.freq_lst[ii])] = self.freq_lst[ii]
            
                mfreq = np.median(farray,axis=0)
                self.freq_dict = dict([('%.6g' % mfreq[ff],ii) 
                               for ii,ff in enumerate(self.freq_step,1) 
                               if mfreq[ff]!=0])
            else:
                self.freq_dict=dict([('%.6g' % ff,ii) 
                            for ii,ff in enumerate(self.freq_step,1)])
        else:
            #find the median frequency list
            maxflen=max([len(ff) for ff in self.freq_lst])
            farray=np.zeros((nstat,maxflen))
            for ii in range(nstat):
                farray[ii,0:len(self.freq_lst[ii])]=self.freq_lst[ii]
            
            mfreq=np.median(farray,axis=0)
        
            #make a dictionary of values        
            self.freq_dict=dict([('%.6g' % ff,ii) for ii,ff in 
                        enumerate(mfreq[range(0,maxflen,self.freq_step)],1) 
                        if ff!=0])
    
        #print the frequencies to look for to make sure its what the user wants
        #make a list of keys that is sorted in descending order
        klst=[float(dd) for dd in self.freq_dict.keys()]
        klst.sort(reverse=True)
        klst=['%.6g' % dd for dd in klst]    
        
        print 'Frequencies to look for are: (# freq(Hz) Period(s)) '
        for key in klst:
            print self.freq_dict[key],key, 1./float(key)
        
        #make lists of parameters to write to file    
        reslst = []
        offsetlst = []
        station_lstsort = []
        for kk in range(nstat):
            #--> set local variable with shorter names
            sresxy = self.survey_lst[kk]['resxy']
            sresxy_err = self.survey_lst[kk]['resxy_err']
            sresyx = self.survey_lst[kk]['resyx']
            sresyx_err = self.survey_lst[kk]['resyx_err']
            
            sphasexy = self.survey_lst[kk]['phasexy']
            sphasexy_err = self.survey_lst[kk]['phasexy_err']
            sphaseyx = self.survey_lst[kk]['phaseyx']+180
            sphaseyx_err = self.survey_lst[kk]['phaseyx_err']
            
            freq = self.survey_lst[kk]['freq']
            offsetlst.append(self.survey_lst[kk]['offset'])  
            station_lstsort.append(self.survey_lst[kk]['station'])
            try:
                tip = self.survey_lst[kk]['tipper']
            except KeyError:
                pass
            
            #loop over frequencies to pick out the ones desired
            for jj, ff in enumerate(freq):
                #jj is the index of edi file frequency list, this index 
                #corresponds to the impedance tensor component index
                #ff is the frequency from the edi file frequency list
                nn = None
                try:
                    #nn is the frequency number out of extracted frequency list
                    nn = self.freq_dict['%.6g' % ff]
                except KeyError:
                    #search around the frequency given by ftol
                    try:
                        for key in self.freq_dict.keys():
                            if ff > float(key)*(1-self.ftol) and \
                               ff < float(key)*(1+self.ftol):
                                nn = self.freq_dict[key]                           
                    except KeyError:
                        pass
                    
                if nn is not None:
                    #calculate apparent resistivity 
                    resxy = sresxy[jj]
                    resyx = sresyx[jj]
            
                    #calculate the phase putting the yx in the 1st quadrant        
                    phasexy = '{0:{1}}'.format(sphasexy[jj], sfmt)
                    phaseyx = '{0:{1}}'.format(sphaseyx[jj], sfmt)
                    
                    #put phases in correct quadrant if should be negative
                    if float(phaseyx) > 180:
                        phaseyx = '{0:{1}}'.format(float(phaseyx)-360, sfmt)
                    if float(phaseyx) < 0:
                        print ('Negative Phase at {0}'.format(
                               self.survey_lst[kk]['station'])+
                               'f={0:.4g}Hz, phi_tm={1: .2f}'.format(ff, 
                                                               float(phaseyx)))    
                    
                    #calculate errors
                    #--> res_xy (TE)
                    if self.resxy_err == 'data':
                        dresxy_err = (sresxy_err[jj]/resxy)/np.log(10)
                        lresxy_err = '{0:{1}}'.format(dresxy_err, sfmt)
                    else:
                        lresxy_err = '{0:{1}}'.format((self.resxy_err/100.)/
                                                        np.log(10), sfmt)
                    
                    #--> Res_yx (TM)
                    if self.resyx_err == 'data':
                        dresyx_err = (sresyx_err[jj]/resyx)/np.log(10)
                        lresyx_err = '{0:{1}}'.format(dresyx_err, sfmt)
                    else:
                        lresyx_err = '{0:{1}}'.format((self.resyx_err/100.)/
                                                        np.log(10), sfmt)
                    
                    #phase_xy(TE)
                    if self.phasexy_err == 'data':
                        dphasexy_err = '{0:{1}}'.format(sphasexy_err[jj], sfmt)
                    else:
                        dphasexy_err = '{0:{1}}'.format((self.phasexy_err/100.)*
                                                         57/2., sfmt)
                        
                    #phase_yx (TM)
                    if self.phaseyx_err == 'data':
                        dphaseyx_err = '{0:{1}}'.format(sphaseyx_err[jj], sfmt)
                    else:
                        dphaseyx_err = '{0:{1}}'.format((self.phaseyx_err/100.)*
                                                        57/2., sfmt)
                    
                    #calculate log10 of resistivity as prescribed by OCCAM
                    lresyx = '{0:{1}}'.format(np.log10(resyx), sfmt)
                    lresxy = '{0:{1}}'.format(np.log10(resxy), sfmt)
                    
                    
                    #if include the tipper
                    if self.tipper_err != None:
                        if tip[jj,0,0].real == 0.0 or tip[jj,0,1] == 0.0:
                            tipyn='n'
                        else:
                            #calculate the projection angle for real and 
                            #imaginary
                            tipphir = np.arctan(tip[jj,0,0].real/
                                                tip[jj,0,1].real)-\
                                                np.deg2rad(self.proj_angle)
                            tipphii = np.arctan(tip[jj,0,0].imag/
                                                tip[jj,0,1].imag)-\
                                                np.deg2rad(self.proj_angle)
                            
                            #project the tipper onto the profile line
                            projtipr = np.sqrt(tip[jj,0,0].real**2+
                                             tip[jj,0,1].real**2)*\
                                             np.cos(tipphir)
                            projtipi = np.sqrt(tip[jj,0,0].imag**2+
                                             tip[jj,0,1].imag**2)*\
                                             np.cos(tipphii)
                                             
                            projtipr = '{0:{1}}'.format(projtipr, sfmt)
                            projtipi = '{0:{1}}'.format(projtipi, sfmt)
                                      
                            #error of tipper is a decimal percentage
                            projtiperr = '{0:{1}}'.format(self.tipper_err/100.,
                                                          sfmt)
                            
                            tipyn = 'y'
                            
                        
                    #make a list of lines to write to the data file
                    if self.model_mode == 'both':
                        reslst.append(ss.join([str(kk+1), str(nn), '1',
                                               lresxy, lresxy_err, '\n']))
                        reslst.append(ss.join([str(kk+1), str(nn), '2',
                                               phasexy, dphasexy_err, '\n']))

                        reslst.append(ss.join([str(kk+1), str(nn), '5',
                                               lresyx, lresyx_err, '\n']))
                        reslst.append(ss.join([str(kk+1),str(nn),'6',
                                               phaseyx, dphaseyx_err, '\n']))
                        if self.tipper_err != None and tipyn == 'y':
                            reslst.append(ss.join([str(kk+1),str(nn),'3',
                                                   projtipr, projtiperr, 
                                                   '\n']))
                            reslst.append(ss.join([str(kk+1),str(nn),'4',
                                                   projtipi, projtiperr, 
                                                   '\n']))
                    elif self.model_mode == 'TM':
                        reslst.append(ss.join([str(kk+1),str(nn),'5',
                                               lresyx, lresyx_err, '\n']))
                        reslst.append(ss.join([str(kk+1),str(nn),'6',
                                               phaseyx, dphaseyx_err, '\n']))
                        if self.tipper_err != None and tipyn == 'y':
                            reslst.append(ss.join([str(kk+1),str(nn),'3',
                                                   projtipr, projtiperr, 
                                                   '\n']))
                            reslst.append(ss.join([str(kk+1),str(nn),'4',
                                                   projtipi, projtiperr, 
                                                   '\n']))
                    elif self.model_mode == 'TE':
                        reslst.append(ss.join([str(kk+1),str(nn),'1',
                                               lresxy, lresxy_err, '\n']))
                        reslst.append(ss.join([str(kk+1),str(nn),'2',
                                               phasexy, dphasexy_err, '\n']))
                        if self.tipper_err != None and tipyn == 'y':
                            reslst.append(ss.join([str(kk+1),str(nn),'3',
                                                   projtipr, projtiperr,
                                                   '\n']))
                            reslst.append(ss.join([str(kk+1),str(nn),'4',
                                                   projtipi, projtiperr, 
                                                   '\n']))
                    else:
                        raise NameError('model_mode {0} not defined'.format(
                                                              self.model_mode))
        
        #======================================================================
        #                             write dat file
        #======================================================================
        if self.save_path!=None:
            if os.path.basename(self.save_path).find('.')>0:
                self.data_fn = self.save_path
            else:
                if not os.path.exists(self.save_path):
                    os.mkdir(self.save_path)
                self.data_fn=os.path.join(self.save_path,'Data.dat')
        else:
            self.data_fn=os.path.join(edipath,'Data.dat')
            
        if self.title==None:
            self.title='Occam Inversion'
            
        datfid=open(self.data_fn,'w')
        datfid.write('FORMAT:{0}OCCAM2MTDATA_1.0\n'.format(' '*11))
        datfid.write('TITLE:{0}{1} proj_angle={2:.4g}\n'.format(' '*12, 
                     self.title,
                     self.proj_angle))
        
        #write station sites
        datfid.write('SITES:{0}{1}\n'.format(' '*12, nstat))
        for station in station_lstsort:
            datfid.write('{0}{1}\n'.format(ss, station))
        
        #write offsets
        datfid.write('OFFSETS (M):\n')
        for offset in offsetlst:
            datfid.write('{0}{1: .2f}\n'.format(ss, offset))
        
        #write frequencies
        datfid.write('FREQUENCIES:{0}{1}\n'.format(' '*8, len(klst)))
        for fkey in klst:
            datfid.write('{0}{1:.5f}\n'.format(ss, float(fkey)))
        
        #write data block
        datfid.write('DATA BLOCKS:{0}{1}\n'.format(' '*10, len(reslst)))
        datfid.write(ss.join(['SITE', 'FREQ', 'TYPE', 'DATUM', 'ERROR', '\n']))
        for ll, datline in enumerate(reslst):
            if datline.find('#IND') >= 0:
                print 'Found #IND on line ',ll
                ndline = datline.replace('#IND','00')
                print 'Replaced with 00'
                datfid.write(ndline)
            elif datline.lower().find('inf') >= 0:
                print 'Found #inf on line ',ll
                ndline = datline.replace('#inf','00')
                print 'Replaced with 00'
                datfid.write(ndline)
            else:
                datfid.write(datline)
        datfid.close()
        
        print 'Wrote Occam2D data file to: ',self.data_fn
        
    def read2DdataFile(self):
        """
            read2DdataFile will read in data from a 2D occam data file.  
            Only supports the first 6 data types of occam2D
        
        Arguments:
        ----------
        
            **OccamPointPicker.data_fn** : full path to data file
        
        Returns:
        --------
            **rp_lst** : list of dictionaries for each station 
                                         with keywords:
                
                *'station'* : string
                              station name
                
                *'offset'* : float
                            relative offset
                
                *'resxy'* : np.array(nf,4)
                          TE resistivity and error as row 0 and 1 ressectively
                
                *'resyx'* : np.array(fn,4)
                          TM resistivity and error as row 0 and 1 respectively
                
                *'phasexy'* : np.array(nf,4)
                            TE phase and error as row 0 and 1 respectively
                
                *'phaseyx'* : np.array(nf,4)
                            Tm phase and error as row 0 and 1 respectively
                
                *'realtip'* : np.array(nf,4)
                            Real Tipper and error as row 0 and 1 respectively
                
                *'imagtip'* : np.array(nf,4)
                            Imaginary Tipper and error as row 0 and 1 
                            respectively
                
                Note: that the resistivity will be in log10 space.  Also, there
                are 2 extra rows in the data arrays, this is to put the 
                response from the inversion. 
                
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/Data.dat"
            >>> ocd.read2DdataFile()
            
        """
        
        dfid = open(self.data_fn,'r')
        
        dlines = dfid.readlines()
        #get format of input data
        self.occamfmt = dlines[0].strip().split(':')[1].strip()
        
        #get title
        self.titlestr = dlines[1].strip().split(':')[1].strip()
    
        if self.titlestr.find('=') > 0:
            tstr=self.titlestr.split('=')
            self.proj_angle = float(tstr[1])
            self.title = tstr[0]
        else:
            self.title = self.titlestr
            self.proj_angle = 0
            print 'Need to figure out angle of profile line'
        #get number of sits
        nsites = int(dlines[2].strip().split(':')[1].strip())
        
        #get station names
        self.station_lst = [dlines[ii].strip() for ii in range(3,nsites+3)]
        
        #get offsets in meters
        offsets = [float(dlines[ii].strip()) 
                    for ii in range(4+nsites,4+2*nsites)]
        
        #get number of frequencies
        nfreq = int(dlines[4+2*nsites].strip().split(':')[1].strip())
    
        #get frequencies
        self.freq = np.array([float(dlines[ii].strip()) 
                                for ii in range(5+2*nsites,5+2*nsites+nfreq)])
        
        #get periods
        self.period = 1./self.freq

        #-----------get data-------------------
        #set zero array size the first row will be the data and second the error
        asize = (4, nfreq)
        #make a list of dictionaries for each station.
        self.rp_lst=[{'station':station,'offset':offsets[ii],
                      'resxy':np.zeros(asize),
                      'resyx':np.zeros(asize),
                      'phasexy':np.zeros(asize),
                      'phaseyx':np.zeros(asize),
                      'realtip':np.zeros(asize),
                      'imagtip':np.zeros(asize)}
                      for ii,station in enumerate(self.station_lst)]
        for line in dlines[7+2*nsites+nfreq:]:
            ls = line.split()
            #station index
            ss = int(float(ls[0]))-1
            #component key
            comp = str(int(float(ls[2])))
            #frequency index        
            ff = int(float(ls[1]))-1
            #put into array
            #input data
            self.rp_lst[ss][occamdict[comp]][0,ff] = float(ls[3]) 
            #error       
            self.rp_lst[ss][occamdict[comp]][1,ff] = float(ls[4])
            
    def rewrite2DdataFile(self, **kwargs):
        """
        rewrite2DDataFile will rewrite an existing data file so you can 
        redefine some of the parameters, such as rotation angle, or errors for 
        the different components or only invert for one mode or add one or add
        tipper or remove tipper.
        
        Arguments:
        ----------
            
            **data_fn** : string
                         full path to data file to rewrite
            
            **thetar** : float
                         rotation angle with positive clockwise
            
            **resxy_err** : float
                           error for TE mode resistivity (percent) or 'data' 
                           for data or 'prev' to take errors from data file.
            
            **resyx_err** : float
                           error for TM mode resistivity (percent) or 'data' 
                           for data or 'prev' to take errors from data file.
                        
            **phasexy_err** : float
                             error for TE mode phase (percent) or 'data' 
                             for data or 'prev' to take errors from data file.
                        
            **phaseyx_err** : float
                             error for TM mode phase (percent) or 'data' 
                             for data or 'prev' to take errors from data file.
                        
            **tipper_err** : float 
                            error for tipper (percent) input only if you want
                            to invert for the tipper or 'data' for data errors
                            or prev to take errors from data file.
                        
            **model_mode** : string can be:
                            * 'both' for both TE and TM
                            * 'TE' for TE
                            * 'TM' for TM
                     
            **new_freq_lst** : list or np.array
                       frequency list in Hz to rewrite, needs to be similar to
                       the datafile, cannot add frequencies
                    
            **remove_station** : list of stations to remove if desired
            
            **save_path** : string
                           full path to save the file to, if None then file 
                           is saved to os.path.dirname(data_fn,'DataRW.dat')
            
            **new_proj_angle** : angle in degree relative to north, positive
                                 clockwise of new angle to project stations 
                                 onto.
            
        Returns:
        --------
        
            **Occam2DData.ndata_fn** : string
                                      full path to new data file
        
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/Data.dat"
            >>>
            >>> #rotate by 10 degrees east of North, remove station MT03 and 
            >>> #increase the resistivity error to 30 percent and put into a new 
            >>> #folder using save_path
            >>> svpath = r"/home/Occam2D/Line1/Inv2"
            >>> edipath = r"/home/EDIfiles"
            >>> ocd.rewrite2DdataFile(edipath=edipath,thetar=10,resxy_err=30,
            >>>                       resyx_err=30,removstation='MT03',
            >>>                       save_path=svpath)
            >>> Rewrote the data file to: /home/Occam2D/Line1/Inv2/DataRW.dat
            
        """
        
        #--> get information from key word arguments        
        self.edipath = kwargs.pop('edipath', self.edipath)
        self.thetar = kwargs.pop('thetar', 0)
        self.resxy_err = kwargs.pop('resxy_err', 'prev')
        self.resyx_err = kwargs.pop('resyx_err', 'prev')
        self.phasexy_err = kwargs.pop('phasexy_err','prev')
        self.phaseyx_err = kwargs.pop('phaseyx_err','prev')
        self.tipper_err = kwargs.pop('tipper_err', None)
        self.model_mode = kwargs.pop('model_mode', 'both')
        new_freq_lst = kwargs.pop('new_freq_lst', None)
        remove_station = kwargs.pop('remove_station', None)
        self.save_path = kwargs.pop('save_path', None)
        new_proj_angle = kwargs.pop('new_proj_angle', 0)
        
        ss = self._ss
        sfmt = self._string_fmt
        
        #load the data for the data file    
        self.read2DdataFile()
        
        #copy the information into local lists in case you want to keep the 
        #original data
        rp_lst = list(self.rp_lst)
        station_lst = list(self.station_lst)
        
        #make a dictionary of rp_lst for easier extraction of data
        rpdict = dict([(station, rp_lst[ii]) for ii,station in 
                                                    enumerate(station_lst)])
    
        #remove stations from rp_lst and station_lst if desired
        if remove_station!=None:
            #if remove_station is not a list make it one
            if type(remove_station) is not list:
                remove_station=[remove_station]
            
            #remove station from station list           
            for rstation in remove_station:        
                try:
                    station_lst.remove(rstation)
                except ValueError:
                    print 'Did not find {0}'.format(rstation)
        self.station_lst = station_lst
        
        #if flst is not the same as freq make freq=flst
        if new_freq_lst != None:
            self.freq = new_freq_lst
        
        #if the rotation angle is not 0 than need to read the original data in
        if self.thetar != 0:
            if self.thetar < 2*np.pi:
                self.thetar = np.rad2deg(self.thetar)
                
            if self.edipath == None:
                raise IOError('Need to input the edipath to original edifiles'
                               ' to get rotations correct')
            
            #get list of edifiles already in data file
            edilst = [os.path.join(self.edipath,edi) for stat in self.station_lst 
                      for edi in os.listdir(self.edipath) if edi[0:len(stat)]==stat]
            reslst = []
            for kk, edifn in enumerate(edilst,1):
                z1 = mtedi.Edi()
                z1.readfile(edifn)
                z1.Z.rotate(self.thetar)
                
                res, phase, res_err, phase_err = z1.Z.res_phase()
                tip = z1.Tipper.tipper
                if tip is None:
                    tip = np.zeros((len(self.freq), 1, 2))
                else:
                    z1.Tipper.rotate(self.thetar)
                    tip = z1.Tipper.tipper

                station = self.station_lst[kk-1]
                freq_dict = dict([('{0:.6g}'.format(fr),ii) for ii,fr in 
                                                    enumerate(z1.freq)])
                                                    
                #loop over frequencies to pick out the ones desired
                for jj, ff in enumerate(self.freq, 1):
                    #jj is the index of edi file frequency list, this index 
                    #corresponds
                    #to the impedance tensor component index
                    #ff is the frequency from the edi file frequency list
                    try:
                        #nn is the frequency number out of extracted frequency
                        #list
                        nn = freq_dict['%.6g' % ff]
                        
                        #calculate resistivity
                        resxy = res[nn, 0, 1]
                        resyx = res[nn, 1, 0]
                
                        #calculate the phase putting the yx in the 1st quadrant
                        phasexy = '{0:{1}}'.format(phase[nn, 0, 1], sfmt)
                        phaseyx = '{0:{1}}'.format(phase[nn, 1, 0]+180, sfmt)
                        #put phases in correct quadrant if should be negative
                        if float(phaseyx) > 180:
                            phaseyx = '{0:{1}}'.format(float(phaseyx)-360, 
                                                        sfmt)
                            print ('Found Negative Phase for station'
                                   '{0} frequency {1}'.format(z1.station, ff))    
                        
                        #calculate errors
                        #res_xy (TE)
                        if self.resxy_err == 'data':
                            lresxy_err = '{0:{1}}'.format(
                                    (res_err[nn, 0, 1]/resxy)/np.log(10), sfmt)
                        #take errors from data file
                        elif self.resxy_err == 'prev':
                            lresxy_err = '{0:{1}}'.format(
                                    rpdict[station]['resxy'][1,jj-1], sfmt)
                        else:
                            lresxy_err = '{0:{1}}'.format(
                                    (self.resxy_err/100.)/np.log(10), sfmt)
                        
                        #Res_yx(TM)
                        if self.resyx_err == 'data':
                            lresxy_err = '{0:{1}}'.format(
                                    (res_err[nn, 1, 0]/resyx)/np.log(10), sfmt)
                        #take errors from data file
                        elif self.resyx_err == 'prev':
                            lresyx_err = '{0:{1}}'.format(
                                    rpdict[station]['resyx'][1,jj-1], sfmt)
                        else:
                            lresyx_err = '{0:{1}}'.format(
                                    (self.resyx_err/100.)/np.log(10), sfmt)
                        
                        #phase_xy(TE)
                        if self.phasexy_err == 'data':
                            dphasexy_err = '{0:{1}}'.format(
                                    phase_err[nn, 0, 1], sfmt)
                            #take errors from data file
                        elif self.phasexy_err == 'prev':
                            dphasexy_err = '{0:{1}}'.format(
                                    rpdict[station]['phasexy'][1,jj-1], sfmt)
                        else:
                            dphasexy_err = '{0:{1}}'.format(
                                    (self.phasexy_err/100.)*57/2., sfmt)
                            
                        #phase_yx (TM)
                        if self.phaseyx_err == 'data':
                            dphaseyx_err = '{0:{1}}'.format(
                                    phase_err[nn, 1, 0], sfmt)
                        elif self.phaseyx_err == 'prev':
                            dphaseyx_err ='{0:{1}}'.format(
                                    rpdict[station]['phaseyx'][1,jj-1], sfmt)
                        else:
                            dphaseyx_err = '{0:{1}}'.format(
                                    (self.phaseyx_err/100.)*57/2., sfmt)
                        
                        #calculate log10 of resistivity as prescribed by OCCAM
                        lresyx = '{0:{1}}'.format(np.log10(resyx), sfmt)
                        lresxy = '{0:{1}}'.format(np.log10(resxy), sfmt)
                        
                        #if include the tipper
                        if self.tipper_err != None:
                            if tip[jj,0,0].real == 0.0 or tip[jj,0,1] == 0.0:
                                tipyn='n'
                            else:
                                #calculate the projection angle for real and 
                                #imaginary
                                tipphir = np.arctan(tip[jj,0,0].real/
                                                    tip[jj,0,1].real)-\
                                                    np.deg2rad(self.proj_angle)
                                tipphii = np.arctan(tip[jj,0,0].imag/
                                                    tip[jj,0,1].imag)-\
                                                    np.deg2rad(self.proj_angle)
                                
                                #project the tipper onto the profile line
                                projtipr = np.sqrt(tip[jj,0,0].real**2+
                                                 tip[jj,0,1].real**2)*\
                                                 np.cos(tipphir)
                                projtipi = np.sqrt(tip[jj,0,0].imag**2+
                                                 tip[jj,0,1].imag**2)*\
                                                 np.cos(tipphii)
                                                 
                                projtipr = '{0:{1}}'.format(projtipr, sfmt)
                                projtipi = '{0:{1}}'.format(projtipi, sfmt)
                                          
                                #error of tipper is a decimal percentage
                                projtiperr = '{0:{1}}'.format(self.tipper_err/
                                                              100., sfmt)
                                
                                tipyn = 'y'
                            
                        
                        #make a list of lines to write to the data file
                        if self.model_mode == 'both':
                            reslst.append(ss.join([str(kk+1), str(nn), '1',
                                                   lresxy, lresxy_err, '\n']))
                            reslst.append(ss.join([str(kk+1), str(nn), '2',
                                                   phasexy, dphasexy_err, 
                                                   '\n']))
    
                            reslst.append(ss.join([str(kk+1), str(nn), '5',
                                                   lresyx, lresyx_err, '\n']))
                            reslst.append(ss.join([str(kk+1),str(nn),'6',
                                                   phaseyx, dphaseyx_err, 
                                                   '\n']))
                            if self.tipper_err != None and tipyn == 'y':
                                reslst.append(ss.join([str(kk+1),str(nn),'3',
                                                       projtipr, projtiperr, 
                                                       '\n']))
                                reslst.append(ss.join([str(kk+1),str(nn),'4',
                                                       projtipi, projtiperr, 
                                                       '\n']))
                        elif self.model_mode == 'TM':
                            reslst.append(ss.join([str(kk+1),str(nn),'5',
                                                   lresyx, lresyx_err, '\n']))
                            reslst.append(ss.join([str(kk+1),str(nn),'6',
                                                   phaseyx, dphaseyx_err,
                                                   '\n']))
                            if self.tipper_err != None and tipyn == 'y':
                                reslst.append(ss.join([str(kk+1),str(nn),'3',
                                                       projtipr, projtiperr, 
                                                       '\n']))
                                reslst.append(ss.join([str(kk+1),str(nn),'4',
                                                       projtipi, projtiperr, 
                                                       '\n']))
                        elif self.model_mode == 'TE':
                            reslst.append(ss.join([str(kk+1),str(nn),'1',
                                                   lresxy, lresxy_err, '\n']))
                            reslst.append(ss.join([str(kk+1),str(nn),'2',
                                                   phasexy, dphasexy_err, 
                                                   '\n']))
                            if self.tipper_err != None and tipyn == 'y':
                                reslst.append(ss.join([str(kk+1),str(nn),'3',
                                                       projtipr, projtiperr,
                                                       '\n']))
                                reslst.append(ss.join([str(kk+1),str(nn),'4',
                                                       projtipi, projtiperr, 
                                                       '\n']))
                        else:
                            raise NameError('model_mode {0} not defined'.format(
                                            self.model_mode))
                    except KeyError:
                        pass
        
        #If no rotation is desired but error bars are than...
        else:
            reslst = []
            for kk, station in enumerate(self.station_lst,1):
                srp = rpdict[station]
                nr = srp['resxy'].shape[1]
                #calculate errors and rewrite
                #res_xy (TE)
                if self.resxy_err != None:
                    if self.resxy_err == 'prev':
                        lresxy_err = rpdict[station]['resxy'][1,:]
                    else:
                        lresxy_err = np.repeat((self.resxy_err/100.)/
                                                np.log(10), nr)
                    srp['resxy'][1,:] = lresxy_err
                
                #Res_yx(TM)
                if self.resyx_err != None:
                    if self.resyx_err == 'prev':
                        lresyx_err = rpdict[station]['resyx'][1,:]
                    else:
                        lresyx_err = np.repeat((self.resyx_err/100.)/
                                                np.log(10), nr)
                    srp['resyx'][1,:] = lresyx_err
                
                #phase_xy(TE)
                if self.phasexy_err != None:
                    if self.phasexy_err == 'prev':
                        dphasexy_err = rpdict[station]['phasexy'][1,:]
                    else:
                        dphasexy_err = np.repeat((self.phasexy_err/100.)*57/2.,
                                                 nr)
                    srp['phasexy'][1,:] = dphasexy_err
                    
                #phase_yx (TM)
                if self.phaseyx_err != None:
                    if self.phaseyx_err == 'prev':
                        dphaseyx_err = rpdict[station]['phaseyx'][1,:]
                    else:
                        dphaseyx_err = np.repeat((self.phaseyx_err/100.)*57/2.,
                                                 nr)
                    srp['phaseyx'][1,:] = dphaseyx_err
                
                if self.tipper_err != None:
                    #error of tipper is a decimal percentage
                    projtiperr = self.tipper_err/100.
                    srp['realtip'][1,:] = np.repeat(projtiperr,nr)
                    srp['imagtip'][1,:] = np.repeat(projtiperr,nr)
                
                for jj, ff in enumerate(self.freq, 1):
                    #make a list of lines to write to the data file
                    if self.model_mode == 'both':
                        if srp['resxy'][0,jj-1] != 0.0:
                            reslst.append(ss.join([str(kk), str(jj), '1',
                                  '{0:{1}}'.format(srp['resxy'][0,jj-1], sfmt),
                                  '{0:{1}}'.format(srp['resxy'][1,jj-1], sfmt),
                                  '\n']))
                        if srp['phasexy'][0,jj-1] != 0.0:
                            reslst.append(ss.join([str(kk), str(jj), '2',
                                 '{0:{1}}'.format(srp['phasexy'][0,jj-1],sfmt),
                                 '{0:{1}}'.format(srp['phasexy'][1,jj-1],sfmt),
                                 '\n']))
                        if srp['resyx'][0,jj-1]!=0.0:
                            reslst.append(ss.join([str(kk), str(jj), '5',
                                  '{0:{1}}'.format(srp['resyx'][0,jj-1],sfmt),
                                  '{0:{1}}'.format(srp['resyx'][1,jj-1],sfmt),
                                  '\n']))
                        if srp['phaseyx'][0,jj-1] != 0.0:
                            reslst.append(ss.join([str(kk), str(jj), '6',
                                 '{0:{1}}'.format(srp['phaseyx'][0,jj-1],sfmt),
                                 '{0:{1}}'.format(srp['phaseyx'][1,jj-1],sfmt),
                                 '\n']))
                        if self.tipper_err != None:
                            if srp['realtip'][0,jj-1] != 0.0:
                                reslst.append(ss.join([str(kk), str(jj), '3',
                                '{0:{1}}'.format(srp['realtip'][0,jj-1], sfmt),
                                '{0:{1}}'.format(srp['realtip'][1,jj-1], sfmt),
                                '\n']))
                            if srp['imagtip'][0,jj-1] != 0.0:
                                reslst.append(ss.join([str(kk), str(jj),'4',
                                 '{0:{1}}'.format(srp['imagtip'][0,jj-1],sfmt),
                                 '{0:{1}}'.format(srp['imagtip'][1,jj-1],sfmt),
                                 '\n']))
                                 
                    elif self.model_mode=='TM':
                        if srp['resyx'][0,jj-1]!=0.0:
                            reslst.append(ss.join([str(kk), str(jj), '5',
                                  '{0:{1}}'.format(srp['resyx'][0,jj-1],sfmt),
                                  '{0:{1}}'.format(srp['resyx'][1,jj-1],sfmt),
                                  '\n']))
                        if srp['phaseyx'][0,jj-1] != 0.0:
                            reslst.append(ss.join([str(kk), str(jj), '6',
                                 '{0:{1}}'.format(srp['phaseyx'][0,jj-1],sfmt),
                                 '{0:{1}}'.format(srp['phaseyx'][1,jj-1],sfmt),
                                 '\n']))
                        if self.tipper_err != None:
                            if srp['realtip'][0,jj-1] != 0.0:
                                reslst.append(ss.join([str(kk), str(jj), '3',
                                '{0:{1}}'.format(srp['realtip'][0,jj-1], sfmt),
                                '{0:{1}}'.format(srp['realtip'][1,jj-1], sfmt),
                                '\n']))
                            if srp['imagtip'][0,jj-1] != 0.0:
                                reslst.append(ss.join([str(kk), str(jj),'4',
                                 '{0:{1}}'.format(srp['imagtip'][0,jj-1],sfmt),
                                 '{0:{1}}'.format(srp['imagtip'][1,jj-1],sfmt),
                                 '\n']))
                                 
                    elif self.model_mode=='TE':
                        if srp['resxy'][0,jj-1] != 0.0:
                            reslst.append(ss.join([str(kk), str(jj), '1',
                                  '{0:{1}}'.format(srp['resxy'][0,jj-1], sfmt),
                                  '{0:{1}}'.format(srp['resxy'][1,jj-1], sfmt),
                                  '\n']))
                        if srp['phasexy'][0,jj-1] != 0.0:
                            reslst.append(ss.join([str(kk), str(jj), '2',
                                 '{0:{1}}'.format(srp['phasexy'][0,jj-1],sfmt),
                                 '{0:{1}}'.format(srp['phasexy'][1,jj-1],sfmt),
                                 '\n']))
                        if self.tipper_err != None:
                            if srp['realtip'][0,jj-1] != 0.0:
                                reslst.append(ss.join([str(kk), str(jj), '3',
                                '{0:{1}}'.format(srp['realtip'][0,jj-1], sfmt),
                                '{0:{1}}'.format(srp['realtip'][1,jj-1], sfmt),
                                '\n']))
                            if srp['imagtip'][0,jj-1] != 0.0:
                                reslst.append(ss.join([str(kk), str(jj),'4',
                                 '{0:{1}}'.format(srp['imagtip'][0,jj-1],sfmt),
                                 '{0:{1}}'.format(srp['imagtip'][1,jj-1],sfmt),
                                 '\n']))
    
        #===========================================================================
        #                             write dat file
        #===========================================================================
        
        #make the file name of the data file
        if self.data_fn.find('RW') > 0:
            if self.save_path == None:
                self.ndata_fn = self.data_fn
            elif os.path.isdir(self.save_path) == True:
                self.ndata_fn = os.path.join(self.save_path, 'DataRW.dat')
            elif os.path.isfile(self.save_path) == True:
                self.ndata_fn = self.save_path
            elif self.save_path.find('.dat') > 0:
                self.ndata_fn = self.save_path
        else:
            if self.save_path == None:
                self.ndata_fn = self.data_fn[:-4]+'RW.dat'
            elif os.path.isdir(self.save_path) == True:
                self.ndata_fn = os.path.join(self.save_path,'DataRW.dat')
            elif os.path.isfile(self.save_path) == True:
                self.ndata_fn = self.save_path
            elif self.save_path.find('.dat') > 0:
                self.ndata_fn = self.save_path
            
        nstat = len(self.station_lst)
            
        if self.titlestr == None:
            self.titlestr = 'Occam Inversion'
            
        datfid = open(self.ndata_fn,'w')
        datfid.write('FORMAT:{0}OCCAM2MTDATA_1.0\n'.format(' '*11))
        if new_proj_angle != 0:
            self.titlestr = 'Occam Inversion proj_angle={0:.1f}'.format(
                                                               new_proj_angle)
                                                               
        datfid.write('TITLE:{0}{1}\n'.format(' '*12,self.titlestr))
        
        #write station sites
        datfid.write('SITES:{0}{1}\n'.format(' '*12, nstat))
        for station in self.station_lst:
            datfid.write('{0}{1}\n'.format(ss, station))
        
        #write offsets
        datfid.write('OFFSETS (M): \n')
        for station in self.station_lst:
            if new_proj_angle != 0:
                new_offset = rpdict[station]['offset']/\
                             np.cos(np.deg2rad(self.proj_angle-new_proj_angle))
                #need to project the stations on to the strike direction
                datfid.write('{0}{1: .2f}\n'.format(ss,new_offset))
            else:
                datfid.write('{0}{1: .2f}\n'.format(ss,
                                                    rpdict[station]['offset']))
        
        #write frequencies
        datfid.write('FREQUENCIES:{0}{1}\n'.format(' '*8, len(self.freq)))
        for ff in self.freq:
            datfid.write('{0}{1:.5f}\n'.format(ss, ff))
        
        #write data block
        datfid.write('DATA BLOCKS:{0}{1}\n'.format(' '*10, len(reslst)))
        datfid.write(ss.join(['SITE', 'FREQ', 'TYPE', 'DATUM', 'ERROR', '\n']))
        for ll, datline in enumerate(reslst):
            if datline.find('#IND') >= 0:
                print 'Found #IND on line {0}'.format(ll)
                ndline = datline.replace('#IND','00')
                print 'Replaced with 00'
                datfid.write(ndline)
            elif datline.lower().find('inf') >= 0:
                print 'Found #inf on line {0}'.format(ll)
                ndline = datline.replace('#inf','00')
                print 'Replaced with 00'
                datfid.write(ndline)
            else:
                datfid.write(datline)
        datfid.close()
        
        #need to reset the projection angle
        if new_proj_angle != 0:
            self.proj_angle = new_proj_angle
        
        print 'Rewrote the data file to: ',self.ndata_fn
        
    def makeModelFiles(self,niter=20,targetrms=1.0,nlayers=100,nlperdec=30,
              z1layer=50,bwidth=200,trigger=.75,save_path=None,rhostart=100,
              occampath=r"c:\Peacock\PHD\OCCAM\MakeFiles"):
        """
        makeModel will make an the input files for occam using Steve Constable's
        MakeModel2DMT.f code.
        
        Inputs:
            data_fn = full path to data file
            niter = maximum number of iterations
            targetrms = target root mean square error
            nlayers = total number of layers in mesh
            nlperdec = number of layers per decade
            z1layer = thickness of the first layer in meters
            bwidth = maximum block width for regularization grid in meters
            trigger = triger point to amalgamate blocks
            save_path = path to save files to
            rhostart = starting resistivity for homogeneous half space in ohm-m
            occampath = path to MakeModel2DMT.exe
            
        Outputs:
            meshfn = mesh file for finite element grid saved ats MESH
            inmodelfn = input model, starting model with rhostart as starting value
                        saved as INMODEL
            startupfn = start up filepath, saved as startup
        """
        #get the base name of data file    
        dfnb=os.path.basename(self.data_fn)
        
        #put data file into the same directory as MakeModel2DMT
        if os.path.dirname(self.data_fn)!=occampath:
            shutil.copy(self.data_fn,os.path.join(occampath,dfnb))
        
        #write input file for MakeModel2DMT
        mmfid=open(os.path.join(occampath,'inputMakeModel.txt'),'w')
        mmfid.write(dfnb, '\n')
        mmfid.write(str(niter), '\n')    
        mmfid.write(str(targetrms), '\n')    
        mmfid.write(str(nlayers), '\n')
        mmfid.write(str(nlperdec), '\n')
        mmfid.write(str(z1layer), '\n')
        mmfid.write(str(bwidth), '\n')
        mmfid.write(str(trigger), '\n')
        mmfid.write('\n')
        mmfid.close()
        
        #get current working directory
        cdir=os.getcwd() 
        
        #change directory path to occam path
        os.chdir(occampath) 
        
        #---call MakeModel2DMT---
        subprocess.os.system("MakeModel2DMT < inputMakeModel.txt")
        
        #change back to original working directory    
        os.chdir(cdir)
        
        if save_path==None:
            save_path=os.path.dirname(self.data_fn)
        
        if not os.path.exists(save_path):
            os.mkdir(save_path)
        
        meshfn=os.path.join(save_path,'MESH')    
        inmodelfn=os.path.join(save_path,'INMODEL')    
        startupfn=os.path.join(save_path,'startup')    
        
        #copy ouput files to save_path
        shutil.copy(os.path.join(occampath,'MESH'),meshfn)
        shutil.copy(os.path.join(occampath,'INMODEL'),inmodelfn)
        shutil.copy(os.path.join(occampath,'startup'),startupfn)
        shutil.copy(os.path.join(occampath,'inputMakeModel.txt'),
                    os.path.join(save_path,'inputMakeModel.txt'))
        if not os.path.exists(os.path.join(save_path,dfnb)):
            shutil.copy(self.data_fn,os.path.join(save_path,dfnb))
        if os.path.getctime(os.path.join(save_path,dfnb))<\
            os.path.getctime(self.data_fn):
            shutil.copy(self.data_fn,os.path.join(save_path,dfnb))
        
        
#        #rewrite mesh so it contains the right number of columns and rows
#        rewriteMesh(meshfn)
        
        #write startup file to have the starting desired starting rho value
        ifid=open(startupfn,'r')
        ilines=ifid.readlines()
        ifid.close()
        
        if rhostart!=100:
            #make startup model a homogeneous half space of rhostart
            rhostart=np.log10(rhostart)
            ifid=open(startupfn,'w')
            for line in ilines:
                if line.find('2.000000')>=0:
                    line=line.replace('2.000000','%.6f' % rhostart)
                ifid.write(line)
        ifid.close()
        
        print 'Be sure to check the INMODEL file for clumped numbers near the bottom.'
        print 'Also, check the MESH and startup files to make sure they are correct.'
        
        self.meshfn=meshfn
        self.inmodelfn=inmodelfn
        self.startupfn=startupfn
    
    def plotMaskPoints(self,plottype=None,reserrinc=.20,phaseerrinc=.05,
                       marker='h',colormode='color',dpi=300,ms=2,
                       res_limits=None,phaselimits=(-5,95)):
        """
        An interactive plotting tool to mask points an add errorbars
        
        Arguments:
        ----------
            **plottype** : string
                           describes the way the responses are plotted. Can be:
                               *list of stations to plot ['mt01','mt02']
                               *one station 'mt01'
                               *None plots all stations *Default*
            
            **reserrinc** : float
                            amount to increase the error bars. Input as a 
                            decimal percentage.  0.3 for 30 percent
                            *Default* is 0.2 (20 percent)
                            
            **phaseerrinc** : float
                              amount to increase the error bars. Input as a 
                              decimal percentage.  0.3 for 30 percent
                              *Default* is 0.05 (5 percent)
                            
            **marker** : string
                         marker that the masked points will be
                         *Default* is 'h' for hexagon
                        
            **colormode** : string
                            defines the color mode to plot the responses in
                            *'color'* for color plots
                            *'bw'* for black and white plots
                            *Default* is 'color'
            
            **dpi** : int
                      dot-per-inch resolution of the plots.
                      *Default* is 300
            
            **ms** : float
                     size of the marker in the response plots
                     *Default* is 2
                     
            **self.res_limits**: tuple (min,max)
                           min and max limits of the resistivity values in 
                           linear scale
                           *Default* is None
                           
            **phaselimits**: tuple (min,max)
                            min and max phase limits 
                            *Default* is (-5,90)
                            
        Returns:
        ---------
            data type **OccamPointPicker**  
                           
        
        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/Data.dat"
            >>> ocd.plotMaskPoints()                   
            
        """
        
        if colormode=='color':
            #color for data
            cted=(0,0,1)
            ctmd=(1,0,0)
            mted='s'
            mtmd='o'
            
        elif colormode=='bw':
            #color for data
            cted=(0,0,0)
            ctmd=(0,0,0)
            mted='s'
            mtmd='o'
            
        #read in data file    
        self.read2DdataFile()
        rp_lst=list(self.rp_lst)
        
        #get periods
        period=self.period
     
        #define some empty lists to put things into
        pstation_lst=[]
        axlst=[]
        linelst=[]
        errlst=[]
        
        #get the stations to plot
        #if none plot all of them
        if plottype==None:
            pstation_lst=range(len(self.station_lst))
            
        #otherwise pick out the stations to plot along with their index number
        elif type(plottype) is not list:
            plottype=[plottype]
            for ii,station in enumerate(self.station_lst):
                for pstation in plottype:
                    if station.find(pstation)>=0:
                        pstation_lst.append(ii) 
        
        #set the subplot grid
        gs=gridspec.GridSpec(6,2,wspace=.1,left=.1,top=.93,bottom=.07)
        for jj,ii in enumerate(pstation_lst):
            fig=plt.figure(ii+1,dpi=dpi)
            plt.clf()
            
            #make subplots
            axrte=fig.add_subplot(gs[:4,0])
            axrtm=fig.add_subplot(gs[:4,1])
            axpte=fig.add_subplot(gs[-2:,0],sharex=axrte)    
            axptm=fig.add_subplot(gs[-2:,1],sharex=axrtm)    
            
            
            #plot resistivity TE Mode
            #cut out missing data points first
            rxy=np.where(rp_lst[ii]['resxy'][0]!=0)[0]
            rte=axrte.errorbar(period[rxy],10**rp_lst[ii]['resxy'][0][rxy],
                            ls=':',marker=mted,ms=ms,mfc=cted,mec=cted,
                            color=cted,
                            yerr=np.log(10)*rp_lst[ii]['resxy'][1][rxy]*\
                                10**rp_lst[ii]['resxy'][0][rxy],
                            ecolor=cted,picker=2)
                            
            #plot Phase TE Mode
            #cut out missing data points first
            pxy=[np.where(rp_lst[ii]['phasexy'][0]!=0)[0]]
            pte=axpte.errorbar(period[pxy],rp_lst[ii]['phasexy'][0][pxy],
                               ls=':',marker=mted,ms=ms,mfc=cted,mec=cted,
                               color=cted,yerr=rp_lst[ii]['phasexy'][1][pxy],
                               ecolor=cted,picker=1) 
            
                           
            #plot resistivity TM Mode
            #cut out missing data points first                
            ryx=np.where(rp_lst[ii]['resyx'][0]!=0)[0]
            rtm=axrtm.errorbar(period[ryx],10**rp_lst[ii]['resyx'][0][ryx],
                            ls=':',marker=mtmd,ms=ms,mfc=ctmd,mec=ctmd,
                            color=ctmd,
                            yerr=np.log(10)*rp_lst[ii]['resyx'][1][ryx]*\
                                10**rp_lst[ii]['resyx'][0][ryx],
                            ecolor=ctmd,picker=2)
            #plot Phase TM Mode
            #cut out missing data points first
            pyx=[np.where(rp_lst[ii]['phaseyx'][0]!=0)[0]]
            ptm=axptm.errorbar(period[pyx],rp_lst[ii]['phaseyx'][0][pyx],
                            ls=':',marker=mtmd,ms=ms,mfc=ctmd,mec=ctmd,
                            color=ctmd,yerr=rp_lst[ii]['phaseyx'][1][pyx],
                            ecolor=ctmd,picker=1)
        
        
            #make the axis presentable
            #set the apparent resistivity scales to log and x-axis to log
            axplst=[axrte,axrtm,axpte,axptm]
            llst=[rte[0],rtm[0],pte[0],ptm[0]]
            elst=[[rte[1][0],rte[1][1],rte[2][0]],
                  [rtm[1][0],rtm[1][1],rtm[2][0]],
                  [pte[1][0],pte[1][1],pte[2][0]],
                  [ptm[1][0],ptm[1][1],ptm[2][0]]]
                
            axlst.append(axplst)
            linelst.append(llst)
            errlst.append(elst)
            
            #set the axes properties for each subplot
            for nn,xx in enumerate(axplst):
                #set xscale to logarithmic in period
                xx.set_xscale('log')
                
                #if apparent resistivity 
                if nn==0 or nn==1:
                    #set x-ticklabels to invisible
                    plt.setp(xx.xaxis.get_ticklabels(),visible=False)
                    
                    #set apparent resistivity scale to logarithmic
                    try:
                        xx.set_yscale('log')
                    except ValueError:
                        pass
                    
                    #if there are resistivity limits set those
                    if self.res_limits!=None:
                        xx.set_ylim(self.res_limits)
                    
                #Set the title of the TE plot                 
                if nn==0:
                    xx.set_title(self.station_lst[ii]+' Obs$_{xy}$ (TE-Mode)',
                                 fontdict={'size':9,'weight':'bold'})
                    xx.yaxis.set_label_coords(-.075,.5)
                    xx.set_ylabel('App. Res. ($\Omega \cdot m$)',
                                  fontdict={'size':9,'weight':'bold'})
                #set the title of the TM plot
                if nn==1:
                    xx.set_title(self.station_lst[ii]+' Obs$_{yx}$ (TM-Mode)',
                                 fontdict={'size':9,'weight':'bold'})
                
                #set the phase axes properties
                if nn==2 or nn==3:
                    #set the phase limits
                    xx.set_ylim(phaselimits)
                    
                    #set label coordinates
                    xx.yaxis.set_label_coords(-.075,.5)
                    
                    #give the y-axis label to the bottom left plot
                    if nn==2:
                        xx.set_ylabel('Phase (deg)',
                                       fontdict={'size':9,'weight':'bold'})
                    #set the x-axis label
                    xx.set_xlabel('Period (s)',
                                  fontdict={'size':9,'weight':'bold'})
                    
                    #set tick marks of the y-axis
                    xx.yaxis.set_major_locator(MultipleLocator(10))
                    xx.yaxis.set_minor_locator(MultipleLocator(2))
                    
                xx.grid(True,alpha=.4,which='both') 
        
        #make points an attribute of self which is a data type OccamPointPicker       
        self.points=OccamPointPicker(axlst,linelst,errlst,reserrinc=reserrinc,
                                     phaseerrinc=phaseerrinc,marker=marker)
        
        #be sure to show the plot
        plt.show()
   
        
    def maskPoints(self):
        """
        maskPoints will take in points found from plotMaskPoints and rewrite 
        the data file to nameRW.dat.  **Be sure to run plotMaskPoints first**
        
        Arguments:
        ---------
            None
            
        Returns:
        ---------
        
            **OccamPointPicker.ndata_fn** : full path to rewritten data file
            
        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/Data.dat"
            >>> ocd.plotMaskPoints()
            >>>
            >>> #after all the points are masked rewrite the data file
            >>> ocd.maskPoints()
            >>>
            >>> Rewrote Occam2D data file to: /home/Occam2D/Line1/Inv1/DataRW.dat 
            
        """
        
        self.read2DdataFile()
        
        rp_lst = list(self.rp_lst)
        #rewrite the data file
        #make a reverse dictionary for locating the masked points in the data 
        #file
        rploc = dict([('{0}'.format(self.points.fndict[key]),int(key)-1) 
                    for key in self.points.fndict.keys()])
                    
        #make a period dictionary to locate points changed
        frpdict = dict([('{0:.5g}'.format(fr),ff) 
                            for ff,fr in enumerate(1./self.freq)])
        
        #loop over the data list
        for dd, dat in enumerate(self.points.data):
            derror = self.points.error[dd]
            #loop over the 4 main entrie
            for ss, skey in enumerate(['resxy', 'resyx', 'phasexy','phaseyx']):
                #rewrite any coinciding points
                for frpkey in frpdict.keys():
                    try:
                        ff = frpdict[frpkey]
                        floc = self.points.fdict[dd][ss][frpkey]
                        
                        #CHANGE APPARENT RESISTIVITY
                        if ss == 0 or ss == 1:
                            #change the apparent resistivity value
                            if rp_lst[rploc[str(dd)]][skey][0][ff] != \
                                                      np.log10(dat[ss][floc]):
                                if dat[ss][floc] == 0:
                                    rp_lst[rploc[str(dd)]][skey][0][ff] = 0.0
                                else:
                                    rp_lst[rploc[str(dd)]][skey][0][ff] = \
                                            np.log10(dat[ss][floc])
                                
                            #change the apparent resistivity error value
                            if dat[ss][floc] == 0.0:
                                rerr = 0.0
                            else:
                                rerr = derror[ss][floc]/dat[ss][floc]/np.log(10)
                            if rp_lst[rploc[str(dd)]][skey][1][ff] != rerr:
                                rp_lst[rploc[str(dd)]][skey][1][ff] = rerr
                        
                        #DHANGE PHASE
                        elif ss == 2 or ss == 3:
                            #change the phase value
                            if rp_lst[rploc[str(dd)]][skey][0][ff] != \
                                                                 dat[ss][floc]:
                                if dat[ss][floc] == 0:
                                    rp_lst[rploc[str(dd)]][skey][0][ff] = 0.0
                                else:
                                    rp_lst[rploc[str(dd)]][skey][0][ff] = \
                                                                  dat[ss][floc]
                                
                            #change the apparent resistivity error value
                            if dat[ss][floc] == 0.0:
                                rerr = 0.0
                            else:
                                rerr = derror[ss][floc]
                            if rp_lst[rploc[str(dd)]][skey][1][ff] != rerr:
                                rp_lst[rploc[str(dd)]][skey][1][ff] = rerr
                    except KeyError:
                        pass
            
            
        #rewrite the data file 
        ss = self._ss
        sfmt = self._string_fmt
        reslst = []
        
        #make a dictionary of rp_lst for easier extraction of data
        rpdict = dict([(station, rp_lst[ii]) 
                    for ii, station in enumerate(self.station_lst)])
        
        #loop over stations in the data file
        for kk, station in enumerate(self.station_lst,1):
            srp = rpdict[station]
            
            #loop over frequencies
            for jj, ff in enumerate(self.freq, 1):
                #make a list of lines to write to the data file
                if srp['resxy'][0,jj-1] != 0.0:
                    reslst.append(ss.join([str(kk), str(jj), '1',
                                  '{0:{1}}'.format(srp['resxy'][0,jj-1], sfmt),
                                  '{0:{1}}'.format(srp['resxy'][1,jj-1], sfmt),
                                  '\n']))
                                  
                if srp['phasexy'][0,jj-1] != 0.0:
                    reslst.append(ss.join([str(kk), str(jj), '2',
                                 '{0:{1}}'.format(srp['phasexy'][0,jj-1],sfmt),
                                 '{0:{1}}'.format(srp['phasexy'][1,jj-1],sfmt),
                                 '\n']))
                                 
                if srp['resyx'][0,jj-1] != 0.0:
                    reslst.append(ss.join([str(kk), str(jj), '5',
                                  '{0:{1}}'.format(srp['resyx'][0,jj-1],sfmt),
                                  '{0:{1}}'.format(srp['resyx'][1,jj-1],sfmt),
                                  '\n']))
                                  
                if srp['phaseyx'][0,jj-1] != 0.0:
                    reslst.append(ss.join([str(kk), str(jj), '6',
                                 '{0:{1}}'.format(srp['phaseyx'][0,jj-1],sfmt),
                                 '{0:{1}}'.format(srp['phaseyx'][1,jj-1],sfmt),
                                 '\n']))
                                 
                if srp['realtip'][0,jj-1] != 0.0:
                    reslst.append(ss.join([str(kk), str(jj), '3',
                                '{0:{1}}'.format(srp['realtip'][0,jj-1], sfmt),
                                '{0:{1}}'.format(srp['realtip'][1,jj-1], sfmt),
                                '\n']))
                                
                if srp['imagtip'][0,jj-1] != 0.0:
                    reslst.append(ss.join([str(kk), str(jj),'4',
                                 '{0:{1}}'.format(srp['imagtip'][0,jj-1],sfmt),
                                 '{0:{1}}'.format(srp['imagtip'][1,jj-1],sfmt),
                                 '\n']))
        
        #======================================================================
        #                             write dat file
        #======================================================================
        #make the file name of the data file
        if self.data_fn.find('RW') > 0:
            self.ndata_fn = self.data_fn
        else:
            self.ndata_fn = self.data_fn[:-4]+'RW.dat'
        
        #get number of stations
        nstat = len(self.station_lst)
        
        #set title string
        if self.titlestr == None:
            self.titlestr = 'Occam Inversion'
            
        datfid = open(self.ndata_fn,'w')
        datfid.write('FORMAT:{0}OCCAM2MTDATA_1.0\n'.format(' '*11))
        datfid.write('TITLE:{0}{1}\n'.format(' '*12,self.titlestr))
        
        #write station sites
        datfid.write('SITES:{0}{1}\n'.format(' '*12, nstat))
        for station in self.station_lst:
            datfid.write('{0}{1}\n'.format(ss, station))
        
        #write offsets
        datfid.write('OFFSETS (M): \n')
        for station in self.station_lst:
            datfid.write('{0}{1: .2f}\n'.format(ss,
                                             rpdict[station]['offset']))
        
        #write frequencies
        datfid.write('FREQUENCIES:{0}{1}\n'.format(' '*8, len(self.freq)))
        for ff in self.freq:
            datfid.write('{0}{1:.5f}\n'.format(ss, ff))
        
        #write data block
        datfid.write('DATA BLOCKS:{0}{1}\n'.format(' '*10, len(reslst)))
        datfid.write(ss.join(['SITE', 'FREQ', 'TYPE', 'DATUM', 'ERROR', '\n']))
        for ll, datline in enumerate(reslst):
            if datline.find('#IND') >= 0:
                print 'Found #IND on line {0}'.format(ll)
                ndline = datline.replace('#IND','00')
                print 'Replaced with 00'
                datfid.write(ndline)
            elif datline.lower().find('inf') >= 0:
                print 'Found #inf on line {0}'.format(ll)
                ndline = datline.replace('#inf','00')
                print 'Replaced with 00'
                datfid.write(ndline)
            else:
                datfid.write(datline)
        datfid.close()
        
        print 'Wrote Occam2D data file to: ',self.ndata_fn
    
    def read2DRespFile(self, resp_fn):
        """
        read2DRespFile will read in a response file and combine the data with info 
        from the data file.
    
        Arguments:
        ----------
            **resp_fn** : full path to the response file
            
            **data_fn** : full path to data file
    
        Returns:
        --------
            for each data array, the rows are ordered as:
                - 0 -> input data
                - 1 -> input error
                - 2 -> model output
                - 3 -> relative error (data-model)/(input error)
                
            **rp_lst** : list of dictionaries for each station with keywords:
                
                *station* : string
                            station name
                
                *offset* : float
                            relative offset
                
                *resxy* : np.array(nf,4)
                          TE resistivity and error as row 0 and 1 ressectively
                
                *resyx* : np.array(fn,4)
                          TM resistivity and error as row 0 and 1 respectively
                
                *phasexy* : np.array(nf,4)
                            TE phase and error as row 0 and 1 respectively
                
                *phaseyx* : np.array(nf,4)
                            Tm phase and error as row 0 and 1 respectively
                
                *realtip* : np.array(nf,4)
                            Real Tipper and error as row 0 and 1 respectively
                
                *imagtip* : np.array(nf,4)
                            Imaginary Tipper and error as row 0 and 1 
                            respectively
                
                Note: that the resistivity will be in log10 space.  Also, there
                are 2 extra rows in the data arrays, this is to put the 
                response from the inversion.  
            
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.read2DRespFile(r"/home/Occam2D/Line1/Inv1/Test_15.resp")
        """
        #make the response file an attribute        
        self.resp_fn = resp_fn
        
        #read in the current data file
        self.read2DdataFile()
        
        if self.resp_fn is None:
            return
        
        rfid = open(self.resp_fn, 'r')
        
        rlines = rfid.readlines()
        for line in rlines:
            ls = line.split()
            #station index
            ss = int(float(ls[0]))-1
            #component key
            comp = str(int(float(ls[2])))
            #frequency index        
            ff = int(float(ls[1]))-1
            #put into array
            #model response
            self.rp_lst[ss][occamdict[comp]][2,ff] = float(ls[5]) 
            #relative error        
            self.rp_lst[ss][occamdict[comp]][3,ff] = float(ls[6])
            
    def plot2DResponses(self, resp_fn=None, **kwargs):
        """
        plotResponse will plot the responses modeled from winglink against the 
        observed data.
        
        Arguments:
        ----------
        
            **resp_fn** : string
                         full path to response file
             
        ==================== ==================================================
        key words            description
        ==================== ==================================================
        color_mode           [ 'color' | 'bw' ] plot figures in color or 
                             black and white ('bw')
        cted                 color of Data TE marker and line
        ctem                 color of Model TE marker and line
        ctewl                color of Winglink Model TE marker and line
        ctmd                 color of Data TM marker and line
        ctmm                 color of Model TM marker and line
        ctmwl                color of Winglink Model TM marker and line
        e_capsize            size of error bar caps in points
        e_capthick           line thickness of error bar caps in points 
        fig_dpi              figure resolution in dots-per-inch 
        fig_lst              list of dictionaries with key words
                             station --> station name
                             fig --> matplotlib.figure instance
                             axrte --> matplotlib.axes instance for TE app.res
                             axrtm --> matplotlib.axes instance for TM app.res
                             axpte --> matplotlib.axes instance for TE phase
                             axptm --> matplotlib.axes instance for TM phase
                 
        fig_num              starting number of figure
        fig_size             size of figure in inches (width, height)
        font_size            size of axes ticklabel font in points
        lw                   line width of lines in points
        ms                   marker size in points
        mted                 marker for Data TE mode
        mtem                 marker for Model TE mode
        mtewl                marker for Winglink Model TE
        mtmd                 marker for Data TM mode
        mtmm                 marker for Model TM mode
        mtmwl                marker for Winglink TM mode
        period               np.ndarray of periods to plot 
        phase_limits         limits on phase plots in degrees (min, max)
        plot_num             [ 1 | 2 ] 
                             1 to plot both modes in a single plot
                             2 to plot modes in separate plots (default)
        plot_type            [ '1' | station_lst]
                             '1' --> to plot all stations in different figures
                             station_lst --> to plot a few stations, give names
                             of stations ex. ['mt01', 'mt07']
        plot_yn              [ 'y' | 'n']
                             'y' --> to plot on instantiation
                             'n' --> to not plot on instantiation
        res_limits           limits on resistivity plot in log scale (min, max)
        rp_lst               list of dictionaries from read2Ddata
        station_lst          station_lst list of stations in rp_lst
        subplot_bottom       subplot spacing from bottom (relative coordinates) 
        subplot_hspace       vertical spacing between subplots
        subplot_left         subplot spacing from left  
        subplot_right        subplot spacing from right
        subplot_top          subplot spacing from top
        subplot_wspace       horizontal spacing between subplots
        wl_fn                Winglink file name (full path)
        ==================== ==================================================
                            
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> rfile = r"/home/Occam2D/Line1/Inv1/Test_15.resp"
            >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
            >>>
            >>> #plot all responses in their own figure and modes in separate 
            >>> #subplots
            >>> rp = ocd.plot2DResponses(resp_fn=rfile)
                      
        """

        self.read2DRespFile(resp_fn)
            
        return PlotOccam2DResponse(self.rp_lst, self.period, **kwargs)
    
    
    def plotPseudoSection(self, resp_fn=None, **kwargs):
        """
        plots a pseudo section of the data and response if input.
        
        Arguments:
        ----------
            **resp_fn** : string
                         full path to response file
        
        ==================== ==================================================
        key words            description
        ==================== ==================================================
        axmpte               matplotlib.axes instance for TE model phase
        axmptm               matplotlib.axes instance for TM model phase
        axmrte               matplotlib.axes instance for TE model app. res 
        axmrtm               matplotlib.axes instance for TM model app. res 
        axpte                matplotlib.axes instance for TE data phase 
        axptm                matplotlib.axes instance for TM data phase
        axrte                matplotlib.axes instance for TE data app. res.
        axrtm                matplotlib.axes instance for TM data app. res.
        cb_pad               padding between colorbar and axes
        cb_shrink            percentage to shrink the colorbar to
        fig                  matplotlib.figure instance
        fig_dpi              resolution of figure in dots per inch
        fig_num              number of figure instance
        fig_size             size of figure in inches (width, height)
        font_size            size of font in points
        label_lst            list to label plots
        ml                   factor to label stations if 2 every other station
                             is labeled on the x-axis
        period               np.array of periods to plot
        phase_cmap           color map name of phase
        phase_limits_te      limits for te phase in degrees (min, max)
        phase_limits_tm      limits for tm phase in degrees (min, max)            
        plot_resp            [ 'y' | 'n' ] to plot response
        plot_yn              [ 'y' | 'n' ] 'y' to plot on instantiation
        res_cmap             color map name for resistivity
        res_limits_te        limits for te resistivity in log scale (min, max)
        res_limits_tm        limits for tm resistivity in log scale (min, max)
        rp_lst               list of dictionaries as made from read2Dresp
        station_id           index to get station name (min, max)
        station_lst          station list got from rp_lst
        subplot_bottom       subplot spacing from bottom (relative coordinates) 
        subplot_hspace       vertical spacing between subplots
        subplot_left         subplot spacing from left  
        subplot_right        subplot spacing from right
        subplot_top          subplot spacing from top
        subplot_wspace       horizontal spacing between subplots
        ==================== ==================================================
                            
       :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> rfile = r"/home/Occam2D/Line1/Inv1/Test_15.resp"
            >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
            >>> ps1 = ocd.plot2PseudoSection(resp_fn=rfile) 
        
        """
        
        self.read2DRespFile(resp_fn)
        
        if resp_fn is not None:
            plot_resp = 'y'
        else:
            plot_resp = 'n'
            
        return PlotPseudoSection(self.rp_lst, self.period, 
                                 plot_resp=plot_resp, **kwargs)
    
    def plotAllResponses(self, station_lst=None, **kwargs):
        """
        Plot all the responses of occam inversion from data file.  This assumes
        the response curves are in the same folder as the datafile.
    
        Arguments:
        ----------
            **station** : string
                          station name to plot
            
            **fignum** : int
                         plot number to put figure into
                         
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
            >>> ocd.plotAllResponses('MT01')
        
        """    
        dir_path = os.path.dirname(self.data_fn)
        resp_fn_lst = [os.path.join(dir_path, rfn) 
                       for rfn in os.listdir(dir_path)
                       if rfn.find('.resp')>0]
        plot_rp_lst = []
        for rfn in resp_fn_lst:
            self.read2DRespFile(rfn)
            plot_rp_lst.append(self.rp_lst)
        
        return PlotAllResponses(plot_rp_lst, self.period, 
                                pstation_lst=station_lst,
                                **kwargs)
        
class Occam2DModel(Occam2DData):
    """
    This class deals with the model side of Occam inversions, including 
    plotting the model, the L-curve, depth profiles.  It will also be able to 
    build a mesh and regularization grid at some point.  
    
    It inherits Occam2DData and the data can be extracted from the method
    get2DData().  After this call you can use all the methods of Occam2DData,
    such as plotting the model responses and pseudo sections.
    
    
    """
    
    def __init__(self,iterfn,meshfn=None,inmodelfn=None):
        self.iterfn=iterfn
    
        self.invpath=os.path.dirname(self.iterfn)
        
        #get meshfile if none is provides assuming the mesh file is named
        #with mesh
        if self.invpath!=None:
            self.meshfn=os.path.join(self.invpath,'MESH')
            if os.path.isfile(self.meshfn)==False:
                for ff in os.listdir(self.invpath):
                    if ff.lower().find('mesh')>=0:
                        self.meshfn=os.path.join(self.invpath,ff)
                if os.path.isfile(self.meshfn)==False:
                    raise NameError('Could not find a mesh file, '+\
                                    'input manually')
            
        #get inmodelfile if none is provides assuming the mesh file is 
        #named with inmodel
        if inmodelfn==None:
            self.inmodelfn=os.path.join(self.invpath,'INMODEL')
            if os.path.isfile(self.inmodelfn)==False:
                for ff in os.listdir(self.invpath):
                    if ff.lower().find('inmodel')>=0:
                        self.inmodelfn=os.path.join(self.invpath,ff)
                if os.path.isfile(self.inmodelfn)==False:
                    raise NameError('Could not find a model file, '+\
                                    'input manually')
        
    def read2DIter(self):
        """
        read2DIter will read an iteration file and combine that info from the 
        data_fn and return a dictionary of variables.
        
        Arguments:
        ----------
            **iterfn** : string
                        full path to iteration file if iterpath=None.  If 
                        iterpath is input then iterfn is just the name
                        of the file without the full path.

        Returns:
        --------
            **Occam2DModel.idict** : dictionary of parameters, 
                                     keys are verbatim from the file, 
                                     except for the key 'model' which is the 
        
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> ocm.read2DIter()
            
        """
    
        #check to see if the file exists
        if os.path.exists(self.iterfn)==False:
            raise IOError('File: '+self.iterfn+' does not exist, check path')
    
        #open file, read lines, close file
        ifid=file(self.iterfn,'r')
        ilines=ifid.readlines()
        ifid.close()
        
        #create dictionary to put things
        self.idict={}
        ii=0
        #put header info into dictionary with similar keys
        while ilines[ii].lower().find('param')!=0:
            iline=ilines[ii].strip().split(':')
            self.idict[iline[0].lower()]=iline[1].strip()
            ii+=1
        
        #get number of parameters
        iline=ilines[ii].strip().split(':')
        nparam=int(iline[1].strip())
        self.idict[iline[0]]=nparam
        self.idict['model']=np.zeros(nparam)
        kk=int(ii+1)
        
        jj=0
        while jj<len(ilines)-kk:
            iline=ilines[jj+kk].strip().split()
            for ll in range(4):
                try:
                    self.idict['model'][jj*4+ll]=float(iline[ll])
                except IndexError:
                    pass
            jj+=1
        
        #get the data file name from the iteration header
        self.data_fn=self.idict['data file']
        if self.data_fn.find(os.sep)==-1:
            self.data_fn=os.path.join(self.invpath,self.data_fn)
        if os.path.isfile(self.data_fn)==False:
            for ff in os.listdir(self.invpath):
                if ff.lower().find('.dat')>=0:
                    self.data_fn=os.path.join(self.invpath,ff)
            if os.path.isfile(self.data_fn)==False:
                raise NameError('Could not find a data file, input manually')
    
    def read2DInmodel(self):
        """
        read an INMODEL file for occam 2D
              
        Arguments:
        ----------
            **inmodelfn** : string
                            full path to INMODEL file
        

        Returns:
        --------
            **Occam2DModel.rows** : list of combined data blocks where first 
                                    number of each list represents the number 
                                    of combined mesh layers for this 
                                    regularization block.  The second number is
                                    the number of columns in the regularization
                                    block layer.
                                    
            **Occam2DModel.cols** : list of combined mesh columns for the 
                                    regularization layer. The sum of this list 
                                    must be equal to the number of mesh columns
                                    
            **Occam2DModel.headerdict** : dictionary of all the header 
                                          information including the binding 
                                          offset.
        
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> ocm.read2DInmodel()
        """
        
        ifid=open(self.inmodelfn,'r')
        
        headerdict={}
        rows=[]
        cols=[]    
        ncols=[]
        
        ilines=ifid.readlines()
        
        for ii,iline in enumerate(ilines):
            if iline.find(':')>0:
                iline=iline.strip().split(':')
                headerdict[iline[0].lower()]=iline[1]
                #append the last line
                if iline[0].lower().find('exception')>0:
                    cols.append(ncols)
            else:
                iline=iline.strip().split()
                iline=[int(jj) for jj in iline]
                if len(iline)==2:
                    if len(ncols)>0:
                        cols.append(ncols)
                    rows.append(iline)
                    ncols=[]
                elif len(iline)>2:
                    ncols=ncols+iline
                    
        self.rows=np.array(rows)
        self.cols=cols
        self.inmodel_headerdict=headerdict
        
    def read2DMesh(self):
        """
        reads an Occam 2D mesh file
        
        Arguments:
        ----------
            **Occam2DModel.meshfn** : string 
                                      full path to mesh file
    
        Returns:
        --------
            **Occam2DModel.hnodes**: array of horizontal nodes 
                                    (column locations (m))
                                    
            **Occam2DModel.vnodes** : array of vertical nodes 
                                      (row locations(m))
                                      
            **Occam2DModel.mdata** : np.array of free parameters
            
        To do:
        ------
            incorporate fixed values
            
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> ocm.read2DMesh()
        """
        
        mfid=file(self.meshfn,'r')
        
        mlines=mfid.readlines()
        
        nh=int(mlines[1].strip().split()[1])-1
        nv=int(mlines[1].strip().split()[2])-1
        
        hnodes=np.zeros(nh)
        vnodes=np.zeros(nv)
        mdata=np.zeros((nh,nv,4),dtype=str)    
        
        #get horizontal nodes
        jj=2
        ii=0
        while ii<nh:
            hline=mlines[jj].strip().split()
            for mm in hline:
                hnodes[ii]=float(mm)
                ii+=1
            jj+=1
        
        #get vertical nodes
        ii=0
        while ii<nv:
            vline=mlines[jj].strip().split()
            for mm in vline:
                vnodes[ii]=float(mm)
                ii+=1
            jj+=1    
        
        #get free parameters        
        for ii,mm in enumerate(mlines[jj+1:]):
            kk=0
            while kk<4:        
                mline=mm.rstrip()
                if mline.lower().find('exception')>0:
                    break
                for jj in range(nh):
                    try:
                        mdata[jj,ii,kk]=mline[jj]
                    except IndexError:
                        pass
                kk+=1
        
        #make the node information an attributes of the occamodel_model class 
        self.hnodes=hnodes
        self.vnodes=vnodes
        self.meshdata=mdata
        
    def get2DData(self):
        """
        get data from data file using the inherited :func:'read2DdataFile' 
        """
        try:
            self.read2DdataFile()
        except AttributeError:
            print 'No Data file defined'
        
    def get2DModel(self):
        """
        get2DModel will create an array based on the FE mesh and fill the 
        values found from the regularization grid.  This way the array can 
        be manipulated as a 2D object and plotted as an image or a mesh.
        
        Returns:
        --------
        
            **Occam2DModel.resmodel** : model array with log resistivity values
            
            **Occam2DModel.plotx** : np.array
                                    horizontal distance of FE mesh (m) blocks
                                    
            **Occam2DModel.ploty** : np.array
                                    depth of vertical nodes of FE mesh (m)
        """
        
        #read iteration file to get model and data file
        self.read2DIter() 
        
        #read in data file as an OccamData type
        print 'Reading data from: ',self.data_fn
        self.get2DData()
        
        #read in MESH file
        print 'Reading mesh from: ',self.meshfn
        self.read2DMesh()
        
        #read in INMODEL
        print 'Reading model from: ',self.inmodelfn
        self.read2DInmodel()
        #get the binding offset which is the right side of the furthest left
        #block, this helps locate the model in relative space
        bndgoff=float(self.inmodel_headerdict['binding offset'])
        
        #make sure that the number of rows and number of columns are the same
        assert len(self.rows)==len(self.cols)
        
        #initiate the resistivity model to the shape of the FE mesh
        resmodel=np.zeros((self.vnodes.shape[0],self.hnodes.shape[0]))
        
        #read in the model and set the regularization block values to map onto
        #the FE mesh so that the model can be plotted as an image or regular 
        #mesh.
        mm=0
        for ii in range(len(self.rows)):
            #get the number of layers to combine
            #this index will be the first index in the vertical direction
            ny1=self.rows[:ii,0].sum()
            #the second index  in the vertical direction
            ny2=ny1+self.rows[ii][0]
            #make the list of amalgamated columns an array for ease
            lc=np.array(self.cols[ii])
            #loop over the number of amalgamated blocks
            for jj in range(len(self.cols[ii])):
                #get first in index in the horizontal direction
                nx1=lc[:jj].sum()
                #get second index in horizontal direction
                nx2=nx1+lc[jj]
                #put the apporpriate resistivity value into all the amalgamated 
                #model blocks of the regularization grid into the forward model
                #grid
                resmodel[ny1:ny2,nx1:nx2]=self.idict['model'][mm]
                mm+=1
        
        #make some arrays for plotting the model
        plotx=np.array([self.hnodes[:ii+1].sum() 
                        for ii in range(len(self.hnodes))])
        ploty=np.array([self.vnodes[:ii+1].sum() 
                        for ii in range(len(self.vnodes))])
        
        #center the grid onto the station coordinates
        x0=bndgoff-plotx[self.cols[0][0]]
        plotx=plotx+x0
        
        #flip the arrays around for plotting purposes
        #plotx=plotx[::-1] and make the first layer start at zero
        ploty=ploty[::-1]-ploty[0]
        
        #make a mesh grid to plot in the model coordinates
        self.meshx,self.meshy=np.meshgrid(plotx,ploty)
        
        #flip the resmodel upside down so that the top is the stations
        resmodel=np.flipud(resmodel)
        
        #make attributes of the class
        self.resmodel=resmodel
        self.plotx=plotx
        self.ploty=ploty
        
        #set the offsets of the stations and station list.
        self.offsetlst=[]
        for rpdict in self.rp_lst:
            self.offsetlst.append(rpdict['offset'])
        
    def plot2DModel(self,data_fn=None,
                    xpad=1.0,ypad=1.0,spad=1.0,ms=10,stationid=None,
                    fdict={'size':8,'rotation':60,'weight':'normal'},
                    dpi=300,ylimits=None,xminorticks=5,yminorticks=1,
                    climits=(0,4), cmap='jet_r',fs=8,femesh='off',
                    regmesh='off',aspect='auto',title='on',meshnum='off',
                    blocknum='off',blkfdict={'size':3},fignum=1,
                    plotdimensions=(10,10),grid='off',yscale='km',
                    xlimits=None):
        """
        plotModel will plot the model output by occam in the iteration file.
        
        Arguments:
        ----------
            
            **data_fn** : string 
                        full path to data file.  If none is input it will use 
                        the data file found in the iteration file.
            
            **xpad** : float (units of km) 
                       padding in the horizontal direction of model.
            
            **ypad** : float (units according to **yscale**)
                       padding in the vertical direction of the top of the 
                       model to fit the station names and markers.
                       
            **spad** : float (units according to **yscale**)
                       padding of station names away from the top of the model,
                       this is kind of awkward at the moment especially if you
                       zoom into the model, it usually looks retarded and 
                       doesn't fit
                    
            **ms** : float  
                     marker size in ambiguous points
            
            **stationid** : tuple (min,max)
                           index of station names to plot -> ex. pb01sdr would 
                           be stationid=(0,4) to plot pb01
                        
            **fdict** : font dictionary for the station names, can have keys:
                
                    *'size'* : font size
                    
                    *'rotation'* : angle of rotation (deg) of font
                    
                    *'weight'* : weight of font 
                    
                    *'color'* : color of font
                    
                    *'style'* : style of font ex. 'italics'
                    
            **plotdimensions** : tuple (x,y)
                               x-y dimensions of the figure in inches
                    
            **dpi** : int 
                      dot per inch of figure, should be 300 for publications
            
            **ylimits** : tuple (min,max)
                          limits of depth scale (km). ex, ylimits=(0,30)
                          
            **xlimits** : tuple (min,max)
                          limits of horizontal scale (km). ex, ylimits=(0,30)
                          
            
            **xminorticks** : int or float
                              location of minor tick marks for the horizontal 
                              axis
            
            **yminorticks** : int or float
                              location of minor tick marks for vertical axis
            
            **climits** : tuple (min,max)
                          limits of log10(resistivity). ex. climits=(0,4)
            
            **cmap** : string
                       color map to plot the model image
                       see matplotlib.cm for all options
            
            **fs** : float
                     font size of axis labels
            
            **femesh** : string ('on','off')
                        'on' to plot finite element forward modeling mesh 
                        (black)
            
            **regmesh** : string ('on','off')
                         'on' to plot regularization mesh (blue)
            
            **aspect** : tuple (width,height)
                        aspect ratio of the figure, depends on your line 
                        length and the depth you want to investigate
            
            **title** : string ('on,'off',input,None)
                        
                        * 'on' to put the RMS and Roughness as the title
                        * input a string that will be added to the RMS and
                          roughness put 
                        * None to not put a title on the plot and print out RMS
                          and roughness
            
            **meshnum** : string ('on','off')
                         'on' to plot FE mesh block numbers
            
            **fignum** : int
                        figure number to plot to
            
            **blocknum** : tuple ('on','off')
                          'on' to plot numbers on the regularization blocks
            
            **blkfdict** : font dictionary for the numbering of regularization 
                           blocks with keys:
                               
                               *'size'* : float font size
                               
                               *'weight'* : font weight
            
            **grid** : string ('major','minor','both')
                        * 'major' for major ticks grid
                        * 'minor' for a grid of the minor ticks
                        * 'both' for a grid with major and minor ticks
            
            **yscale** : string ('km','m')
                        * 'km' for depth in km 
                        * 'm' for depth in meters
        
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> ocm.plot2DModel(ms=20,ylimits=(0,.350),yscale='m',spad=.10,
            >>>                 ypad=.125,xpad=.025,climits=(0,2.5),
            >>>                 aspect='equal')
        """   
                    
        #set the scale of the plot
        if yscale=='km':
            dfactor=1000.
            pfactor=1.0
        elif yscale=='m':
            dfactor=1.
            pfactor=1000.
        else:
            dfactor=1000.
            pfactor=1.0
        
        #get the model
        self.get2DModel()
        
        #set some figure properties to use the maiximum space 
        plt.rcParams['font.size']=int(dpi/40.)
        plt.rcParams['figure.subplot.left']=.08
        plt.rcParams['figure.subplot.right']=.99
        plt.rcParams['figure.subplot.bottom']=.1
        plt.rcParams['figure.subplot.top']=.92
        plt.rcParams['figure.subplot.wspace']=.01
#        plt.rcParams['text.usetex']=True
        
        #plot the model as a mesh
        fig=plt.figure(fignum,plotdimensions,dpi=dpi)
        plt.clf()
        
        #add a subplot to the figure with the specified aspect ratio
        ax=fig.add_subplot(1,1,1,aspect=aspect)
        
        #plot the model as a pcolormesh so the extents are constrained to 
        #the model coordinates
        ax.pcolormesh(self.meshx/dfactor,self.meshy/dfactor,self.resmodel,
                      cmap=cmap,vmin=climits[0],vmax=climits[1])
        
        #make a colorbar for the resistivity
        cbx=mcb.make_axes(ax,shrink=.8,pad=.01)
        cb=mcb.ColorbarBase(cbx[0],cmap=cmap,norm=Normalize(vmin=climits[0],
                        vmax=climits[1]))
        cb.set_label('Resistivity ($\Omega \cdot$m)',
                     fontdict={'size':fs,'weight':'bold'})
        cb.set_ticks(np.arange(int(climits[0]),int(climits[1])+1))
        cb.set_ticklabels(['10$^{0}$'.format(nn) for nn in 
                            np.arange(int(climits[0]),int(climits[1])+1)])
        
        #set the offsets of the stations and plot the stations
        #need to figure out a way to set the marker at the surface in all
        #views.
        for rpdict in self.rp_lst:
            #plot the station marker
            #plots a V for the station cause when you use scatter the spacing
            #is variable if you change the limits of the y axis, this way it
            #always plots at the surface.
            ax.text(rpdict['offset']/dfactor,self.ploty.min(),'V',
                    horizontalalignment='center',
                    verticalalignment='baseline',
                    fontdict={'size':ms,'weight':'bold','color':'black'})
                    
            #put station id onto station marker
            #if there is a station id index
            if stationid!=None:
                ax.text(rpdict['offset']/dfactor,-spad*pfactor,
                        rpdict['station'][stationid[0]:stationid[1]],
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict=fdict)
            #otherwise put on the full station name found form data file
            else:
                ax.text(rpdict['offset']/dfactor,-spad*pfactor,
                        rpdict['station'],
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict=fdict)
        
        #set the initial limits of the plot to be square about the profile line  
        if ylimits==None:  
            ax.set_ylim(abs(max(self.offsetlst)-min(self.offsetlst))/dfactor,
                        -ypad*pfactor)
        else:
            ax.set_ylim(ylimits[1]*pfactor,(ylimits[0]-ypad)*pfactor)
        ax.set_xlim(min(self.offsetlst)/dfactor-(xpad*pfactor),
                     (max(self.offsetlst)/dfactor+(xpad*pfactor)))
        #set the axis properties
        ax.xaxis.set_minor_locator(MultipleLocator(xminorticks*pfactor))
        ax.yaxis.set_minor_locator(MultipleLocator(yminorticks*pfactor))
        if yscale=='km':
            ax.set_xlabel('Horizontal Distance (km)',
                          fontdict={'size':fs,'weight':'bold'})
            ax.set_ylabel('Depth (km)',fontdict={'size':fs,'weight':'bold'})
        elif yscale=='m':
            ax.set_xlabel('Horizontal Distance (m)',
                          fontdict={'size':fs,'weight':'bold'})
            ax.set_ylabel('Depth (m)',fontdict={'size':fs,'weight':'bold'})
        
        #put a grid on if one is desired    
        if grid=='major':
            ax.grid(alpha=.3,which='major')
        if grid=='minor':
            ax.grid(alpha=.3,which='minor')
        if grid=='both':
            ax.grid(alpha=.3,which='both')
        else:
            pass
        
        #set title as rms and roughness
        if type(title) is str:
            if title=='on':
                titlestr=os.path.join(os.path.basename(os.path.dirname(self.iterfn)),
                                      os.path.basename(self.iterfn))
                ax.set_title(titlestr+\
                            ': RMS {0:.2f}, Roughness={1:.0f}'.format(
                            float(self.idict['misfit value']),
                            float(self.idict['roughness value'])),
                            fontdict={'size':fs+1,'weight':'bold'})
            else:
                ax.set_title(title+'; RMS {0:.2f}, Roughness={1:.0f}'.format(
                         float(self.idict['misfit value']),
                         float(self.idict['roughness value'])),
                         fontdict={'size':fs+1,'weight':'bold'})
        else:
            print 'RMS {0:.2f}, Roughness={1:.0f}'.format(
                         float(self.idict['misfit value']),
                         float(self.idict['roughness value'])) 
        
        #plot forward model mesh    
        if femesh=='on':
            for xx in self.plotx/dfactor:
                ax.plot([xx,xx],[0,self.ploty[0]/dfactor],color='k',lw=.5)
            for yy in self.ploty/dfactor:
                ax.plot([self.plotx[0]/dfactor,self.plotx[-1]/dfactor],
                        [yy,yy],color='k',lw=.5)
        
        #plot the regularization mesh
        if regmesh=='on':
            linelst=[]
            for ii in range(len(self.rows)):
                #get the number of layers to combine
                #this index will be the first index in the vertical direction
                ny1=self.rows[:ii,0].sum()
                #the second index  in the vertical direction
                ny2=ny1+self.rows[ii][0]
                #make the list of amalgamated columns an array for ease
                lc=np.array(self.cols[ii])
                yline=ax.plot([self.plotx[0]/dfactor,self.plotx[-1]/dfactor],
                              [self.ploty[-ny1]/dfactor,
                               self.ploty[-ny1]/dfactor],
                              color='b',lw=.5)
                linelst.append(yline)
                #loop over the number of amalgamated blocks
                for jj in range(len(self.cols[ii])):
                    #get first in index in the horizontal direction
                    nx1=lc[:jj].sum()
                    #get second index in horizontal direction
                    nx2=nx1+lc[jj]
                    try:
                        if ny1==0:
                            ny1=1
                        xline=ax.plot([self.plotx[nx1]/dfactor,
                                       self.plotx[nx1]/dfactor],
                                      [self.ploty[-ny1]/dfactor,
                                       self.ploty[-ny2]/dfactor],
                                      color='b',lw=.5)
                        linelst.append(xline)
                    except IndexError:
                        pass
                    
        ##plot the mesh block numbers
        if meshnum=='on':
            kk=1
            for yy in self.ploty[::-1]/dfactor:
                for xx in self.plotx/dfactor:
                    ax.text(xx,yy,'{0}'.format(kk),fontdict={'size':3})
                    kk+=1
                    
        ##plot regularization block numbers
        if blocknum=='on':
            kk=1
            for ii in range(len(self.rows)):
                #get the number of layers to combine
                #this index will be the first index in the vertical direction
                ny1=self.rows[:ii,0].sum()
                #the second index  in the vertical direction
                ny2=ny1+self.rows[ii][0]
                #make the list of amalgamated columns an array for ease
                lc=np.array(self.cols[ii])
                #loop over the number of amalgamated blocks
                for jj in range(len(self.cols[ii])):
                    #get first in index in the horizontal direction
                    nx1=lc[:jj].sum()
                    #get second index in horizontal direction
                    nx2=nx1+lc[jj]
                    try:
                        if ny1==0:
                            ny1=1
                        #get center points of the blocks
                        yy=self.ploty[-ny1]-(self.ploty[-ny1]-
                                                self.ploty[-ny2])/2
                        xx=self.plotx[nx1]-(self.plotx[nx1]-self.plotx[nx2])/2
                        #put the number
                        ax.text(xx/dfactor,yy/dfactor,'{0}'.format(kk),
                                fontdict=blkfdict,
                                horizontalalignment='center',
                                verticalalignment='center')
                        kk+=1
                    except IndexError:
                        pass
                    
        plt.show()
    
    def plotL2Curve(self,fnstem=None,fignum=1,dpi=300):
        """
        PlotL2Curve will plot the RMS vs iteration number for the given 
        inversion folder and roughness vs iteration number
        
        Arguments:
        ----------
            **fnstem** : string
                         filename stem to look for in case multiple inversions 
                         were run in the same folder.  If none then searches 
                         for anything ending in .iter
            
            **fignum** : int
                         figure number

            dpi** : int
                    dots per inch resolution of the figure
                    
        
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> ocm.plotL2Curve(fignum=2)
        """ 

        invpath=os.path.dirname(self.iterfn)        
        
        if fnstem==None:
            iterlst=[os.path.join(invpath,itfile) 
                    for itfile in os.listdir(invpath) if itfile.find('.iter')>0]
        else:
            iterlst=[os.path.join(invpath,itfile) 
                    for itfile in os.listdir(invpath) if itfile.find('.iter')>0 and
                    itfile.find(fnstem)>0]
                    
        nr=len(iterlst)
        
        rmsarr=np.zeros((nr,2))
        
        for itfile in iterlst:
            self.iterfn=itfile
            self.read2DIter()
            ii=int(self.idict['iteration'])
            rmsarr[ii,0]=float(self.idict['misfit value'])
            rmsarr[ii,1]=float(self.idict['roughness value'])
        
        #set the dimesions of the figure
        plt.rcParams['font.size']=int(dpi/40.)
        plt.rcParams['figure.subplot.left']=.08
        plt.rcParams['figure.subplot.right']=.90
        plt.rcParams['figure.subplot.bottom']=.1
        plt.rcParams['figure.subplot.top']=.90
        plt.rcParams['figure.subplot.wspace']=.01
        
        #make figure instance
        fig=plt.figure(fignum,[6,5],dpi=dpi)
        plt.clf()
        
        #make a subplot for RMS vs Iteration
        ax1=fig.add_subplot(1,1,1)
        
        #plot the rms vs iteration
        l1,=ax1.plot(np.arange(1,nr,1),rmsarr[1:,0],'-k',lw=1,marker='d',ms=5)
        
        #plot the median of the RMS
        m1,=ax1.plot(np.arange(0,nr,1),np.repeat(np.median(rmsarr[1:,0]),nr),
                     '--r',lw=.75)
        
        #plot the mean of the RMS
        m2,=ax1.plot(np.arange(0,nr,1),np.repeat(np.mean(rmsarr[1:,0]),nr),
                     ls='--',color='orange',lw=.75)
    
        #make subplot for RMS vs Roughness Plot
        ax2=ax1.twiny()
        
        #plot the rms vs roughness 
        l2,=ax2.plot(rmsarr[1:,1],rmsarr[1:,0],'--b',lw=.75,marker='o',ms=7,
                     mfc='white')
        for ii,rms in enumerate(rmsarr[1:,0],1):
            ax2.text(rmsarr[ii,1],rms,'{0}'.format(ii),
                     horizontalalignment='center',
                     verticalalignment='center',
                     fontdict={'size':6,'weight':'bold','color':'blue'})
        
        #make a legend
        ax1.legend([l1,l2,m1,m2],['RMS','Roughness',
                   'Median_RMS={0:.2f}'.format(np.median(rmsarr[1:,0])),
                    'Mean_RMS={0:.2f}'.format(np.mean(rmsarr[1:,0]))],
                    ncol=4,loc='upper center',columnspacing=.25,markerscale=.75,
                    handletextpad=.15)
                    
        #set the axis properties for RMS vs iteration
        ax1.yaxis.set_minor_locator(MultipleLocator(.1))
        ax1.xaxis.set_minor_locator(MultipleLocator(1))
        ax1.set_ylabel('RMS',fontdict={'size':8,'weight':'bold'})                                   
        ax1.set_xlabel('Iteration',fontdict={'size':8,'weight':'bold'})
        ax1.grid(alpha=.25,which='both')
        ax2.set_xlabel('Roughness',fontdict={'size':8,'weight':'bold',
                                             'color':'blue'})
        for t2 in ax2.get_xticklabels():
            t2.set_color('blue')          
        
        plt.show()
                
    def plotDepthModel(self,dpi=300,depthmm=(1,10000),plottype='1',
                       yscale='log',plotdimensions=(3,6),plot_num=1,fignum=1):
        """
        Plots a depth section profile for a given set of stations.
        
        Arguments:
        ----------
            
            **dpi** : int
                      dots-per-inch resolution of figure
            
            **depthmm** : tuple (min,max)
                          minimum and maximum depth to plot in meters
            
            **plottype** : input as:
                            * '1' to plot all stations found
                            * [station list] to plot only a few stations
            
            **plot_num** : input as:
                          * 1 to plot in different figures
                          * 'all' to plot in all into one figure.
            
            **yscale** : 'log' for logarithmic or 'linear' for linear
            
            **plotdimensions** : tuple (width,height)
                                 figure dimensions in inches
            
            **fignum** : int
                         figure number
                         
        :Example: ::
            
            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> #plot just a few stations depth profile in one figure
            >>> ocm.plotDepthModel(plottype=['MT01','MT05'],plot_num='all')
                         
            
        """

        try:
            self.offsetlst
        except AttributeError:
            self.get2DModel()
        #get stations to plot
        if plottype=='1':
            pstation_lst=np.arange(len(self.station_lst))
        else:
            pstation_lst=[]
            if type(plottype) is not list:
                plottype=[plottype]
            for ps in plottype:
                for ii,ss in enumerate(self.station_lst):
                    if ss.find(ps)==0:
                        pstation_lst.append(ii)
                                  
        #get the average x-spacing within the station region, occam pads by 
        #7 cells by default        
        xavg=np.floor(np.mean([abs(self.plotx[ii]-self.plotx[ii+1]) 
                        for ii in range(7,len(self.plotx)-7)]))
        
        #get the station indices to extract from the model
        slst=[]
        for ff in pstation_lst:
            offset=self.offsetlst[ff]
            for ii,xx in enumerate(self.plotx):
                if offset>=xx-xavg/2. and offset<=xx+xavg/2.:
                    slst.append(ii)
        
        #get depth limits
        if depthmm==None:
            depthmm=(self.ploty.min(),self.ploty.max())
        if depthmm[0]==0:
            depthmm[0]=1
        
        #set the dimesions of the figure
        plt.rcParams['font.size']=int(dpi/40.)
        plt.rcParams['figure.subplot.left']=.15
        plt.rcParams['figure.subplot.right']=.95
        plt.rcParams['figure.subplot.bottom']=.15
        plt.rcParams['figure.subplot.top']=.90
        plt.rcParams['figure.subplot.wspace']=.05
        
        if plot_num=='all':
            #set the dimesions of the figure
            plt.rcParams['font.size']=int(dpi/60.)
            plt.rcParams['figure.subplot.left']=.09
            plt.rcParams['figure.subplot.right']=.95
            plt.rcParams['figure.subplot.bottom']=.15
            plt.rcParams['figure.subplot.top']=.90
            plt.rcParams['figure.subplot.wspace']=.1
            
            fig=plt.figure(fignum,plotdimensions,dpi=dpi)
            plt.clf()
            ns=len(slst)
            #plot the depth section for each station        
            for ii,ss in enumerate(slst):
                ax=fig.add_subplot(1,ns,ii+1)
                
                #plot resistivity vs depth
                if yscale=='linear':
                    p1,=ax.semilogx(10**self.resmodel[:,ss],self.ploty,
                                    ls='steps-')
                elif yscale=='log':
                    if self.ploty[-1]==0.0:
                        self.ploty[-1]=1
                    p1,=ax.loglog(10**self.resmodel[:,ss],self.ploty,
                                  ls='steps-')
                ax.set_ylim(depthmm[1],depthmm[0])
                
                ax.set_title(self.data.station_lst[pstation_lst[ii]],
                             fontdict={'size':10,'weight':'bold'})
                if ii==0:
                    ax.set_ylabel('Depth (m)',
                                  fontdict={'size':8,'weight':'bold'})
                else:
                    plt.setp(ax.yaxis.get_ticklabels(),visible=False)
                if ii==np.round(ns/2.):
                    ax.set_xlabel('Resistivity ($\Omega \cdot$m)',
                                  fontdict={'size':8,'weight':'bold'})
                ax.grid(True,alpha=.3,which='both')
                ax.set_xlim(10**self.resmodel.min(),10**self.resmodel.max())
        else:
            #plot the depth section for each station        
            for ii,ss in enumerate(slst):
                fig=plt.figure(ii+1,plotdimensions,dpi=dpi)
                plt.clf()
                ax=fig.add_subplot(1,1,1)
                
                #plot resistivity vs depth
                if yscale=='linear':
                    p1,=ax.semilogx(10**self.resmodel[:,ss],self.ploty,
                                    ls='steps-')
                elif yscale=='log':
                    if self.ploty[-1]==0.0:
                        self.ploty[-1]=1
                    p1,=ax.loglog(10**self.resmodel[:,ss],self.ploty,
                                  ls='steps-')
                ax.set_ylim(depthmm[1],depthmm[0])
                
                ax.set_title(self.station_lst[pstation_lst[ii]],
                             fontdict={'size':10,'weight':'bold'})    
                ax.set_ylabel('Depth (m)',fontdict={'size':8,'weight':'bold'})
                ax.set_xlabel('Resistivity ($\Omega \cdot$m)',
                              fontdict={'size':8,'weight':'bold'})
                ax.grid(True,alpha=.3,which='both')       
            
class PlotOccam2DResponse():
    """
    class to deal with plotting the occam response
    
    """
    
    def __init__(self, rp_lst, period, **kwargs):
        self.rp_lst = rp_lst
        self.period = period
        self.wl_fn = kwargs.pop('wl_fn', None)

        self.color_mode = kwargs.pop('color_mode', 'color')
        
        self.ms = kwargs.pop('ms', 1.5)
        self.lw = kwargs.pop('lw', .5)
        self.e_capthick = kwargs.pop('e_capthick', .5)
        self.e_capsize = kwargs.pop('e_capsize', 2)

        #color mode
        if self.color_mode == 'color':
            #color for data
            self.cted = kwargs.pop('cted', (0, 0, 1))
            self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')
            
            #color for occam model
            self.ctem = kwargs.pop('ctem', (0, .6, .3))
            self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', '+')
            
            #color for Winglink model
            self.ctewl = kwargs.pop('ctewl', (0, .6, .8))
            self.ctmwl = kwargs.pop('ctmwl', (.8, .7, 0))
            self.mtewl = kwargs.pop('mtewl', 'x')
            self.mtmwl = kwargs.pop('mtmwl', 'x')
         
        #black and white mode
        elif self.color_mode == 'bw':
            #color for data
            self.cted = kwargs.pop('cted', (0, 0, 0))
            self.ctmd = kwargs.pop('ctmd', (0, 0, 0))
            self.mted = kwargs.pop('mted', '*')
            self.mtmd = kwargs.pop('mtmd', 'v')
            
            #color for occam model
            self.ctem = kwargs.pop('ctem', (0.6, 0.6, 0.6))
            self.ctmm = kwargs.pop('ctmm', (0.6, 0.6, 0.6))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', 'x')
            
            #color for Winglink model
            self.ctewl = kwargs.pop('ctewl', (0.3, 0.3, 0.3))
            self.ctmwl = kwargs.pop('ctmwl', (0.3, 0.3, 0.3))
            self.mtewl = kwargs.pop('mtewl', '|')
            self.mtmwl = kwargs.pop('mtmwl', '_')
            
        self.phase_limits = kwargs.pop('phase_limits', (-5, 95))
        self.res_limits = kwargs.pop('res_limits', None)

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        
        self.subplot_wspace = .1
        self.subplot_hspace = .15
        self.subplot_right = .98
        self.subplot_left = .085
        self.subplot_top = .93
        self.subplot_bottom = .1

        self.font_size = kwargs.pop('font_size', 6)
        
        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_num = kwargs.pop('plot_num', 2)
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        
        self.fig_lst = []
        
        if self.plot_yn == 'y':
            self.plot()
        
    def plot(self):
        """
        plot the data and model response if given in individual plots.
         
        """
        
        #make a local copy of the rp_lst    
        rp_lst = list(self.rp_lst)
        
        #create station lst
        self.station_lst = [rp['station'] for rp in rp_lst]
        
        #boolean for adding winglink output to the plots 0 for no, 1 for yes
        addwl = 0
        #read in winglink data file
        if self.wl_fn != None:
            addwl = 1
            self.subplot_hspace+.1
            wld, wlrp_lst, wlplst, wlslst, wltlst = wlt.readOutputFile(
                                                                   self.wl_fn)
            sdict = dict([(ostation, wlstation) for wlstation in wlslst 
                          for ostation in self.station_lst 
                          if wlstation.find(ostation)>=0])
        
        #set a local parameter period for less typing
        period = self.period
                                          
        #---------------plot each respones in a different figure---------------
        if self.plot_type == '1':
            pstation_lst = self.station_lst

        else:
            if type(self.plot_type) is not list:
                self.plot_type = [self.plot_type]
            
            pstation_lst = []
            for ii, station in enumerate(self.station_lst):
                for pstation in self.plot_type:
                    if station.find(pstation)>=0:
                        pstation_lst.append(ii)  
            
        #set the grid of subplots
        gs = gridspec.GridSpec(6, 2,
                               wspace=self.subplot_wspace,
                               left=self.subplot_left,
                               top=self.subplot_top,
                               bottom=self.subplot_bottom, 
                               right=self.subplot_right, 
                               hspace=self.subplot_hspace)
        
        #--> set default font size                           
        plt.rcParams['font.size'] = self.font_size
                                 
        #loop over each station to plot
        for ii, station in enumerate(pstation_lst):
            
            #empty lists for legend marker and label
            rlstte = []
            llstte = []
            rlsttm = []
            llsttm = []
            
            #calculate rms's
            rmslstte = np.hstack((rp_lst[ii]['resxy'][3],
                                  rp_lst[ii]['phasexy'][3]))
            rmslsttm = np.hstack((rp_lst[ii]['resyx'][3],
                                  rp_lst[ii]['phaseyx'][3]))
            rmste = np.sqrt(np.sum([rms**2 for rms in rmslstte])/
                            len(rmslstte))
            rmstm = np.sqrt(np.sum([rms**2 for rms in rmslsttm])/
                            len(rmslsttm))
            
            fig = plt.figure(ii+1, self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            
            
            #--> set subplot instances
            #---plot both TE and TM in same subplot---
            if self.plot_num == 1:
                axrte = fig.add_subplot(gs[:4,:])
                axrtm = axrte
                axpte = fig.add_subplot(gs[-2:,:],sharex=axrte)
                axptm = axpte
                
            #---plot TE and TM in separate subplots---
            elif self.plot_num == 2:
                axrte = fig.add_subplot(gs[:4,0])
                axrtm = fig.add_subplot(gs[:4,1])
                axpte = fig.add_subplot(gs[-2:,0], sharex=axrte)
                axptm = fig.add_subplot(gs[-2:,1], sharex=axrtm)
            
            #------------Plot Resistivity----------------------------------
            #cut out missing data points first
            #--> data
            rxy = np.where(rp_lst[ii]['resxy'][0]!=0)[0]
            ryx = np.where(rp_lst[ii]['resyx'][0]!=0)[0]
            
            #--> response
            mrxy = np.where(rp_lst[ii]['resxy'][2]!=0)[0]
            mryx = np.where(rp_lst[ii]['resyx'][2]!=0)[0]
            
            #--> TE mode Data 
            if len(rxy) > 0:
                rte = axrte.errorbar(period[rxy],
                                     10**rp_lst[ii]['resxy'][0][rxy],
                                     ls=':',
                                     marker=self.mted,
                                     ms=self.ms,
                                     mfc=self.cted,
                                     mec=self.cted,
                                     color=self.cted,
                                     yerr=np.log(10)*\
                                         rp_lst[ii]['resxy'][1][rxy]*\
                                         10**rp_lst[ii]['resxy'][0][rxy],
                                    ecolor=self.cted,
                                    picker=2,
                                    lw=self.lw,
                                    elinewidth=self.lw,
                                    capsize=self.e_capsize,
                                    capthick=self.e_capthick)
                rlstte.append(rte[0])
                llstte.append('$Obs_{TE}$')
            else:
                pass
            
            #--> TE mode Model Response
            if len(mrxy) > 0:
                yerrxy = 10**(rp_lst[ii]['resxy'][3][mrxy]*\
                         rp_lst[ii]['resxy'][2][mrxy]/np.log(10))
                r3 = axrte.errorbar(period[mrxy],
                                    10**rp_lst[ii]['resxy'][2][mrxy],
                                    ls='--',
                                    marker=self.mtem,
                                    ms=self.ms,
                                    mfc=self.ctem,
                                    mec=self.ctem,
                                    color=self.ctem,
                                    yerr=yerrxy,
                                    ecolor=self.ctem,
                                    lw=self.lw,
                                    elinewidth=self.lw,
                                    capsize=self.e_capsize,
                                    capthick=self.e_capthick)
                rlstte.append(r3[0])
                llstte.append('$Mod_{TE}$')
            else:
                pass
            
            #--> TM mode data
            if len(ryx) > 0:
                rtm = axrtm.errorbar(period[ryx],
                                     10**rp_lst[ii]['resyx'][0][ryx],
                                     ls=':',
                                     marker=self.mtmd,
                                     ms=self.ms,
                                     mfc=self.ctmd,
                                     mec=self.ctmd,
                                     color=self.ctmd,
                                     yerr=np.log(10)*\
                                          rp_lst[ii]['resyx'][1][ryx]*\
                                          10**rp_lst[ii]['resyx'][0][ryx],
                                     ecolor=self.ctmd,
                                     picker=2,
                                     lw=self.lw,
                                     elinewidth=self.lw,
                                     capsize=self.e_capsize,
                                     capthick=self.e_capthick)
                rlsttm.append(rtm[0])
                llsttm.append('$Obs_{TM}$')
            else:
                pass 

            #--> TM mode model response
            if len(mryx)>0:
                yerryx = 10**(rp_lst[ii]['resyx'][3][mryx]*\
                             rp_lst[ii]['resyx'][2][mryx]/np.log(10))
                r4 = axrtm.errorbar(period[mryx],
                                    10**rp_lst[ii]['resyx'][2][mryx],
                                    ls='--',
                                    marker=self.mtmm,
                                    ms=self.ms,
                                    mfc=self.ctmm,
                                    mec=self.ctmm,
                                    color=self.ctmm,
                                    yerr=yerryx,
                                    ecolor=self.ctmm,
                                    lw=self.lw,
                                    elinewidth=self.lw,
                                    capsize=self.e_capsize,
                                    capthick=self.e_capthick)
                rlsttm.append(r4[0])
                llsttm.append('$Mod_{TM}$')
            else:
                pass

            #--------------------plot phase--------------------------------
            #cut out missing data points first
            #--> data
            pxy = np.where(rp_lst[ii]['phasexy'][0]!=0)[0]
            pyx = np.where(rp_lst[ii]['phaseyx'][0]!=0)[0]
            
            #--> reponse
            mpxy = np.where(rp_lst[ii]['phasexy'][2]!=0)[0]
            mpyx = np.where(rp_lst[ii]['phaseyx'][2]!=0)[0]
            
            #--> TE mode data
            if len(pxy) > 0:
                axpte.errorbar(period[pxy],
                               rp_lst[ii]['phasexy'][0][pxy],
                               ls=':',
                               marker=self.mted,
                               ms=self.ms,
                               mfc=self.cted,
                               mec=self.cted,
                               color=self.cted,
                               yerr=rp_lst[ii]['phasexy'][1][pxy],
                               ecolor=self.cted,
                               picker=1,
                               lw=self.lw,
                               elinewidth=self.lw,
                               capsize=self.e_capsize,
                               capthick=self.e_capthick)
            else:
                pass
            
            #--> TE mode response
            if len(mpxy) > 0:
                axpte.errorbar(period[mpxy],
                               rp_lst[ii]['phasexy'][2][mpxy],
                               ls='--',
                               marker=self.mtem,
                               ms=self.ms,
                               mfc=self.ctem,
                               mec=self.ctem,
                               color=self.ctem,
                               yerr=rp_lst[ii]['phasexy'][3][mpxy],
                               ecolor=self.ctem,
                               lw=self.lw,
                               elinewidth=self.lw,
                               capsize=self.e_capsize,
                               capthick=self.e_capthick)
            else:
                pass
            
            #--> TM mode data
            if len(pyx)>0:
                axptm.errorbar(period[pyx],
                               rp_lst[ii]['phaseyx'][0][pyx],
                               ls=':',
                               marker=self.mtmd,
                               ms=self.ms,
                               mfc=self.ctmd,
                               mec=self.ctmd,
                               color=self.ctmd,
                               yerr=rp_lst[ii]['phaseyx'][1][pyx],
                               ecolor=self.ctmd,
                               picker=1,
                               lw=self.lw,
                               elinewidth=self.lw,
                               capsize=self.e_capsize,
                               capthick=self.e_capthick)
            else:
                pass
            
            #--> TM mode response
            if len(mpyx) > 0:
                axptm.errorbar(period[mpyx],
                               rp_lst[ii]['phaseyx'][2][mpyx],
                               ls='--',
                               marker=self.mtmm,
                               ms=self.ms,
                               mfc=self.ctmm,
                               mec=self.ctmm,
                               color=self.ctmm,
                               yerr=rp_lst[ii]['phaseyx'][3][mpyx],
                               ecolor=self.ctmm,
                               lw=self.lw,
                               elinewidth=self.lw,
                               capsize=self.e_capsize,
                               capthick=self.e_capthick)
            else:
                pass
            
            
            #--------------add in winglink responses------------------------
            if addwl==1:
                try:
                    wlrms=wld[sdict[station]]['rms']
                    axrte.set_title(self.station_lst[ii]+
                                   '\n rms_occ_TE={0:.2f}'.format(rmste)+
                                   'rms_occ_TM={0:.2f}'.format(rmstm)+
                                   'rms_wl={0:.2f}'.format(wlrms),
                                   fontdict={'size':self.font_size,
                                             'weight':'bold'})
                    for ww, wlstation in enumerate(wlslst):
                        if wlstation.find(station)==0:
                            print '{0} was Found {0} in winglink file'.format(
                                                            station, wlstation)
                            wlrpdict = wlrp_lst[ww]
                    
                    zrxy = [np.where(wlrpdict['resxy'][0]!=0)[0]]
                    zryx = [np.where(wlrpdict['resyx'][0]!=0)[0]]
                    
                     #plot winglink resistivity
                    r5 = axrte.loglog(wlplst[zrxy],
                                      wlrpdict['resxy'][1][zrxy],
                                      ls='-.',
                                      marker=self.mtewl,
                                      ms=self.ms,
                                      color=self.ctewl,
                                      mfc=self.ctewl,
                                      lw=self.lw)
                    r6 = axrtm.loglog(wlplst[zryx],
                                      wlrpdict['resyx'][1][zryx],
                                      ls='-.',
                                      marker=self.mtmwl,
                                      ms=self.ms,
                                      color=self.ctmwl,
                                      mfc=self.ctmwl,
                                      lw=self.lw)
                    
                    #plot winglink phase
                    axpte.semilogx(wlplst[zrxy],
                                   wlrpdict['phasexy'][1][zrxy],
                                   ls='-.',
                                   marker=self.mtewl,
                                   ms=self.ms,
                                   color=self.ctewl,
                                   mfc=self.ctewl,
                                   lw=self.lw)
                                   
                    axptm.semilogx(wlplst[zryx],
                                   wlrpdict['phaseyx'][1][zryx],
                                   ls='-.',
                                   marker=self.mtmwl,
                                   ms=self.ms,
                                   color=self.ctmwl,
                                   mfc=self.ctmwl,
                                   lw=self.lw)
                    
                    rlstte.append(r5[0])
                    rlsttm.append(r6[0])
                    llstte.append('$WLMod_{TE}$')
                    llsttm.append('$WLMod_{TM}$')
                except IndexError:
                    print 'Station not present'
            else:
                if self.plot_num == 1:
                    axrte.set_title(self.station_lst[ii]+\
                    ' rms_TE={0:.2f}, rms_TM={1:.2f}'.format(rmste,rmstm),
                              fontdict={'size':self.font_size+2,
                                        'weight':'bold'})
                elif self.plot_num == 2:
                    axrte.set_title(self.station_lst[ii]+\
                                    ' rms_TE={0:.2f}'.format(rmste),
                                    fontdict={'size':self.font_size+2,
                                              'weight':'bold'})
                    axrtm.set_title(self.station_lst[ii]+\
                                    ' rms_TM={0:.2f}'.format(rmstm),
                                    fontdict={'size':self.font_size+2,
                                              'weight':'bold'})
            
            #set the axis properties
            for aa, axr in enumerate([axrte, axrtm]):
                #set both axes to logarithmic scale
                axr.set_xscale('log')
                axr.set_yscale('log')
                
                #put on a grid
                axr.grid(True, alpha=.3, which='both')
                axr.yaxis.set_label_coords(-.12, .5)
                
                #set resistivity limits if desired
                if self.res_limits != None:
                    axr.set_ylim(10**self.res_limits[0],
                                 10**self.res_limits[1])
                    
                #set the tick labels to invisible
                plt.setp(axr.xaxis.get_ticklabels(), visible=False)
                if aa == 0:
                    axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                                   fontdict={'size':self.font_size+2,
                                             'weight':'bold'})
                           
                #set legend based on the plot type
                if self.plot_num == 1:
                    if aa == 0:
                        axr.legend(rlstte+rlsttm,llstte+llsttm,
                                   loc=2,markerscale=1,
                                   borderaxespad=.05,
                                   labelspacing=.08,
                                   handletextpad=.15,
                                   borderpad=.05,
                                   prop={'size':self.font_size+1})
                elif self.plot_num == 2:
                    if aa == 0:
                        axr.legend(rlstte,
                                   llstte,
                                   loc=2,markerscale=1,
                                   borderaxespad=.05,
                                   labelspacing=.08,
                                   handletextpad=.15,
                                   borderpad=.05,
                                   prop={'size':self.font_size+1}) 
                                           
                    if aa==1:
                        axr.legend(rlsttm,
                                   llsttm,
                                   loc=2,markerscale=1,
                                   borderaxespad=.05,
                                   labelspacing=.08,
                                   handletextpad=.15,
                                   borderpad=.05,
                                   prop={'size':self.font_size+1})
                      
            #set Properties for the phase axes
            for aa, axp in enumerate([axpte, axptm]):
                #set the x-axis to log scale
                axp.set_xscale('log')
                
                #set the phase limits
                axp.set_ylim(self.phase_limits)
                
                #put a grid on the subplot
                axp.grid(True,alpha=.3,which='both')
                
                #set the tick locations
                axp.yaxis.set_major_locator(MultipleLocator(10))
                axp.yaxis.set_minor_locator(MultipleLocator(2))
                
                #set the x axis label
                axp.set_xlabel('Period (s)',
                               fontdict={'size':self.font_size+2,
                                         'weight':'bold'})
                
                #put the y label on the far left plot
                axp.yaxis.set_label_coords(-.12,.5)
                if aa==0:
                    axp.set_ylabel('Phase (deg)',
                                   fontdict={'size':self.font_size+2,
                                             'weight':'bold'})
            
            #make sure the axis and figure are accessible to the user
            self.fig_lst.append({'station':station, 'fig':fig, 'axrte':axrte,
                                 'axrtm':axrtm, 'axpte':axpte, 'axptm':axptm})
        
        #set the plot to be full screen well at least try
        plt.show()
        
    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        """
        plt.close('all')
        self.plot()
        
    def save_figures(self, save_path, fig_fmt='pdf', fig_dpi=None, 
                     close_fig='y'):
        """
        save all the figure that are in self.fig_lst
        """
        
        if not os.path.exists(save_path):
            os.mkdir(save_path)
            
        for fdict in self.fig_lst:
            svfn = '{0}_resp.{1}'.format(fdict['station'], fig_fmt)
            fdict['fig'].savefig(os.path.join(save_path, svfn), 
                                 dpi=self.fig_dpi)
            if close_fig == 'y':
                plt.close(fdict['fig'])
                
class PlotPseudoSection(object):
    """
    plot a pseudo section of the data and response if given
    
    """
    
    def __init__(self, rp_lst, period, **kwargs):
        
        self.rp_lst = rp_lst
        self.period = period
        self.station_lst = [rp['station'] for rp in self.rp_lst]

        self.plot_resp = kwargs.pop('plot_resp', 'y')
        
        self.label_lst = ['$r_{TE-Data}$','$r_{TE-Model}$',
                          '$r_{TM-Data}$','$r_{TM-Model}$',
                          '$\phi_{TE-Data}$','$\phi_{TE-Model}$',
                          '$\phi_{TM-Data}$','$\phi_{TM-Model}$']
        
        self.phase_limits_te = kwargs.pop('phase_limits_te', (-5, 95))
        self.phase_limits_tm = kwargs.pop('phase_limits_tm', (-5, 95))
        self.res_limits_te = kwargs.pop('res_limits_te', (0,3))
        self.res_limits_tm = kwargs.pop('res_limits_tm', (0,3))
        
        self.phase_cmap = kwargs.pop('phase_cmap', 'jet')
        self.res_cmap = kwargs.pop('res_cmap', 'jet_r')
        
        self.ml = kwargs.pop('ml', 2)
        self.station_id = kwargs.pop('station_id', [0,4])

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        
        self.subplot_wspace = .05
        self.subplot_hspace = .0
        self.subplot_right = .95
        self.subplot_left = .085
        self.subplot_top = .97
        self.subplot_bottom = .1

        self.font_size = kwargs.pop('font_size', 6)
        
        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_num = kwargs.pop('plot_num', 2)
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        
        self.cb_shrink = .7
        self.cb_pad = .015
        
        self.axrte = None
        self.axrtm = None
        self.axpte = None
        self.axptm = None
        self.axmrte = None
        self.axmrtm = None
        self.axmpte = None
        self.axmptm = None
        
        self.fig = None
        
        if self.plot_yn == 'y':
            self.plot()
                        
    def plot(self):
        """
        plot pseudo section of data and response if given
        
        """
        if self.plot_resp == 'y':
            nr = 2
        else:
            nr = 1
            
        ns = len(self.station_lst)
        nf = len(self.period)
        ylimits = (self.period.max(), self.period.min())
        
        #make a grid for pcolormesh so you can have a log scale
        #get things into arrays for plotting
        offset_lst = np.zeros(ns)
        resxy_arr = np.zeros((nf, ns, nr))    
        resyx_arr = np.zeros((nf, ns, nr))    
        phasexy_arr = np.zeros((nf, ns, nr))    
        phaseyx_arr = np.zeros((nf, ns, nr))
    
        for ii, rpdict in enumerate(self.rp_lst):
            offset_lst[ii] = rpdict['offset']     
            resxy_arr[:, ii, 0] = rpdict['resxy'][0]
            resyx_arr[:, ii, 0] = rpdict['resyx'][0]
            phasexy_arr[:, ii, 0] = rpdict['phasexy'][0]
            phaseyx_arr[:, ii, 0] = rpdict['phaseyx'][0]
            if self.plot_resp == 'y':
                resxy_arr[:,ii,1]=rpdict['resxy'][2]
                resyx_arr[:,ii,1]=rpdict['resyx'][2]
                phasexy_arr[:,ii,1]=rpdict['phasexy'][2]
                phaseyx_arr[:,ii,1]=rpdict['phaseyx'][2]
                
                
        #make a meshgrid for plotting
        #flip frequency so bottom corner is long period
        dgrid, fgrid = np.meshgrid(offset_lst, self.period[::-1])
    
        #make list for station labels
        slabel = [self.station_lst[ss][self.station_id[0]:self.station_id[1]] 
                    for ss in range(0, ns, self.ml)]

        xloc = offset_lst[0]+abs(offset_lst[0]-offset_lst[1])/5
        yloc = 1.10*self.period[1]
        
        
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top        
        
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()
           
        if self.plot_resp == 'y':
            
            gs1 = gridspec.GridSpec(1, 2,
                                    left=self.subplot_left,
                                    right=self.subplot_right,
                                    wspace=.125)
        
            gs2 = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                   hspace=self.subplot_hspace,
                                                   wspace=self.subplot_wspace,
                                                   subplot_spec=gs1[0])
            gs3 = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                   hspace=self.subplot_hspace,
                                                   wspace=self.subplot_wspace,
                                                   subplot_spec=gs1[1])
            
            #plot TE resistivity data
            
            self.axrte = plt.Subplot(self.fig, gs2[0, 0])
            self.fig.add_subplot(self.axrte)
            self.axrte.pcolormesh(dgrid, 
                                  fgrid, 
                                  np.flipud(resxy_arr[:, :, 0]),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_te[0],
                                  vmax=self.res_limits_te[1])
            
            #plot TE resistivity model
            self.axmrte = plt.Subplot(self.fig, gs2[0, 1])
            self.fig.add_subplot(self.axmrte)
            self.axmrte.pcolormesh(dgrid,
                                   fgrid,
                                   np.flipud(resxy_arr[:, :, 1]),
                                   cmap=self.res_cmap,
                                   vmin=self.res_limits_te[0],
                                   vmax=self.res_limits_te[1])
            
            #plot TM resistivity data
            self.axrtm = plt.Subplot(self.fig, gs3[0, 0])
            self.fig.add_subplot(self.axrtm)   
            self.axrtm.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(resyx_arr[:,:,0]),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_tm[0],
                                  vmax=self.res_limits_tm[1])
            
            #plot TM resistivity model
            self.axmrtm = plt.Subplot(self.fig, gs3[0, 1])
            self.fig.add_subplot(self.axmrtm)
            self.axmrtm.pcolormesh(dgrid,
                                   fgrid,
                                   np.flipud(resyx_arr[:,:,1]),
                                   cmap=self.res_cmap,
                                   vmin=self.res_limits_tm[0],
                                   vmax=self.res_limits_tm[1])
    
            #plot TE phase data
            self.axpte = plt.Subplot(self.fig, gs2[1, 0])
            self.fig.add_subplot(self.axpte)
            self.axpte.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(phasexy_arr[:,:,0]),
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_te[0],
                                  vmax=self.phase_limits_te[1])
            
            #plot TE phase model
            self.axmpte = plt.Subplot(self.fig, gs2[1, 1])
            self.fig.add_subplot(self.axmpte)
            self.axmpte.pcolormesh(dgrid,
                                   fgrid,
                                   np.flipud(phasexy_arr[:,:,1]),
                                   cmap=self.phase_cmap,
                                   vmin=self.phase_limits_te[0],
                                   vmax=self.phase_limits_te[1])
            
            #plot TM phase data 
            self.axptm = plt.Subplot(self.fig, gs3[1, 0])
            self.fig.add_subplot(self.axptm)              
            self.axptm.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(phaseyx_arr[:,:,0]),
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_tm[0],
                                  vmax=self.phase_limits_tm[1])
            
            #plot TM phase model
            self.axmptm = plt.Subplot(self.fig, gs3[1, 1])
            self.fig.add_subplot(self.axmptm)
            self.axmptm.pcolormesh(dgrid,
                                   fgrid,
                                   np.flipud(phaseyx_arr[:,:,1]),
                                   cmap=self.phase_cmap,
                                   vmin=self.phase_limits_tm[0],
                                   vmax=self.phase_limits_tm[1])
            
            axlst=[self.axrte, self.axmrte, self.axrtm, self.axmrtm, 
                   self.axpte, self.axmpte, self.axptm, self.axmptm]
            
            #make everthing look tidy
            for xx, ax in enumerate(axlst):
                ax.semilogy()
                ax.set_ylim(ylimits)
                ax.xaxis.set_ticks(offset_lst[np.arange(0, ns, self.ml)])
                ax.xaxis.set_ticks(offset_lst, minor=True)
                ax.xaxis.set_ticklabels(slabel)
                ax.set_xlim(offset_lst.min(),offset_lst.max())
                if np.remainder(xx, 2.0) == 1:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cbx = mcb.make_axes(ax, 
                                        shrink=self.cb_shrink, 
                                        pad=self.cb_pad)
                                        
                if xx < 4:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    if xx == 1:
                        cb = mcb.ColorbarBase(cbx[0],cmap=self.res_cmap,
                                norm=Normalize(vmin=self.res_limits_te[0],
                                               vmax=self.res_limits_te[1]))
                    if xx == 3:
                        cb = mcb.ColorbarBase(cbx[0],cmap=self.res_cmap,
                                norm=Normalize(vmin=self.res_limits_tm[0],
                                               vmax=self.res_limits_tm[1]))
                        cb.set_label('App. Res. ($\Omega \cdot$m)',
                                     fontdict={'size':self.font_size+1,
                                               'weight':'bold'})
                else:
                    if xx == 5:
                        cb = mcb.ColorbarBase(cbx[0],cmap=self.phase_cmap,
                                norm=Normalize(vmin=self.phase_limits_te[0],
                                               vmax=self.phase_limits_te[1]))
                    if xx == 7:
                        cb = mcb.ColorbarBase(cbx[0],cmap=self.phase_cmap,
                                norm=Normalize(vmin=self.phase_limits_tm[0],
                                               vmax=self.phase_limits_tm[1]))
                        cb.set_label('Phase (deg)', 
                                     fontdict={'size':self.font_size+1,
                                               'weight':'bold'})
                ax.text(xloc, yloc, self.label_lst[xx],
                        fontdict={'size':self.font_size+1},
                        bbox={'facecolor':'white'},
                        horizontalalignment='left',
                        verticalalignment='top')
                if xx == 0 or xx == 4:
                    ax.set_ylabel('Period (s)',
                                  fontdict={'size':self.font_size+2, 
                                  'weight':'bold'})
                if xx>3:
                    ax.set_xlabel('Station',fontdict={'size':self.font_size+2,
                                                      'weight':'bold'})
                
                    
            plt.show()
            
        else: 
            gs1 = gridspec.GridSpec(2, 2,
                        left=self.subplot_left,
                        right=self.subplot_right,
                        hspace=self.subplot_hspace,
                        wspace=self.subplot_wspace)
            
            #plot TE resistivity data
            self.axrte = self.fig.add_subplot(gs1[0, 0])
            self.axrte.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(resxy_arr[:, :, 0]),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_te[0],
                                  vmax=self.res_limits_te[1])
            
            #plot TM resistivity data               
            self.axrtm = self.fig.add_subplot(gs1[0, 1])
            self.axrtm.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(resyx_arr[:, :, 0]),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_tm[0],
                                  vmax=self.res_limits_tm[1])
            
            #plot TE phase data
            self.axpte = self.fig.add_subplot(gs1[1, 0])
            self.axpte.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(phasexy_arr[:, :, 0]),
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_te[0],
                                  vmax=self.phase_limits_te[1])
            
            #plot TM phase data               
            self.axptm = self.fig.add_subplot(gs1[1, 1])
            self.axptm.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(phaseyx_arr[:,:, 0]),
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_tm[0],
                                  vmax=self.phase_limits_tm[1])
            
            
            axlst=[self.axrte, self.axrtm, self.axpte, self.axptm]
            
            #make everything look tidy
            for xx,ax in enumerate(axlst):
                ax.semilogy()
                ax.set_ylim(ylimits)
                ax.xaxis.set_ticks(offset_lst[np.arange(0, ns, self.ml)])
                ax.xaxis.set_ticks(offset_lst, minor=True)
                ax.xaxis.set_ticklabels(slabel)
                ax.grid(True, alpha=.25)
                ax.set_xlim(offset_lst.min(),offset_lst.max())
                cbx = mcb.make_axes(ax, 
                                    shrink=self.cb_shrink, 
                                    pad=self.cb_pad)
                if xx == 0:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(cbx[0],cmap=self.res_cmap,
                                    norm=Normalize(vmin=self.res_limits_te[0],
                                                   vmax=self.res_limits_te[1]))
                    cb.set_ticks(np.arange(self.res_limits_te[0], 
                                           self.res_limits_te[1]+1))
                elif xx == 1:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    
                    cb = mcb.ColorbarBase(cbx[0],cmap=self.res_cmap,
                                    norm=Normalize(vmin=self.res_limits_tm[0],
                                                   vmax=self.res_limits_tm[1]))
                    cb.set_label('App. Res. ($\Omega \cdot$m)',
                                 fontdict={'size':self.font_size+1,
                                           'weight':'bold'})
                    cb.set_ticks(np.arange(self.res_limits_tm[0], 
                                           self.res_limits_tm[1]+1))
                elif xx == 2:
                    cb = mcb.ColorbarBase(cbx[0],cmap=self.phase_cmap,
                                    norm=Normalize(vmin=self.phase_limits_te[0],
                                                   vmax=self.phase_limits_te[1]))
                    cb.set_ticks(np.arange(self.phase_limits_te[0], 
                                           self.phase_limits_te[1]+1, 15))
                elif xx == 3:
                    cb = mcb.ColorbarBase(cbx[0],cmap=self.phase_cmap,
                                    norm=Normalize(vmin=self.phase_limits_tm[0],
                                                   vmax=self.phase_limits_tm[1]))
                    cb.set_label('Phase (deg)',
                                 fontdict={'size':self.font_size+1,
                                           'weight':'bold'})
                    cb.set_ticks(np.arange(self.phase_limits_te[0], 
                                           self.phase_limits_te[1]+1, 15))
                ax.text(xloc, yloc, self.label_lst[xx],
                        fontdict={'size':self.font_size+1},
                        bbox={'facecolor':'white'},
                        horizontalalignment='left',
                        verticalalignment='top')
                if xx == 0 or xx == 2:
                    ax.set_ylabel('Period (s)',
                                  fontdict={'size':self.font_size+2,
                                            'weight':'bold'})
                if xx>1:
                    ax.set_xlabel('Station',fontdict={'size':self.font_size+2,
                                                      'weight':'bold'})
                
                    
            plt.show()
            
    def redraw_plot(self):
        """
        redraw plot if parameters change
        """
        
        plt.close(self.fig)
        self.plot()
        
    def save_figure(self, save_fn, file_format='pdf', orientation='portrait', 
                  fig_dpi=None, close_plot='y'):
        """
        save_plot will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format
                            
                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path
                            
            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 
                              
            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait
                              
            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.
                          
            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open
                          
        :Example: ::
            
            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')
            
        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi
            
        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation)
            
        else:
            save_fn = os.path.join(save_fn, 'OccamPseudoSection.'+
                                    file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation)
        
        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)
        
        else:
            pass
        
        self.fig_fn = save_fn
        print 'Saved figure to: '+self.fig_fn
        
    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  
        
        Use this if you change an of the .fig or axes properties
        
        :Example: ::
            
            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()
        
        """

        self.fig.canvas.draw()
                          
    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """
        
        return ("Plots a pseudo section of TE and TM modes for data and "
                "response if given.") 

#==============================================================================
# plot all response from a given folder                             
#==============================================================================
                                              
class PlotAllResponses(object):
    """
    plot all responses for all iterations
    
    need to input a nested list
    """
    
    def __init__(self, rp_lst, period, pstation_lst=None, **kwargs):
        self.rp_lst = rp_lst
        self.period = period
        
        if pstation_lst == None:
            self.pstation_lst = [rp['station'] for rp in self.rp_lst[0]]
        else:
            self.pstation_lst = pstation_lst
            
        self.station_lst = [rp['station'] for rp in self.rp_lst[0]]
        self.ms = kwargs.pop('ms', 1.5)
        self.lw = kwargs.pop('lw', .5)
        self.e_capthick = kwargs.pop('e_capthick', .5)
        self.e_capsize = kwargs.pop('e_capsize', 2)
            
        self.phase_limits = kwargs.pop('phase_limits', (-5, 95))
        self.res_limits = kwargs.pop('res_limits', None)

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        
        self.subplot_wspace = .2
        self.subplot_hspace = .15
        self.subplot_right = .98
        self.subplot_left = .085
        self.subplot_top = .93
        self.subplot_bottom = .1
        
        self.axrte = None
        self.axrtm = None
        self.axpte = None
        self.axptm = None

        self.font_size = kwargs.pop('font_size', 6)
        
        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_num = kwargs.pop('plot_num', 2)
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        
        self.fig_lst = []
        
        if self.plot_yn == 'y':
            self.plot()
        
        
        
    def plot(self):
        """
        plot all responses into one figure with subplots for TE and TM
        
        """
        
        gs = gridspec.GridSpec(6, 2, wspace=.20)
        
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        
        nresp = len(self.rp_lst)
        
        #--> make a list of colors to fade for each iteration
        color_lst = [(cc, 0, 1-cc) for cc in np.arange(0,1,1./nresp)]
        
        for ss, station in enumerate(self.pstation_lst, 1):
            fig = plt.figure(self.fig_num+ss, self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            
            axrte = fig.add_subplot(gs[:4, 0])
            axrtm = fig.add_subplot(gs[:4, 1])
            axpte = fig.add_subplot(gs[-2:, 0], sharex=axrte)
            axptm = fig.add_subplot(gs[-2:, 1], sharex=axrtm)
            
            rmstelst = []
            rmstmlst = []
            rmstestr = []
            rmstmstr = []
            #read responses
            for jj, rp in enumerate(self.rp_lst):
                
                ii = np.where(np.array(self.station_lst)==station)[0][0]
                
                rmslstte=np.hstack((rp[ii]['resxy'][3],
                                    rp[ii]['phasexy'][3]))
                rmslsttm=np.hstack((rp[ii]['resyx'][3],
                                    rp[ii]['phaseyx'][3]))
                rmste=np.sqrt(np.sum(ms**2 for ms in rmslstte)/len(rmslstte))
                rmstm=np.sqrt(np.sum(ms**2 for ms in rmslsttm)/len(rmslsttm))
                rmstelst.append('{0} rms={1:.3f}'.format(jj, rmste))
                rmstmlst.append('{0} rms={1:.3f}'.format(jj, rmstm))
                rmstestr.append(rmste)
                rmstmstr.append(rmstm)
                #plot resistivity
                
                
                if jj == 0:
                    #cut out missing data points first
                    rxy = np.where(rp[ii]['resxy'][0]!=0)[0]
                    ryx = np.where(rp[ii]['resyx'][0]!=0)[0]
                    r1, = axrte.loglog(self.period[rxy],
                                       10**rp[ii]['resxy'][0][rxy],
                                       ls=':',
                                       marker='s',
                                       ms=self.ms,
                                       color='k',
                                       mfc='k')
                    r2,=axrtm.loglog(self.period[ryx],
                                     10**rp[ii]['resyx'][0][ryx],
                                     ls=':',
                                     marker='o',
                                     ms=self.ms,
                                     color='k',
                                     mfc='k')
                    rlstte = [r1]
                    rlsttm = [r2]
            
                mrxy = [np.where(rp[ii]['resxy'][2]!=0)[0]]
                mryx = [np.where(rp[ii]['resyx'][2]!=0)[0]]
                r3, = axrte.loglog(self.period[mrxy],
                                   10**rp[ii]['resxy'][2][mrxy],
                                   ls='-',
                                   color=color_lst[jj])
                r4, = axrtm.loglog(self.period[mryx],
                                   10**rp[ii]['resyx'][2][mryx],
                                    ls='-',
                                    color=color_lst[jj])
            
                rlstte.append(r3)
                rlsttm.append(r4)
                                    
                #plot phase
                #cut out missing data points first
                pxy = [np.where(rp[ii]['phasexy'][0]!=0)[0]]
                pyx = [np.where(rp[ii]['phaseyx'][0]!=0)[0]]
                
                if jj == 0:            
                    axpte.semilogx(self.period[pxy],
                                   rp[ii]['phasexy'][0][pxy],
                                   ls=':',
                                   marker='s',
                                   ms=self.ms,
                                   color='k',
                                   mfc='k')
                                   
                    axptm.semilogx(self.period[pyx],
                                   rp[ii]['phaseyx'][0][pyx],
                                   ls=':',
                                   marker='o',
                                   ms=self.ms,
                                   color='k',
                                   mfc='k')
                                 
                mpxy = [np.where(rp[ii]['phasexy'][2]!=0)[0]]
                mpyx = [np.where(rp[ii]['phaseyx'][2]!=0)[0]]
                axpte.semilogx(self.period[mpxy],
                               rp[ii]['phasexy'][2][mpxy],
                             ls='-',color=color_lst[jj])
                axptm.semilogx(self.period[mpyx],
                               rp[ii]['phaseyx'][2][mpyx],
                                 ls='-',
                                 color=color_lst[jj])
                           
            axrte.grid(True, alpha=.4)
            axrtm.grid(True, alpha=.4)
            
            
            axrtm.set_xticklabels(['' for ii in range(10)])
            axrte.set_xticklabels(['' for ii in range(10)])
            
            rmstestr = np.median(np.array(rmstestr)[1:])
            rmstmstr = np.median(np.array(rmstmstr)[1:])
            axrte.set_title('TE rms={0:.2f}'.format(rmstestr),
                            fontdict={'size':self.font_size+2,'weight':'bold'})
            axrtm.set_title('TM rms={0:.2f}'.format(rmstmstr),
                            fontdict={'size':self.font_size+2,'weight':'bold'})
            
            axpte.grid(True, alpha=.4)
            axpte.yaxis.set_major_locator(MultipleLocator(10))
            axpte.yaxis.set_minor_locator(MultipleLocator(1))
            
            axrte.set_ylabel('App. Res. ($\Omega \cdot m$)',
                           fontdict={'size':self.font_size+2,'weight':'bold'})
            axpte.set_ylabel('Phase (deg)',
                           fontdict={'size':self.font_size+2,'weight':'bold'})
            axpte.set_xlabel('Period (s)',
                             fontdict={'size':self.font_size+2,
                                       'weight':'bold'})
        
            axrte.yaxis.set_label_coords(-.08,.5)
            axpte.yaxis.set_label_coords(-.08,.5)
            
            axrtm.set_xticklabels(['' for ii in range(10)])
            axptm.grid(True,alpha=.4)
            axptm.yaxis.set_major_locator(MultipleLocator(10))
            axptm.yaxis.set_minor_locator(MultipleLocator(1))
            
            axrtm.set_ylabel('App. Res. ($\Omega \cdot m$)',
                           fontdict={'size':self.font_size+2,'weight':'bold'})
            axptm.set_ylabel('Phase (deg)',
                           fontdict={'size':self.font_size+2,'weight':'bold'})
            axptm.set_xlabel('Period (s)',
                             fontdict={'size':self.font_size+2,
                                       'weight':'bold'})
        
            axrtm.yaxis.set_label_coords(-.08,.5)
            axptm.yaxis.set_label_coords(-.08,.5)
            plt.suptitle(station,fontsize=self.font_size+2,fontweight='bold')
            plt.show()
            
        self.fig_lst.append({'station':station, 'fig':fig, 'axrte':axrte,
                             'axrtm':axrtm, 'axpte':axpte, 'axptm':axptm})
            
    def redraw_plot(self):
        """
        redraw plot if parameters change
        """
        
        plt.close(self.fig)
        self.plot()
        
    def save_figure(self, save_path, file_format='pdf', orientation='portrait', 
                  fig_dpi=None, close_fig='y'):
        """
        save_plot will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format
                            
                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path
                            
            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 
                              
            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait
                              
            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.
                          
            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open
                          
        :Example: ::
            
            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')
            
        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi
            
        if not os.path.exists(save_path):
            os.mkdir(save_path)
            
        for fdict in self.fig_lst:
            svfn = '{0}_responses.{1}'.format(fdict['station'], file_format)
            fdict['fig'].savefig(os.path.join(save_path, svfn), 
                                 dpi=self.fig_dpi, orientation=orientation)
        if close_fig == 'y':
            plt.close(fdict['fig'])
        
    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  
        
        Use this if you change an of the .fig or axes properties
        
        :Example: ::
            
            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()
        
        """

        self.fig.canvas.draw()
                          
    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """
        
        return ("Plots all the responses for all iterations in a given folder") 