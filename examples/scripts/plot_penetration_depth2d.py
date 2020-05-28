#! /usr/bin/env python
"""
Description:
    Plot 2D penetration depth for a folder of EDI files.

CreationDate:   23/03/2018
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     23/03/2018   FZ

    brenainn.moushall@ga.gov.au 03-04-2020 15:39:33 AEDT:
        - Add selection of periods by specifying periods in seconds
        - Clean up script
"""
from mtpy.imaging import penetration_depth2d as pen2d

edidir = '/path/to/edi/files'

# selected_periods: the periods in seconds to plot depth for across each
# station.
selected_periods = [10., 100., 500., 600.]
# ptol: tolerance to use when finding nearest period to each selected
# period. If abs(selected period - nearest period) is greater than
# selected period * ptol, then the period is discarded and will appear
# as a gap in the plot.
ptol = 0.20
# zcomponent: component to plot. Valid parameters are 'det, 'zxy' and
# 'zyx'
zcomponent = 'det'  # 'zxy', 'zyx' also options

pen2d.plot2Dprofile(edi_dir=edidir,
                    selected_periods=selected_periods,
                    ptol=ptol,
                    zcomponent=zcomponent,
                    save=True,
                    savepath='/tmp/Depth2D.png')

# selected_period_indices: indices of periods to plot.
# 'ptol' ins't required if using indices.
selected_period_indices = [0, 10, 20]

# p_index: needs to be set to True if using indices.
period_by_index = True

pen2d.plot2Dprofile(edi_dir=edidir,
                    selected_periods=selected_period_indices,
                    period_by_index=period_by_index,
                    zcomponent=zcomponent,
                    save=True,
                    savepath='/tmp/Depth2D_by_index.png')
