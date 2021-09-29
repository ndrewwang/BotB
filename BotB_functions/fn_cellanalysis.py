def gravimetric_energy(cell):
    unit = cell.unit
    prop = cell.energy/cell.mass.total
    prop.ito(unit.W*unit.hr/unit.kg)
    return prop

def volumetric_energy(cell):
    unit = cell.unit
    prop = cell.energy/cell.volume
    prop.ito(unit.W*unit.hr/unit.L)
    return prop

#print
def print_cellresults(cell):
    print('\033[1m' + str(cell.name) + '\033[0m' + '   (' + str(cell.format) + ')')
    poscc = str(cell.cellstack.positive.currentcollector.name)
    pos = str(cell.cellstack.positive.composite.active.name)
    sep = str(cell.cellstack.separator.name)
    elyte = str(cell.cellstack.electrolyte.name)
    neg = str(cell.cellstack.negative.composite.active.name)
    negcc = str(cell.cellstack.negative.currentcollector.name)
    
    print(poscc + ' | ' + pos + ' | ' + sep + ' | ' + elyte + ' | ' + neg + ' | ' + negcc )
    print('============================================================')
    print('grav. energy dens.:' + str(gravimetric_energy(cell)))
    print('vol energy dens.:' + str(volumetric_energy(cell)))
    print('------------------------------------------------------------')
    print('cell energy: '   + str(cell.energy))
    print('cell capacity: '   + str(cell.capacity))
    print('cell mass: ' + str(cell.mass.total))
    print('np ratio: ' + str(cell.NPratio))
    
    
#PLOT ELECTRODE STACK THICKNESSES 
def plot_thickbreakdown(cell):
    import numpy as np 
    import matplotlib.pyplot as plt
    from dotmap import DotMap
    
    figsz = (11,1)
    #     figsz = (20,70)
    dpi = 100

    def get_nominal(value):
        try:
            nominal = value.n
        except AttributeError:
            nominal = value.magnitude
        return nominal

    unit = cell.unit
    if 'cellstack' in list(cell.keys()):
        cellstack = cell.cellstack
    else:
        cellstack = cell
        
    cellstack.separator.thick.ito(unit.um)
    sep = get_nominal(cellstack.separator.thick)

    cellstack.positive.composite.thick.ito(unit.um)
    pos = get_nominal(cellstack.positive.composite.thick)

    cellstack.positive.currentcollector.thick.ito(unit.um)
    poscc = get_nominal(cellstack.positive.currentcollector.thick)

    cellstack.negative.composite.thick.ito(unit.um)
    neg = get_nominal(cellstack.negative.composite.thick)

    cellstack.negative.currentcollector.thick.ito(unit.um)
    negcc = get_nominal(cellstack.negative.currentcollector.thick)

    thicks = np.array([sep, pos, poscc, pos, sep, neg, negcc, neg])

    cumthick = thicks.cumsum()

    category_colors = plt.get_cmap('RdYlGn')(
            np.linspace(0.15, 0.85, np.size(thicks)))

    labels = ['']
    height = 1

    category_names = ['separator','positive electrode','positive c.c.',
                      'positive electrode','separator','negative electrode',
                      'negative c.c.', 'negative electrode'] 


    fig, ax = plt.subplots(figsize=figsz, dpi=dpi) 


    i = 0
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='w',
            hatch='xxxx')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 1
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='#001a33',
            hatch='..',
            edgecolor='k')
    ax.barh(labels,thicc, left=starts, height=height, 
            color='none',edgecolor='k')

    i = 2
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='#7B7B7B')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 3
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            color='#001a33',
            hatch='..',
            edgecolor='k')
    ax.barh(labels,thicc, left=starts, height=height, 
            color='none',edgecolor='k')

    i = 4
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            color='w',
            hatch='xxxx')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 5
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='k',
            hatch='--',
            edgecolor='#7B7B7B')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 6
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    if cell.cellstack.negative.currentcollector.name=='Cu':
        cc_color = '#D08840'
    else: 
        cc_color = '#7B7B7B'

    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color=cc_color)
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 7
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            color='k',
            hatch='--',
            edgecolor='#7B7B7B')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    ax.set_ylim((-0.5,0.5))
    ax.set_xlim((0,cumthick[-1]))
    ax.set_xlabel('μm')
 
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.8),
              fancybox=True, ncol=5, frameon=False)  
    
def plot_massbreakdown(cell):
    import numpy as np 
    import matplotlib.pyplot as plt

    fontsize = 10
    figsz = (11,1)
    #     figsz = (20,70)
    dpi = 100

    def get_nominal(value):
        try:
            nominal = value.n
        except AttributeError:
            nominal = value.magnitude
        return nominal

    unit = cell.unit
    cellmass = cell.mass.total
    canmass = cell.mass.case
    posccmass = cell.mass.positivecc
    posmass = cell.mass.positive
    sepmass = cell.mass.separator
    elytemass = cell.mass.electrolyte
    negmass = cell.mass.negative
    negccmass = cell.mass.negativecc

    labels = ['']
    height = 1

    category_names = ['casing', 'positive c.c.','positive electrode','separator','electrolyte','negative electrode','negative c.c'] 
    thicks = np.array([100*get_nominal(canmass/cellmass),
                       100*get_nominal(posccmass/cellmass), 
                       100*get_nominal(posmass/cellmass), 
                       100*get_nominal(sepmass/cellmass), 
                       100*get_nominal(elytemass/cellmass), 
                       100*get_nominal(negmass/cellmass),
                       100*get_nominal(negccmass/cellmass)]) #array for mass
    cumthick = thicks.cumsum() #cumulative array for mass
    # width = 1
    # fig = figure()
    fig, ax = plt.subplots(figsize=figsz, dpi=dpi) 


    i = 0
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='#bfbfbf',
            hatch='')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 1
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='#7B7B7B')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 2
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='#001a33',
            hatch='..',
            edgecolor='k')
    ax.barh(labels,thicc, left=starts, height=height, 
            color='none',edgecolor='k')

    i = 3
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='w',
            hatch='xxxx')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 4
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='#b3e7ff',
            hatch='')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')


    i = 5
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='k',
            hatch='--',
            edgecolor='#7B7B7B')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 6
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    if cell.cellstack.negative.currentcollector.name=='Cu':
        cc_color = '#D08840'
    else: 
        cc_color = '#7B7B7B'



    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color=cc_color)
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')


    ax.set_ylim((-0.5,0.5))
    ax.set_xlim((0,cumthick[-1]))
    ax.set_xlabel('mass fraction (%)')



    # # ax.set_title('matplotlib.axes.Axes.barh Example') 
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1),
    #           fancybox=True, ncol=4)  
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.8),
              fancybox=True, ncol=4, frameon=False)  
    ax.set_xticks(np.linspace(0,100,11))
    plt.show() 
    
    
    
    
def plot_3D_cylindrical(cell):
    from matplotlib.pyplot import figure
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    fontsize = 10
    figsz = (5,5)
#     figsz = (20,70)
    dpi = 150
    
    def get_nominal(value):
        try:
            nominal = value.n
        except AttributeError:
            nominal = value.magnitude

        return nominal

    unit = cell.unit
    
    r_cell = cell.diameter/2
    r_cell.ito(unit.cm)
    r_cell = get_nominal(r_cell)
    
    h_cell = cell.height
    h_cell.ito(unit.cm)
    h_cell = get_nominal(h_cell)
    
    r_jroll_outer = (cell.diameter-2*cell.canthick)/2
    r_jroll_outer.ito(unit.cm)
    r_jroll_outer = get_nominal(r_jroll_outer)
    
    r_jroll_inner = cell.mandreldiam/2 #cm
    r_jroll_inner.ito(unit.cm)
    r_jroll_inner = get_nominal(r_jroll_inner)
    
    h_jroll = cell.height-cell.headspace-2*cell.canthick
    h_jroll.ito(unit.cm)
    h_jroll = get_nominal(h_jroll)
 
    layerthick = cell.stackthick
    layerthick.ito(unit.cm)
    layerthick = get_nominal(layerthick)

    def data_for_cylinder(center_x,center_y,radius,height_z):
        gridspacing = 100
        z = np.linspace(0, height_z, gridspacing)
        theta = np.linspace(0, 2*np.pi, gridspacing)
        theta_grid, z_grid=np.meshgrid(theta, z)
        x_grid = radius*np.cos(theta_grid) + center_x
        y_grid = radius*np.sin(theta_grid) + center_y
        return x_grid,y_grid,z_grid

    def data_for_spiral(outer_radii,inner_radii,layerthick,h_jroll):
        turns = (outer_radii-inner_radii)/layerthick
        th = 2*turns*np.pi
        Th = np.linspace(0,th,10000)
        x = (inner_radii + layerthick*Th/(2*np.pi))*np.cos(Th)
        y = (inner_radii + layerthick*Th/(2*np.pi))*np.sin(Th)
        z1 = np.linspace(0,0,np.size(x))
        z2 = np.linspace(h_jroll,h_jroll,np.size(x))
        return x,y,z1,z2

    #get data
    Xo,Yo,Zo = data_for_cylinder(0,0,r_cell,h_cell)
    Xi,Yi,Zi = data_for_cylinder(0,0,r_jroll_outer,h_jroll)
    x,y,z1,z2 = data_for_spiral(r_jroll_outer,r_jroll_inner,layerthick,h_jroll)

    fig = figure(figsize=figsz, dpi=dpi)
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(Xo, Yo, Zo, alpha=0.3,shade=True,facecolors=cm.bone(Xi),linewidth=0)
    ax.plot_surface(Xi, Yi, Zi, alpha=0.3,shade=True,facecolors=cm.bone(Xo),linewidth=0)

    ax.plot(x,y,z1,'k-',linewidth=0.3,alpha=1)
    ax.plot(x,y,z2,'k-',linewidth=0.3,alpha=1)
    
    d_lim = max(h_cell,2*r_cell)
    ax.set_xlim((-d_lim/2,d_lim/2))
    ax.set_ylim((-d_lim/2,d_lim/2))
    ax.set_zlim((0,d_lim))
    
#     ax.view_init(elev=90, azim=0)
#     ax.view_init(elev=0, azim=0)
    ax.view_init(elev=20, azim=0)
    
    ax.set_axis_off()
    
    
    #ANNOTATING ==================
    anno_width = 0.5
    anno_alpha = 0.6
    anno_fontsize = 6

    l1 = 0.2*r_cell #Longest line distance from cell edge
    lcap = abs(0.5*l1)

    #Cell height
    ax.plot([0,0],[l1+r_cell,l1+r_cell],[0,h_cell],'k-',linewidth=anno_width,alpha=anno_alpha)
    ax.plot([0,0],[-lcap+l1+r_cell,lcap+l1+r_cell],[h_cell,h_cell],'k-',linewidth=anno_width,alpha=anno_alpha)
    ax.plot([0,0],[-lcap+l1+r_cell,lcap+l1+r_cell],[0,0],'k-',linewidth=anno_width,alpha=anno_alpha)
    (cell.height.ito(unit.mm))
    txt = cell.height
    txt = (f'{txt:~}')
    ax.text(0,2*l1+r_cell,h_cell/2,
            txt,color='k',
           fontsize=anno_fontsize,
           zdir = ('z'),
           ha='center', va='center')
    
    #Cell Width
    ax.plot([l1+r_cell,l1+r_cell],[-r_cell,r_cell],[0,0],'k-',linewidth=anno_width,alpha=anno_alpha)
    ax.plot([-lcap+l1+r_cell,lcap+l1+r_cell],[r_cell,r_cell],[0,0],'k-',linewidth=anno_width,alpha=anno_alpha)
    ax.plot([-lcap+l1+r_cell,lcap+l1+r_cell],[-r_cell,-r_cell],[0,0],'k-',linewidth=anno_width,alpha=anno_alpha)
    (cell.diameter.ito(unit.mm))
    txt = cell.diameter
    txt = (f'{txt:~}')
    ax.text(3*l1+r_cell,0,0,
            txt,color='k',
           fontsize=anno_fontsize,
           zdir = ('y'),
           ha='center', va='center')
    
    
#     plt.savefig('pic.png',
#                 bbox_inches='tight',
#                 pad_inches=0,
#                 format='png',
#                 dpi=300)
    plt.show()

# 2D SPIRAL PLOT

def plot_2D_cylindrical(cell):
    from matplotlib.pyplot import figure
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    fontsize = 10
    figsz = (5,5)
    dpi = 100
    
    def get_nominal(value):
        try:
            nominal = value.n
        except AttributeError:
            nominal = value.magnitude
        return nominal

    unit = cell.unit
    
    #Get layer thicknesses
    r_cell = cell.diameter/2
    r_cell.ito(unit.mm)
    r_cell = get_nominal(r_cell)
    
    h_cell = cell.height
    h_cell.ito(unit.mm)
    h_cell = get_nominal(h_cell)
    
    r_jroll_outer = (cell.diameter-2*cell.canthick)/2
    r_jroll_outer.ito(unit.mm)
    r_jroll_outer = get_nominal(r_jroll_outer)
    
    r_jroll_inner = cell.mandreldiam/2 #cm
    r_jroll_inner.ito(unit.mm)
    r_jroll_inner = get_nominal(r_jroll_inner)
    
    h_jroll = cell.height-cell.headspace-2*cell.canthick
    h_jroll.ito(unit.mm)
    h_jroll = get_nominal(h_jroll)
 
    layerthick = cell.stackthick
    layerthick.ito(unit.mm)
    layerthick = get_nominal(layerthick)
    
    canthick = cell.canthick
    canthick.ito(unit.mm)
    canthick = get_nominal(canthick)

    def data_for_spiral(outer_radii,inner_radii,layerthick,h_jroll):
        turns = (outer_radii-inner_radii)/layerthick
        th = 2*turns*np.pi
        Th = np.linspace(0,th,10000)
        x = (inner_radii + layerthick*Th/(2*np.pi))*np.cos(Th)
        y = (inner_radii + layerthick*Th/(2*np.pi))*np.sin(Th)
        z1 = np.linspace(0,0,np.size(x))
        z2 = np.linspace(h_jroll,h_jroll,np.size(x))
        return x,y,z1,z2

    #get data
    x,y,z1,z2 = data_for_spiral(r_jroll_outer,r_jroll_inner,layerthick,h_jroll)
    
    theta = np.linspace(0, 2*np.pi, 100)

    fig = figure(figsize=figsz, dpi=dpi)
    ax = fig.add_subplot()
    ax.plot(x,y,'r-',linewidth=1,alpha=0.5)
    
    r = r_cell
    x1 = r*np.cos(theta)
    x2 = r*np.sin(theta)
    ax.plot(x1,x2,'k-',linewidth=1,alpha=1)
    
    r = r_cell-canthick
    x1 = r*np.cos(theta)
    x2 = r*np.sin(theta)
    ax.plot(x1,x2,'k-',linewidth=1,alpha=1)

    ax.get_yaxis().set_visible(False)
    ax.set_xlabel('mm')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
#     plt.savefig('chod.png',
#                 bbox_inches='tight',
#                 pad_inches=0,
#                 format='png',
#                 dpi=300)
    plt.show()
    
    
    
    
#PLOT POUCH CELL
def plot_3D_pouch(cell):
    from matplotlib.pyplot import figure
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    fontsize = 10
    figsz = (5,5)
    dpi = 150

    def get_nominal(value):
        try:
            nominal = value.n
        except AttributeError:
            nominal = value.magnitude
        return nominal

    unit = cell.unit

    #get data
    cell.height.ito(unit.cm)
    h_cell = get_nominal(cell.height)
    cell.width.ito(unit.cm)
    w_cell = get_nominal(cell.width)
    cell.depth.ito(unit.cm)
    t_cell = get_nominal(cell.depth)
    cell.pouchclearance.ito(unit.cm)
    celledge = get_nominal(cell.pouchclearance)/2
    
    fig = figure(figsize=figsz, dpi=dpi)
    ax = fig.add_subplot(111, projection='3d')

    facealpha = 0.2
    fcol = np.ones((10,10))*0.01

    # Big face 1
    (x, y) = np.meshgrid(np.linspace(0, h_cell, 10), np.linspace(0, w_cell, 10))
    z = np.ones(np.shape(x))*0
    ax.plot_surface(x, y, z, alpha=facealpha/2,shade=True,linewidth=0,facecolors=cm.bone(fcol))
    xx = [0,0,h_cell,h_cell,0]
    yy = [0,w_cell,w_cell,0,0]
    zz = [0,0,0,0,0]
    ax.plot(xx,yy,zz,'k-',linewidth=0.1,alpha=1)

    # Big face 2
    (x, y) = np.meshgrid(np.linspace(0, h_cell, 10), np.linspace(0, w_cell, 10))
    z = np.ones(np.shape(x))*t_cell
    ax.plot_surface(x, y, z, alpha=facealpha/2,shade=True,linewidth=0,facecolors=cm.bone(fcol))
    xx = [0,0,h_cell,h_cell,0]
    yy = [0,w_cell,w_cell,0,0]
    zz = [t_cell,t_cell,t_cell,t_cell,t_cell]
    ax.plot(xx,yy,zz,'k-',linewidth=0.1,alpha=1)

    # Edge Rings
    (z, y) = np.meshgrid(np.linspace(0, t_cell, 10), np.linspace(0, w_cell, 10))
    x = np.ones(np.shape(z))*0
    ax.plot_surface(x, y, z, alpha=facealpha,shade=True,linewidth=0,facecolors=cm.bone(fcol))
    xx = [0,0,0,0,0]
    yy = [0,w_cell,w_cell,0,0]
    zz = [0,0,t_cell,t_cell,0]
    ax.plot(xx,yy,zz,'k-',linewidth=0.1,alpha=1)
    (z, y) = np.meshgrid(np.linspace(0, t_cell, 10), np.linspace(0, w_cell, 10))
    x = np.ones(np.shape(z))*h_cell
    ax.plot_surface(x, y, z, alpha=facealpha,shade=True,linewidth=0,facecolors=cm.bone(fcol))
    xx = [h_cell,h_cell,h_cell,h_cell,h_cell]
    yy = [0,w_cell,w_cell,0,0]
    zz = [0,0,t_cell,t_cell,0]
    ax.plot(xx,yy,zz,'k-',linewidth=0.1,alpha=1)
    (z, x) = np.meshgrid(np.linspace(0, t_cell, 10), np.linspace(0, h_cell, 10))
    y = np.ones(np.shape(z))*0
    ax.plot_surface(x, y, z, alpha=facealpha,shade=True,linewidth=0,facecolors=cm.bone(fcol))
    (z, x) = np.meshgrid(np.linspace(0, t_cell, 10), np.linspace(0, h_cell, 10))
    y = np.ones(np.shape(z))*w_cell
    ax.plot_surface(x, y, z, alpha=facealpha,shade=True,linewidth=0,facecolors=cm.bone(fcol))

    #FLAT CRIMP EDGES
    (x, y) = np.meshgrid(np.linspace(0-celledge, h_cell+celledge, 10), np.linspace(0-celledge, w_cell+celledge, 10))
    z = np.ones(np.shape(x))*0.5*t_cell
    ax.plot_surface(x, y, z, alpha=facealpha,shade=True,linewidth=0,facecolors=cm.bone(fcol))
    xx = [-celledge,-celledge,h_cell+celledge,h_cell+celledge,-celledge]
    yy = [-celledge,w_cell+celledge,w_cell+celledge,-celledge,-celledge]
    zz = [0.5*t_cell,0.5*t_cell,0.5*t_cell,0.5*t_cell,0.5*t_cell]
    ax.plot(xx,yy,zz,'k-',linewidth=0.3,alpha=1)

    #LAYERS but you cant see them
    for i in range(cell.nlayers):
        z = (i/cell.nlayers)*cell.depth
        z = get_nominal(z)
        xx = [0,0,h_cell,h_cell,0]
        yy = [0,w_cell,w_cell,0,0]
        zz = [z,z,z,z,z]
        ax.plot(xx,yy,zz,'k-',linewidth=0.3,alpha=0.5)


    #TAB
    if cell.tabloc == 'top':
        
        tabh = get_nominal(cell.tabh)
        tabw = get_nominal(cell.tabw)
        (x, y) = np.meshgrid(np.linspace(h_cell, h_cell+tabh, 10), np.linspace((w_cell/3)-(tabw/2), (w_cell/3)+(tabw/2), 10))
        z = np.ones(np.shape(x))*0.5*t_cell
        ax.plot_surface(x, y, z, alpha=facealpha,shade=True,linewidth=0,facecolors=cm.bone(fcol))
        xx = [h_cell,h_cell,h_cell+tabh,h_cell+tabh,h_cell]
        yy = [(w_cell/3)-(tabw/2),(w_cell/3)+(tabw/2),(w_cell/3)+(tabw/2),(w_cell/3)-(tabw/2),(w_cell/3)-(tabw/2)]
        zz = [0.5*t_cell,0.5*t_cell,0.5*t_cell,0.5*t_cell,0.5*t_cell]
        ax.plot(xx,yy,zz,'k-',linewidth=0.5,alpha=1)

        (x, y) = np.meshgrid(np.linspace(h_cell, h_cell+tabh, 10), np.linspace((2*w_cell/3)-(tabw/2), (2*w_cell/3)+(tabw/2), 10))
        z = np.ones(np.shape(x))*0.5*t_cell
        ax.plot_surface(x, y, z, alpha=facealpha,shade=True,linewidth=0,facecolors=cm.bone(fcol))
        xx = [h_cell,h_cell,h_cell+tabh,h_cell+tabh,h_cell]
        yy = [(2*w_cell/3)-(tabw/2),(2*w_cell/3)+(tabw/2),(2*w_cell/3)+(tabw/2),(2*w_cell/3)-(tabw/2),(2*w_cell/3)-(tabw/2)]
        zz = [0.5*t_cell,0.5*t_cell,0.5*t_cell,0.5*t_cell,0.5*t_cell]
        ax.plot(xx,yy,zz,'k-',linewidth=0.5,alpha=1)
        
    elif cell.tabloc == 'sides':
        
        tabh = get_nominal(cell.tabh)
        tabw = get_nominal(cell.tabw)
        (x, y) = np.meshgrid(np.linspace(h_cell/2-tabw/2, h_cell/2+tabw/2, 10), np.linspace(w_cell, w_cell+tabh, 10))
        z = np.ones(np.shape(x))*0.5*t_cell
        ax.plot_surface(x, y, z, alpha=facealpha,shade=True,linewidth=0,facecolors=cm.bone(fcol))
        xx = [h_cell/2-tabw/2,h_cell/2+tabw/2,h_cell/2+tabw/2,h_cell/2-tabw/2,h_cell/2-tabw/2]
        yy = [w_cell,w_cell, w_cell+tabh, w_cell+tabh, w_cell]
        zz = [0.5*t_cell,0.5*t_cell,0.5*t_cell,0.5*t_cell,0.5*t_cell]
        ax.plot(xx,yy,zz,'k-',linewidth=0.5,alpha=1)

        (x, y) = np.meshgrid(np.linspace(h_cell/2-tabw/2, h_cell/2+tabw/2, 10), np.linspace(0, -tabh, 10))
        z = np.ones(np.shape(x))*0.5*t_cell
        ax.plot_surface(x, y, z, alpha=facealpha,shade=True,linewidth=0,facecolors=cm.bone(fcol))
        xx = [h_cell/2-tabw/2,h_cell/2+tabw/2,h_cell/2+tabw/2,h_cell/2-tabw/2,h_cell/2-tabw/2]
        yy = [0,0, -tabh, -tabh, 0]
        zz = [0.5*t_cell,0.5*t_cell,0.5*t_cell,0.5*t_cell,0.5*t_cell]
        ax.plot(xx,yy,zz,'k-',linewidth=0.5,alpha=1)

    
    d_lim = max(h_cell,w_cell)
    pad = 0.1
    ax.set_xlim((-pad*d_lim,(1+pad)*d_lim))
    ax.set_ylim((-pad*d_lim,(1+pad)*d_lim))
    ax.set_zlim(((-0.5-pad)*d_lim,(pad+0.5)*d_lim))

#     ax.view_init(elev=220, azim=0)
    ax.view_init(elev=235, azim=10)
    ax.set_axis_off()


    #     plt.savefig('test.png',
    #                 bbox_inches='tight',
    #                 pad_inches=0,
    #                 format='png',
    #                 dpi=300)
    
    #ANNOTATING ==================
    anno_width = 0.5
    anno_alpha = 0.6
    anno_fontsize = 6
    
    l1 = -0.15*h_cell #Longest line distance from cell edge
    if w_cell > h_cell:
        l1 = -0.15*w_cell
    lcap = abs(0.1*l1)
    
    #Cell Width
    ax.plot([l1,l1],[0,w_cell],[0.5*t_cell,0.5*t_cell],'k-',linewidth=anno_width,alpha=anno_alpha)
    ax.plot([l1-lcap,l1+lcap],[0,0],[0.5*t_cell,0.5*t_cell],'k-',linewidth=anno_width,alpha=anno_alpha)
    ax.plot([l1-lcap,l1+lcap],[w_cell,w_cell],[0.5*t_cell,0.5*t_cell],'k-',linewidth=anno_width,alpha=anno_alpha)
    (cell.width.ito(unit.mm))
    txt = cell.width
    txt = (f'{txt:~}')
    ax.text(1.25*l1,0.5*w_cell,0.5*t_cell,
            txt,color='k',
           fontsize=anno_fontsize,
           zdir = ('y'),
           ha='center', va='center')
    
    #Cell Height
    
    if cell.tabloc == 'sides':
        l1 = -0.15*w_cell-tabh
        lcap = abs(0.1*l1)
        
    ax.plot([0,h_cell],[l1,l1],[0.5*t_cell,0.5*t_cell],'k-',linewidth=anno_width,alpha=anno_alpha)
    ax.plot([h_cell,h_cell],[l1-lcap,l1+lcap],[0.5*t_cell,0.5*t_cell],'k-',linewidth=anno_width,alpha=anno_alpha)
    ax.plot([0,0],[l1-lcap,l1+lcap],[0.5*t_cell,0.5*t_cell],'k-',linewidth=anno_width,alpha=anno_alpha)
    (cell.height.ito(unit.mm))
    txt = cell.height
    txt = (f'{txt:~}')
    ax.text(0.5*h_cell,1.25*l1,0.5*t_cell,
            txt,color='k',
           fontsize=anno_fontsize,
           zdir = ('x'),
           ha='center', va='center')
    #     plt.show()
