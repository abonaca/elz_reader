from orbit_fits import *
from astropy.io import fits

def elz():
    """"""
    
    # streams
    names = get_names()
    ts = Table.read('../../disrupted_gc/data/overall_summary.fits')
    ts['lperp'] = np.sqrt(ts['lx']**2 + ts['ly']**2)
    
    # H3 stars
    th = Table.read('../data/rcat_giants.fits')
    
    # gloublar clusters
    tgc = Table.read('../../disrupted_gc/data/gc_orbits.fits')
    tgc['lperp'] = np.sqrt(tgc['ly']**2 + tgc['lx']**2)
    
    gc_progenitor = ['NGC 3201', 'NGC 4590', 'NGC 5024', 'NGC 5139', 'NGC 5272', 'NGC 5824']
    ind_progenitor = np.array([True if x in gc_progenitor else False for x in tgc['name']])
    indices = [~ind_progenitor, ind_progenitor]
    
    ms = [6, 18]
    marker = ['o', '*']
    color = ['none', '#11d5c7']
    labels = ['Globular clusters', 'Tentative stream\nprogenitors']
    
    # plotting
    plt.close()
    fig, ax = plt.subplots(1,2,figsize=(16,7))
    
    # background H3 stars
    plt.sca(ax[0])
    plt.plot(th['Lz'], th['E_tot_pot1'], '.', color='0.4', alpha=0.5, mew=0, ms=3, label='')
    plt.plot(10, 10, '.', color='0.4', alpha=1, mew=0, ms=3, label='Giant stars [H3]')
    
    # background globular clusters
    plt.sca(ax[1])
    plt.scatter(tgc['lz'], tgc['etot'], marker='o', s=40, c=tgc['lperp'], edgecolors='none', linewidths=2, cmap='magma', vmin=0, vmax=6, label='')
    plt.plot(tgc['lz'][ind_progenitor], tgc['etot'][ind_progenitor], 'o', color=color[1], zorder=0, ms=8, label='')
    
    # legend entries
    f_color = 0.6
    plt.plot(10, 10, '*', color=mpl.cm.magma(f_color), ms=18, mec='k', label='Streams')
    plt.plot(10, 10, 'o', color=mpl.cm.magma(f_color), mec='none', ms=10, label='Globular clusters')
    plt.plot(10, 10, 'o', color=mpl.cm.magma(f_color), mec=color[1], ms=10, mew=2, label='Tentative stream\nprogenitors')
    
    offsets = get_offsets()
    
    for i, name in enumerate(names[:]):
        t = Table.read('../data/output/orbit_props_{:s}.fits'.format(name))
        
        lperp = np.nanmedian(np.sqrt(t['lx']**2 + t['ly']**2))
        color = lperp/6.
        
        # plot samples from the posterior
        plt.sca(ax[0])
        plt.plot(t['lz'], t['etot'], '.', alpha=0.5, ms=3, mew=0, color=mpl.cm.magma(color), rasterized=True, label='')
        
        # place labels
        lz = ts['lz'][i] + offsets[name][0]
        etot = ts['etot'][i] + offsets[name][1]
        
        for j in range(2):
            plt.sca(ax[j])
            #plt.text(lz, etot, name, fontsize='x-small')
            plt.text(lz, etot, '${:s}$'.format(get_properties(name)['label']), fontsize='x-small', alpha=1, va='center', ha=offsets[name][2])
    
    
    plt.sca(ax[1])
    plt.scatter(ts['lz'], ts['etot'], c=ts['lperp'], s=400, marker='*', edgecolors='k', cmap='magma', vmin=0, vmax=6, label='')
    
    for i in range(2):
        plt.sca(ax[i])
        plt.xlim(-4.5,4.5)
        plt.ylim(-0.265, -0.03)
        plt.ylim(-0.18, -0.03)
        
        plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
        plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
        
        plt.legend(loc=3, fontsize='small', handlelength=0.5)
    
        # custom colorbar
        sm = plt.cm.ScalarMappable(cmap=plt.cm.magma, norm=plt.Normalize(vmin=0, vmax=6))
        # fake up the array of the scalar mappable. Urgh…
        sm._A = []
        plt.colorbar(sm, label='$L_\perp$ [kpc$^2$ Myr$^{-1}$]', pad=0.02, aspect=30)

    plt.tight_layout()
    plt.savefig('../paper/elz_streams.pdf')

def get_offsets():
    """Returns a dictionary for placing stream labels in the ELz"""
    
    offsets = dict()
    offsets['leiptr'] = [0.0, -0.005, 'left']
    offsets['gjoll'] = [0.15, -0.002, 'left']
    offsets['gd1'] = [0.15, -0.002, 'left']
    offsets['phlegethon'] = [0.0, 0.005, 'center']
    offsets['ylgr'] = [0.15, -0.002, 'left']
    offsets['wambelong'] = [0.0, -0.005, 'left']
    offsets['fimbulthul'] = [0.15, -0.002, 'left']
    offsets['ophiuchus'] = [0.0, -0.005, 'center']
    offsets['elqui'] = [0.15, -0.002, 'left']
    offsets['svol'] = [0.0, -0.004, 'right']
    offsets['ravi'] = [-0.1, 0.002, 'right']
    offsets['sylgr'] = [0.15, -0.002, 'left']
    offsets['jhelum'] = [0.15, -0.002, 'left']
    offsets['indus'] = [0.15, -0.002, 'left']
    offsets['phoenix'] = [0.0, -0.004, 'right']
    offsets['slidr'] = [0.15, 0.002, 'left']
    offsets['atlas'] = [0.1, -0.003, 'left']
    offsets['aliqa_uma'] = [0.15, -0.003, 'left']
    offsets['turbio'] = [-0.15, 0.00, 'right']
    offsets['turranburra'] = [-0.0, -0.003, 'right']
    offsets['fjorm'] = [0.0, -0.004, 'right']
    offsets['triangulum'] = [0.2, -0.005, 'center']
    offsets['willka_yaku'] = [-0.2, 0.005, 'center']
    
    return offsets

def elz_progenitors():
    """Figure 2: streams in ELz with relevant Naidu+20 substructure contours"""
    
    # streams
    ts = Table.read('../../disrupted_gc/data/overall_summary.fits')
    names = get_names()
    
    # plot configs
    labels = ['gse', 'sgr', 'thamnos', 'sequoia', 'iitoi', 'arjuna', 'helmi',  'wukong']
    titles = ['Gaia Enceladus', 'Sagittarius', 'Thamnos', 'Sequoia', "I'itoi", 'Arjuna', 'Helmi', 'Wukong']
    colors = ['#ffa22f', '#223195', '#3dabdb', '#ba532e', '#9a4526', '#6c301b', '#ca5a87', '#752c84']
    
    lw = 2
    alpha = 0.8
    
    plt.close()
    fig, ax = plt.subplots(1,1,figsize=(6.5,6))
    
    # plot progenitor contours
    for i, label in enumerate(labels):
        kde = np.load('../data/kde_{:s}.npz'.format(label))
        color = colors[i]

        plt.contour(kde['X'][0], kde['Y'][0], kde['Z'][0], colors=color, levels=[0.2*np.max(kde['Z'][0])], linewidths=lw, alpha=alpha)
        plt.plot([10,20], [10,10], '-', color=color, lw=lw, alpha=alpha, label=titles[i])

    # plot streams
    plt.plot(ts['lz'], ts['etot'], '*', ms=16, mec='k', mew=1.2, alpha=0.7, color='none', label='')
    
    # label streams
    offsets = get_offsets()
    
    for i, name in enumerate(names[:]):
        lz = ts['lz'][i] + offsets[name][0]
        etot = ts['etot'][i] + offsets[name][1]
        
        plt.text(lz, etot, '${:s}$'.format(get_properties(name)['label']), fontsize='x-small', alpha=0.7, va='center', ha=offsets[name][2])
    
    plt.legend(handlelength=1.5, ncol=3, loc=9, fontsize='x-small', framealpha=0.8)
    
    plt.xlim(-4.5,4.5)
    plt.ylim(-0.17, -0.041)
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
    
    plt.tight_layout()
    #plt.savefig('../paper/stream_hosts.pdf')

def elz_origin():
    """Figure 2: streams in ELz with relevant Naidu+20 substructure contours"""
    
    # streams
    ts = Table.read('../../disrupted_gc/data/overall_summary.fits')
    names = get_names()
    
    to = Table.read('../data/stream_origin.fits')
    ind_dg = to['type']=='DG'
    ind_gc = ~ind_dg
    ind_type = [ind_gc, ind_dg]
    marker_type = ['*', 'p']
    ms_type = [16, 12]
    label_type = ['GC', 'DG']
    
    # plot configs
    labels = ['gse', 'sgr', 'thamnos', 'sequoia', 'iitoi', 'arjuna', 'helmi',  'wukong']
    titles = ['Gaia Enceladus', 'Sagittarius', 'Thamnos', 'Sequoia', "I'itoi", 'Arjuna', 'Helmi', 'Wukong']
    colors = ['#ffa22f', '#223195', '#3dabdb', '#ba532e', '#9a4526', '#6c301b', '#ca5a87', '#752c84']
    
    lw = 2
    alpha = 0.8
    
    labels_origin = ['insitu', 'exsitu']
    titles_origin = ['In-situ', 'Ex-situ']
    colors_origin = ['0.6', '0.2']
    ls_origin = ['--', '-']
    
    
    plt.close()
    fig, ax = plt.subplots(1,2,figsize=(13,6))
    
    plt.sca(ax[0])
    
    # plot origin contours
    for i, label in enumerate(labels_origin[:1]):
        kde = np.load('../data/kde_{:s}.npz'.format(label))
        color = colors_origin[i]

        plt.contour(kde['X'][0], kde['Y'][0], kde['Z'][0], linestyles=ls_origin[i], colors=color, levels=[0.2*np.max(kde['Z'][0])], linewidths=lw, alpha=alpha)
        plt.plot([10,20], [10,10], ls=ls_origin[i], color=color, lw=lw, alpha=alpha, label=titles_origin[i])

    
    plt.sca(ax[1])
    
    # plot progenitor contours
    for i, label in enumerate(labels):
        kde = np.load('../data/kde_{:s}.npz'.format(label))
        color = colors[i]

        plt.contour(kde['X'][0], kde['Y'][0], kde['Z'][0], colors=color, levels=[0.2*np.max(kde['Z'][0])], linewidths=lw, alpha=alpha)
        plt.plot([10,20], [10,10], '-', color=color, lw=lw, alpha=alpha, label=titles[i])

    for j in range(2):
        plt.sca(ax[j])
        
        # plot streams
        for e in range(2):
            if j==0:
                plt.plot(ts['lz'][ind_type[e]], ts['etot'][ind_type[e]], marker_type[e], ms=ms_type[e], mec='k', mew=1.2, alpha=0.7, color='none', label='{:s} stream'.format(label_type[e]))
            else:
                plt.plot(ts['lz'][ind_type[e]], ts['etot'][ind_type[e]], marker_type[e], ms=ms_type[e], mec='k', mew=1.2, alpha=0.7, color='none', label='')
        #plt.plot(ts['lz'], ts['etot'], '*', ms=16, mec='k', mew=1.2, alpha=0.7, color='none', label='')
        
        # label streams
        offsets = get_offsets()
        
        for i, name in enumerate(names[:]):
            lz = ts['lz'][i] + offsets[name][0]
            etot = ts['etot'][i] + offsets[name][1]
            
            plt.text(lz, etot, '${:s}$'.format(get_properties(name)['label']), fontsize='x-small', alpha=0.7, va='center', ha=offsets[name][2])
            
        #if j>0:
        plt.legend(handlelength=1.5, ncol=3, loc=9, fontsize='x-small', framealpha=0.8)
        
        plt.text(0.07,0.07, titles_origin[j], fontsize='medium', transform=plt.gca().transAxes)
        
        plt.xlim(-4.5,4.5)
        plt.ylim(-0.17, -0.041)
        plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
        plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
    
    plt.tight_layout()
    plt.savefig('../paper/stream_hosts.pdf')
    
def sky_orbits():
    """Connect globular clusters to streams"""
    
    t = Table.read('/home/ana/data/baumgardt_positions.fits')
    
    clusters = ['NGC 3201', 'NGC 4590', 'NGC 5824', 'NGC 5139', 'NGC 5272', 'NGC 5024']
    N = len(clusters)
    
    match = dict()
    match['NGC 3201'] = dict(streams=['gjoll'], direction=[-1], nstep=[35])
    match['NGC 4590'] = dict(streams=['fjorm'], direction=[1], nstep=[100])
    match['NGC 5024'] = dict(streams=['sylgr', 'ravi'], direction=[-1, 1], nstep=[300,500])
    match['NGC 5139'] = dict(streams=['fimbulthul'], direction=[-1], nstep=[70])
    match['NGC 5272'] = dict(streams=['svol'], direction=[1], nstep=[70])
    match['NGC 5824'] = dict(streams=['triangulum', 'turbio'], direction=[1,1], nstep=[600,1])
    
    dt = 0.5*u.Myr
    wangle = 180*u.deg
    ra_off = 120*u.deg
    l_off = 0*u.deg
    
    colors = [mpl.cm.plasma(0.95*x/N) for x in range(N)]
    
    plt.close()
    fig = plt.figure(figsize=(12,12))
    
    ax0 = fig.add_subplot(211, projection='mollweide')
    ax1 = fig.add_subplot(212, projection='mollweide')
    ax = [ax0, ax1]
    
    for i in range(N):
        #ind = t['Name']== clusters[i]
        ind = t['Name']==clusters[i]
        t_ = t[ind]
        
        c = coord.SkyCoord(ra=t_['RAJ2000'], dec=t_['DEJ2000'], distance=t_['Rsun'], pm_ra_cosdec=t_['pmRA_'], pm_dec=t_['pmDE'], radial_velocity=t_['RV'], frame='icrs')
        cgal = c.transform_to(coord.Galactic)
        w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)[0]
        
        color = colors[i]
        alpha_text = 0.8
        
        plt.sca(ax[0])
        plt.plot((c.ra+ra_off).wrap_at(wangle).rad, c.dec.rad, '+', color=color, mew=2, ms=10, label=t_['Name'][0])
        plt.text((c.ra+ra_off).wrap_at(wangle).rad, (c.dec+3*u.deg).rad, t_['Name'][0], fontsize='small', ha='center', alpha=alpha_text)
        
        plt.sca(ax[1])
        plt.plot((cgal.l+l_off).wrap_at(wangle).rad, cgal.b.rad, '+', color=color, mew=2, ms=10, label=t_['Name'][0])
        plt.text((cgal.l+l_off).wrap_at(wangle).rad, (cgal.b+3*u.deg).rad, t_['Name'][0], fontsize='small', ha='center', alpha=alpha_text)
        

        for j in range(len(match[clusters[i]]['direction'])):
            orbit = ham.integrate_orbit(w0, dt=dt*match[clusters[i]]['direction'][j], n_steps=match[clusters[i]]['nstep'][j])
            orbit_eq = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=gc_frame)
            orbit_gal = orbit.to_coord_frame(coord.Galactic, galactocentric_frame=gc_frame)
            
                
            plt.sca(ax[0])
            plt.plot((orbit_eq.ra+ra_off).wrap_at(wangle).rad, orbit_eq.dec.rad, '-', color=color, lw=1, label='')
            
            plt.sca(ax[1])
            dl = orbit_gal.l.wrap_at(wangle)[1:] - orbit_gal.l.wrap_at(wangle)[:-1]
            if np.any(np.abs(dl)>180*u.deg):
                pos_break = dl>180*u.deg
                ind_break = np.argmax(pos_break)
                ipad = 1
                plt.plot((orbit_gal.l+l_off).wrap_at(wangle).rad[:ind_break-ipad], orbit_gal.b.rad[:ind_break-ipad], '-', color=color, lw=1, label='')
                plt.plot((orbit_gal.l+l_off).wrap_at(wangle).rad[ind_break+ipad:], orbit_gal.b.rad[ind_break+ipad:], '-', color=color, lw=1, label='')
            else:
                plt.plot((orbit_gal.l+l_off).wrap_at(wangle).rad, orbit_gal.b.rad, '-', color=color, lw=1, label='')
            
            # add streams
            pkl = pickle.load(open('../data/streams/data_{:s}.pkl'.format(match[clusters[i]]['streams'][j]), 'rb'))
            cs = coord.SkyCoord(ra=pkl['dec'][0], dec=pkl['dec'][1], frame='icrs')
            cs_gal = cs.transform_to(coord.Galactic)
            
            plt.sca(ax[0])
            plt.plot((cs.ra+ra_off).wrap_at(wangle).rad, cs.dec.rad, 'o', color=color, ms=8, label=match[clusters[i]]['streams'][j])
            plt.text((cs.ra+ra_off).wrap_at(wangle).rad[-1], cs.dec.rad[-1], match[clusters[i]]['streams'][j], fontsize='small', alpha=alpha_text)
            
            plt.sca(ax[1])
            plt.plot((cs_gal.l+l_off).wrap_at(wangle).rad, cs_gal.b.rad, 'o', color=color, ms=8, label=match[clusters[i]]['streams'][j])
            plt.text((cs_gal.l+l_off).wrap_at(wangle).rad[-1], cs_gal.b.rad[-1], match[clusters[i]]['streams'][j], fontsize='small', alpha=alpha_text)
    
    
    plt.sca(ax[0])
    plt.grid(ls=':')
    plt.xlabel('R.A. + {:.0f} [deg]'.format(ra_off.value))
    plt.ylabel('Dec [deg]')
    
    plt.sca(ax[1])
    plt.grid(ls=':')
    plt.xlabel('Galactic longitude [deg]')
    plt.ylabel('Galactic latitude [deg]')
    
    
    plt.tight_layout(h_pad=2)
    plt.savefig('../paper/sky_orbits.pdf')

def table_summary():
    """create a table with summary of the stream origin"""
    
    t = dict()
    t['name'] = get_names()
    t['Name'] = [get_properties(name)['label'] for name in t['name']]
    N = len(t['name'])
    
    # host
    t['host'] = ['Sagittarius', 'Sagittarius', 'itself', 'Gaia Enceladus', 'Sagittarius', 'Sequoia / Arjuna / I\'itoi', 'Sequoia', '(Gaia Enceladus)', '(Gaia Enceladus)', 'Sequoia', 'Gaia Enceladus', 'Sequoia', 'Helmi/Wukong', 'Helmi/Wukong', 'Sagittarius', 'in situ', 'Helmi/Wukong', 'Cetus', 'Cetus', 'Sagittarius', 'Sequoia / Arjuna / I\'itoi', 'Cetus', 'Sequoia / Arjuna / I\'itoi']
    
    # progenitor
    t['progenitor'] = [np.nan, np.nan, 'itself', 'NGC 5139', 'NGC 4590', np.nan, 'NGC 3201', 'itself', 'itself', np.nan, np.nan, np.nan, np.nan, 'NGC 5024', 'NGC 5272', np.nan, 'NGC 5024', 'NGC 5824', 'NGC 5824', np.nan, np.nan, np.nan, np.nan]
    
    # progenitor type
    t['type'] = ['DG'  if name in ['elqui', 'indus', 'jhelum'] else 'GC' for name in t['name']]
    
    # metallicity
    t['feh'] = [-2.4, -2.4, -2.2, -1.5, -2.16, -2.3, -1.5, -2.1, -2.1, -1.6, -1.95, -2.7, -1.6, -1.1, -1.7, -1.6, -2.7, -1.9, np.nan, np.nan, -2.2, np.nan, -1.9]
    
    tout = Table(t)
    tout.pprint()
    tout.write('../data/stream_origin.fits', overwrite=True)
    
def latex_table():
    """Create a latex table"""
    
    t = Table.read('../data/stream_origin.fits')
    N = len(t)
    
    f = open('../paper/stream_origin.tex', 'w')
    for i in range(N):
        t_ = t[i]
        f.write('{:s} & {:s} & {:s} & {:s} & {:.1f}\\\\ \n'.format(t_['Name'], t_['host'], t_['progenitor'], t_['type'], t_['feh']))
        
    f.close()

