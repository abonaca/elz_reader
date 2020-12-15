from orbit_fits import *
from astropy.io import fits
from scipy.stats import gaussian_kde

def elz_rand():
    """"""
    
    # streams
    ts = Table.read('../../disrupted_gc/data/overall_summary.fits')
    
    names = get_names()
    
    # add globular clusters
    tgc = Table.read('../../disrupted_gc/data/gc_orbits.fits')

    gc_sequoia = np.array(['NGC 5466', 'NGC 7006', 'IC 4499'])
    ind_sequoia = np.array([True if x in gc_sequoia else False for x in tgc['name']])
    
    gc_ge = np.array(['NGC 288', 'NGC 362', 'NGC 1261', 'NGC 1851', 'NGC 1904', 'NGC 2298', 'NGC 2808', 'NGC 4147', 'NGC 4833', 'NGC 5286', 'NGC 5897', 'NGC 6205', 'NGC 6235', 'NGC 6284', 'NGC 6341', 'NGC 6779', 'NGC 6864', 'NGC 7089', 'NGC 7099', 'NGC 7492'])
    ind_ge = np.array([True if x in gc_ge else False for x in tgc['name']])
    
    gc_helmi = np.array(['NGC 4590', 'NGC 5024', 'NGC 5053', 'NGC 5272', 'NGC 6981'])
    ind_helmi = np.array([True if x in gc_helmi else False for x in tgc['name']])
    
    gc_sgr = np.array(['NGC 2419', 'NGC 5824', 'NGC 6715', 'Pal 12', 'Ter 7', 'Ter 8', 'Arp 2', 'Whiting 1'])
    ind_sgr = np.array([True if x in gc_sgr else False for x in tgc['name']])

    gc_kraken = np.array(['NGC 5946', 'NGC 5986', 'NGC 6093', 'NGC 6121', 'NGC 6144', 'NGC 6254', 'NGC 6273', 'NGC 6287', 'NGC 6541', 'NGC 6544', 'NGC 6681', 'NGC 6712', 'NGC 6809'])
    ind_kraken = np.array([True if x in gc_kraken else False for x in tgc['name']])
    
    ind_smooth = ~ind_sequoia & ~ind_ge & ~ind_helmi & ~ind_sgr & ~ind_kraken
    
    labels = ['Sequoia', 'Gaia - Enceladus', 'Helmi', 'Sagittarius', 'Kraken', 'Unassociated']
    indices = [ind_sequoia, ind_ge, ind_helmi, ind_sgr, ind_kraken, ind_smooth]
    markers = ['^', 'H', 'p', 'D', 's', 'o']
    sizes = [60, 80, 80, 40, 60, 20]
    logm = np.array([7.9, 8.43, 7.96, 8.44, 8.28, -1])
    masses = 10**logm
    
    labels_masses = ['{:s} ({:.0f}$\cdot$10$^{:.0f}$ M$_\odot$)'.format(labels[i], masses[i]*10**-int(logm[i]), int(logm[i])) for i in range(5)]
    #masses = [7.9e7, 2.7e8, 9.1e7, 2.7e8, 1.9e8, 0]
    labels_masses += ['Unassociated']
    uscale = 1.

    np.random.seed(182)

    # plotting
    
    plt.close()
    plt.figure(figsize=(9,9))
    
    for i in range(6):
        plt.scatter(tgc['lz'][indices[i]], tgc['etot'][indices[i]]*uscale, c=0*tgc['lz'][indices[i]], vmin=-2, vmax=2, zorder=0, cmap='PuOr', s=1.5*sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label='')
    
        # legend entries
        plt.scatter(0, 0, facecolors='w', s=sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label=labels[i])
    
    ## label globular clusters
    #if labeled:
        #for i in range(len(tgc)):
            #if (tgc['lz'][i]>2) & (tgc['lz'][i]<4):
                #plt.text(tgc['lz'][i], tgc['etot'][i]*uscale, '{:s}'.format(tgc['name'][i]), fontsize='x-small')

    plt.legend(loc=4, frameon=False, fontsize='x-small', handlelength=0.6, title='Globular clusters', title_fontsize='x-small')

    
    for name in names[:]:
        t = Table.read('../data/output/orbit_props_{:s}.fits'.format(name))
        plt.plot(t['lz'], t['etot'], '.', alpha=0.2, ms=3, mew=0, label=name, color=mpl.cm.magma(np.random.rand()*0.9))
        
        lz = np.nanmedian(t['lz']) + 0.01
        etot = np.nanmedian(t['etot']) + 0.001
        plt.text(lz, etot, name, fontsize='x-small')
    
    #plt.plot(ts['Lz'][:imax], ts['Etot'][:imax], 'ko')
    plt.xlim(-4,4)
    plt.ylim(-0.265, -0.01)
    
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
    
    plt.tight_layout()
    plt.savefig('../plots/elz_streams_dr3_rand.png')

def h3_giants():
    """"""
    
    t = Table(fits.getdata('../data/rcat.fits'))
    ind = (t['FLAG']==0) & (t['logg']<3.5)
    t = t[ind]
    
    t['Lx'] = (t['Lx']*u.km/u.s*u.kpc).to(u.kpc**2/u.Myr)
    t['Ly'] = (t['Ly']*u.km/u.s*u.kpc).to(u.kpc**2/u.Myr)
    t['Lz'] = (t['Lz']*u.km/u.s*u.kpc).to(u.kpc**2/u.Myr)
    t['Lperp'] = np.sqrt(t['Lx']**2 + t['Ly']**2)
    t['E_tot_pot1'] = (t['E_tot_pot1']*u.km**2*u.s**-2).to(u.kpc**2*u.Myr**-2)
    
    t.write('../data/rcat_giants.fits', overwrite=True)

def h3_substructure():
    """"""
    t = Table.read('/home/ana/data/substructure_rcat_full.fits')
    
    t['Lx'] = (t['Lx']*u.km/u.s*u.kpc*1e3).to(u.kpc**2/u.Myr)
    t['Ly'] = (t['Ly']*u.km/u.s*u.kpc*1e3).to(u.kpc**2/u.Myr)
    t['Lz'] = (t['Lz']*u.km/u.s*u.kpc*1e3).to(u.kpc**2/u.Myr)
    t['Lperp'] = np.sqrt(t['Lx']**2 + t['Ly']**2)
    t['E_tot_pot1'] = (t['E_tot_pot1']*u.km**2*u.s**-2*1e5).to(u.kpc**2*u.Myr**-2)
    print(t['Lperp'])
    
    t.write('../data/rcat_substructure.fits', overwrite=True)


def elz_h3():
    """"""
    
    # streams
    names = get_names()
    
    th = Table.read('../data/rcat_giants.fits')
    
    # plotting
    
    plt.close()
    plt.figure(figsize=(9.5,9))
    
    plt.plot(th['Lz'], th['E_tot_pot1'], '.', color='0.4', alpha=0.5, mew=0, ms=3)
    
    for name in names[:]:
        t = Table.read('../data/output/orbit_props_{:s}.fits'.format(name))
        
        lperp = np.nanmedian(np.sqrt(t['lx']**2 + t['ly']**2))
        color = lperp/6.
        plt.plot(t['lz'], t['etot'], '.', alpha=0.5, ms=3, mew=0, label=name, color=mpl.cm.magma(color))
        
        lz = np.nanmedian(t['lz']) + 0.03
        etot = np.nanmedian(t['etot']) - 0.003
        plt.text(lz, etot, name, fontsize='x-small')
    
    #plt.plot(ts['Lz'][:imax], ts['Etot'][:imax], 'ko')
    plt.xlim(-4.5,4.5)
    plt.ylim(-0.265, -0.03)
    plt.ylim(-0.18, -0.03)
    
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
    

    # custom colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.magma, norm=plt.Normalize(vmin=0, vmax=6))
    # fake up the array of the scalar mappable. Urgh…
    sm._A = []
    plt.colorbar(sm, label='$L_\perp$ [kpc$^2$ Myr$^{-1}$]', pad=0.02, aspect=30)

    plt.tight_layout()
    
    plt.savefig('../plots/elz_streams_h3_dr3.png')

def elz_h3_ly():
    """"""
    
    # streams
    names = get_names()
    
    th = Table.read('../data/rcat_giants.fits')
    
    # plotting
    cmap = mpl.cm.viridis
    
    plt.close()
    fig, ax = plt.subplots(1,1,figsize=(9.5,9))
    
    plt.plot(th['Lz'], th['E_tot_pot1'], '.', color='0.4', alpha=0.5, mew=0, ms=3)
    
    for name in names[:]:
        t = Table.read('../data/output/orbit_props_{:s}.fits'.format(name))
        
        lperp = np.nanmedian(np.sqrt(t['lx']**2 + t['ly']**2))
        ly = np.nanmedian(t['ly']) + 4
        color = ly/8.
        if color<0: color = 0
        if color>1: color = 1
        plt.plot(t['lz'], t['etot'], '.', alpha=0.5, ms=3, mew=0, label=name, color=cmap(color))
        
        lz = np.nanmedian(t['lz']) + 0.03
        etot = np.nanmedian(t['etot']) - 0.003
        plt.text(lz, etot, name, fontsize='x-small')
    
    #plt.plot(ts['Lz'][:imax], ts['Etot'][:imax], 'ko')
    plt.xlim(-4.5,4.5)
    plt.ylim(-0.265, -0.03)
    plt.ylim(-0.18, -0.03)
    
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
    

    # custom colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-4, vmax=4))
    # fake up the array of the scalar mappable. Urgh…
    sm._A = []
    plt.colorbar(sm, label='$L_y$ [kpc$^2$ Myr$^{-1}$]', pad=0.02, aspect=30)

    plt.tight_layout()
    
    plt.savefig('../plots/elz_streams_ly.png')

def rapo_ecc():
    """"""
    
    # streams
    names = get_names()
    
    # plotting
    
    plt.close()
    fig, ax = plt.subplots(1,2,figsize=(10,5.5), sharey=True)
    
    
    for name in names[:]:
        t = Table.read('../data/output/orbit_props_{:s}.fits'.format(name))
        lperp = np.nanmedian(np.sqrt(t['lx']**2 + t['ly']**2))
        color = lperp/6.
        
        plt.sca(ax[0])
        plt.plot(t['rperi'], t['ecc'], '.', alpha=0.5, ms=3, mew=0, label=name, color=mpl.cm.magma(color))
        
        xname = np.nanmedian(t['rperi']) + 0.5
        yname = np.nanmedian(t['ecc']) + 0.01
        plt.text(xname, yname, name, fontsize='x-small')
        
        plt.sca(ax[1])
        plt.plot(t['rapo'], t['ecc'], '.', alpha=0.5, ms=3, mew=0, label=name, color=mpl.cm.magma(color))
        
        xname = np.nanmedian(t['rapo']) + 1
        yname = np.nanmedian(t['ecc']) - 0.01
        plt.text(xname, yname, name, fontsize='x-small')
    
    plt.sca(ax[0])
    plt.xlim(0,30)
    plt.ylim(0,1)
    plt.xlabel('$r_{peri}$ [kpc]')
    plt.ylabel('Eccentricity')
    
    plt.sca(ax[1])
    plt.xlim(0,70)
    plt.xlabel('$r_{apo}$ [kpc]')
    
    # custom colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.magma, norm=plt.Normalize(vmin=0, vmax=6))
    # fake up the array of the scalar mappable. Urgh…
    sm._A = []
    #plt.colorbar(sm, label='$L_\perp$ [kpc$^2$ Myr$^{-1}$]', aspect=30)

    plt.tight_layout(w_pad=0.1)
    
    pos_l = ax[0].get_position()
    pos_r = ax[1].get_position()
    cax = plt.axes([pos_l.x0, pos_l.y1+0.015, pos_r.x1-pos_l.x0, 0.025])
    plt.colorbar(sm, cax=cax, label='$L_\perp$ [kpc$^2$ Myr$^{-1}$]', orientation='horizontal', ticklocation='top')

    
    plt.savefig('../plots/rapo_ecc.png')



def elz():
    """"""
    
    # streams
    ts = Table.read('../../disrupted_gc/data/overall_summary.fits')
    
    names = get_names()
    
    # add globular clusters
    tgc = Table.read('../../disrupted_gc/data/gc_orbits.fits')

    gc_sequoia = np.array(['NGC 5466', 'NGC 7006', 'IC 4499'])
    ind_sequoia = np.array([True if x in gc_sequoia else False for x in tgc['name']])
    
    gc_ge = np.array(['NGC 288', 'NGC 362', 'NGC 1261', 'NGC 1851', 'NGC 1904', 'NGC 2298', 'NGC 2808', 'NGC 4147', 'NGC 4833', 'NGC 5286', 'NGC 5897', 'NGC 6205', 'NGC 6235', 'NGC 6284', 'NGC 6341', 'NGC 6779', 'NGC 6864', 'NGC 7089', 'NGC 7099', 'NGC 7492'])
    ind_ge = np.array([True if x in gc_ge else False for x in tgc['name']])
    
    gc_helmi = np.array(['NGC 4590', 'NGC 5024', 'NGC 5053', 'NGC 5272', 'NGC 6981'])
    ind_helmi = np.array([True if x in gc_helmi else False for x in tgc['name']])
    
    gc_sgr = np.array(['NGC 2419', 'NGC 5824', 'NGC 6715', 'Pal 12', 'Ter 7', 'Ter 8', 'Arp 2', 'Whiting 1'])
    ind_sgr = np.array([True if x in gc_sgr else False for x in tgc['name']])

    gc_kraken = np.array(['NGC 5946', 'NGC 5986', 'NGC 6093', 'NGC 6121', 'NGC 6144', 'NGC 6254', 'NGC 6273', 'NGC 6287', 'NGC 6541', 'NGC 6544', 'NGC 6681', 'NGC 6712', 'NGC 6809'])
    ind_kraken = np.array([True if x in gc_kraken else False for x in tgc['name']])
    
    ind_smooth = ~ind_sequoia & ~ind_ge & ~ind_helmi & ~ind_sgr & ~ind_kraken
    
    labels = ['Sequoia', 'Gaia - Enceladus', 'Helmi', 'Sagittarius', 'Kraken', 'Unassociated']
    indices = [ind_sequoia, ind_ge, ind_helmi, ind_sgr, ind_kraken, ind_smooth]
    markers = ['^', 'H', 'p', 'D', 's', 'o']
    sizes = [60, 80, 80, 40, 60, 20]
    logm = np.array([7.9, 8.43, 7.96, 8.44, 8.28, -1])
    masses = 10**logm
    
    labels_masses = ['{:s} ({:.0f}$\cdot$10$^{:.0f}$ M$_\odot$)'.format(labels[i], masses[i]*10**-int(logm[i]), int(logm[i])) for i in range(5)]
    #masses = [7.9e7, 2.7e8, 9.1e7, 2.7e8, 1.9e8, 0]
    labels_masses += ['Unassociated']
    uscale = 1.

    # plotting
    
    plt.close()
    plt.figure(figsize=(9.5,9))
    
    for i in range(6):
        plt.scatter(tgc['lz'][indices[i]], tgc['etot'][indices[i]]*uscale, c=0*tgc['lz'][indices[i]], vmin=-2, vmax=2, zorder=0, cmap='PuOr', s=1.5*sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label='')
    
        # legend entries
        plt.scatter(0, 0, facecolors='w', s=sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label=labels[i])
    
    ## label globular clusters
    #if labeled:
        #for i in range(len(tgc)):
            #if (tgc['lz'][i]>2) & (tgc['lz'][i]<4):
                #plt.text(tgc['lz'][i], tgc['etot'][i]*uscale, '{:s}'.format(tgc['name'][i]), fontsize='x-small')

    plt.legend(loc=4, frameon=False, fontsize='x-small', handlelength=0.6, title='Globular clusters', title_fontsize='x-small')

    
    for name in names[:]:
        t = Table.read('../data/output/orbit_props_{:s}.fits'.format(name))
        
        lperp = np.nanmedian(np.sqrt(t['lx']**2 + t['ly']**2))
        color = lperp/6.
        plt.plot(t['lz'], t['etot'], '.', alpha=0.2, ms=3, mew=0, label=name, color=mpl.cm.magma(color))
        
        lz = np.nanmedian(t['lz']) + 0.03
        etot = np.nanmedian(t['etot']) - 0.003
        plt.text(lz, etot, name, fontsize='x-small')
    
    #plt.plot(ts['Lz'][:imax], ts['Etot'][:imax], 'ko')
    plt.xlim(-4.5,4.5)
    plt.ylim(-0.265, -0.03)
    plt.ylim(-0.18, -0.03)
    
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
    

    # custom colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.magma, norm=plt.Normalize(vmin=0, vmax=6))
    # fake up the array of the scalar mappable. Urgh…
    sm._A = []
    plt.colorbar(sm, label='$L_\perp$ [kpc$^2$ Myr$^{-1}$]', pad=0.02, aspect=30)

    plt.tight_layout()
    
    plt.savefig('../plots/elz_streams_dr3.png')

def lylz():
    """"""
    
    # streams
    ts = Table.read('../../disrupted_gc/data/overall_summary.fits')
    names = get_names()
    
    # add globular clusters
    tgc = Table.read('../../disrupted_gc/data/gc_orbits.fits')

    gc_sequoia = np.array(['NGC 5466', 'NGC 7006', 'IC 4499'])
    ind_sequoia = np.array([True if x in gc_sequoia else False for x in tgc['name']])
    
    gc_ge = np.array(['NGC 288', 'NGC 362', 'NGC 1261', 'NGC 1851', 'NGC 1904', 'NGC 2298', 'NGC 2808', 'NGC 4147', 'NGC 4833', 'NGC 5286', 'NGC 5897', 'NGC 6205', 'NGC 6235', 'NGC 6284', 'NGC 6341', 'NGC 6779', 'NGC 6864', 'NGC 7089', 'NGC 7099', 'NGC 7492'])
    ind_ge = np.array([True if x in gc_ge else False for x in tgc['name']])
    
    gc_helmi = np.array(['NGC 4590', 'NGC 5024', 'NGC 5053', 'NGC 5272', 'NGC 6981'])
    ind_helmi = np.array([True if x in gc_helmi else False for x in tgc['name']])
    
    gc_sgr = np.array(['NGC 2419', 'NGC 5824', 'NGC 6715', 'Pal 12', 'Ter 7', 'Ter 8', 'Arp 2', 'Whiting 1'])
    ind_sgr = np.array([True if x in gc_sgr else False for x in tgc['name']])

    gc_kraken = np.array(['NGC 5946', 'NGC 5986', 'NGC 6093', 'NGC 6121', 'NGC 6144', 'NGC 6254', 'NGC 6273', 'NGC 6287', 'NGC 6541', 'NGC 6544', 'NGC 6681', 'NGC 6712', 'NGC 6809'])
    ind_kraken = np.array([True if x in gc_kraken else False for x in tgc['name']])
    
    ind_smooth = ~ind_sequoia & ~ind_ge & ~ind_helmi & ~ind_sgr & ~ind_kraken
    
    labels = ['Sequoia', 'Gaia - Enceladus', 'Helmi', 'Sagittarius', 'Kraken', 'Unassociated']
    indices = [ind_sequoia, ind_ge, ind_helmi, ind_sgr, ind_kraken, ind_smooth]
    markers = ['^', 'H', 'p', 'D', 's', 'o']
    sizes = [60, 80, 80, 40, 60, 20]
    logm = np.array([7.9, 8.43, 7.96, 8.44, 8.28, -1])
    masses = 10**logm
    
    labels_masses = ['{:s} ({:.0f}$\cdot$10$^{:.0f}$ M$_\odot$)'.format(labels[i], masses[i]*10**-int(logm[i]), int(logm[i])) for i in range(5)]
    #masses = [7.9e7, 2.7e8, 9.1e7, 2.7e8, 1.9e8, 0]
    labels_masses += ['Unassociated']
    uscale = 1.

    np.random.seed(182)

    # plotting
    
    plt.close()
    plt.figure(figsize=(9.5,9))
    
    for i in range(6):
        plt.scatter(tgc['lz'][indices[i]], tgc['ly'][indices[i]]*uscale, c=0*tgc['lz'][indices[i]], vmin=-2, vmax=2, zorder=0, cmap='PuOr', s=1.5*sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label='')
    
        # legend entries
        plt.scatter(-100, -100, facecolors='w', s=sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label=labels[i])
    
    ## label globular clusters
    #if labeled:
        #for i in range(len(tgc)):
            #if (tgc['lz'][i]>2) & (tgc['lz'][i]<4):
                #plt.text(tgc['lz'][i], tgc['etot'][i]*uscale, '{:s}'.format(tgc['name'][i]), fontsize='x-small')

    plt.legend(loc=1, frameon=False, fontsize='x-small', handlelength=0.6, title='Globular clusters', title_fontsize='x-small')

    
    for name in names[:]:
        t = Table.read('../data/output/orbit_props_{:s}.fits'.format(name))
        
        lperp = np.nanmedian(np.sqrt(t['lx']**2 + t['ly']**2))
        color = lperp/6.
        plt.plot(t['lz'], t['ly'], '.', alpha=0.2, ms=3, mew=0, label=name, color=mpl.cm.magma(color))
        
        lz = np.nanmedian(t['lz']) + 0.01
        ly = np.nanmedian(t['ly']) + 0.01
        plt.text(lz, ly, name, fontsize='x-small')
    
    lz_ = np.linspace(-10,10,100)
    ly_ = -2.5 - 0.3*lz_
    plt.plot(lz_, ly_, 'k:', lw=1)
    
    plt.xlim(-7,7)
    plt.ylim(-7,7)
    
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$L_y$ [kpc$^2$ Myr$^{-1}$]')
    
    # custom colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.magma, norm=plt.Normalize(vmin=0, vmax=6))
    # fake up the array of the scalar mappable. Urgh…
    sm._A = []
    plt.colorbar(sm, label='$L_\perp$ [kpc$^2$ Myr$^{-1}$]', pad=0.02, aspect=30)
    
    plt.tight_layout()
    plt.savefig('../plots/lylz_streams_dr3.png')

def lylx():
    """"""
    
    # streams
    ts = Table.read('../../disrupted_gc/data/overall_summary.fits')
    names = get_names()
    
    # add globular clusters
    tgc = Table.read('../../disrupted_gc/data/gc_orbits.fits')

    gc_sequoia = np.array(['NGC 5466', 'NGC 7006', 'IC 4499'])
    ind_sequoia = np.array([True if x in gc_sequoia else False for x in tgc['name']])
    
    gc_ge = np.array(['NGC 288', 'NGC 362', 'NGC 1261', 'NGC 1851', 'NGC 1904', 'NGC 2298', 'NGC 2808', 'NGC 4147', 'NGC 4833', 'NGC 5286', 'NGC 5897', 'NGC 6205', 'NGC 6235', 'NGC 6284', 'NGC 6341', 'NGC 6779', 'NGC 6864', 'NGC 7089', 'NGC 7099', 'NGC 7492'])
    ind_ge = np.array([True if x in gc_ge else False for x in tgc['name']])
    
    gc_helmi = np.array(['NGC 4590', 'NGC 5024', 'NGC 5053', 'NGC 5272', 'NGC 6981'])
    ind_helmi = np.array([True if x in gc_helmi else False for x in tgc['name']])
    
    gc_sgr = np.array(['NGC 2419', 'NGC 5824', 'NGC 6715', 'Pal 12', 'Ter 7', 'Ter 8', 'Arp 2', 'Whiting 1'])
    ind_sgr = np.array([True if x in gc_sgr else False for x in tgc['name']])

    gc_kraken = np.array(['NGC 5946', 'NGC 5986', 'NGC 6093', 'NGC 6121', 'NGC 6144', 'NGC 6254', 'NGC 6273', 'NGC 6287', 'NGC 6541', 'NGC 6544', 'NGC 6681', 'NGC 6712', 'NGC 6809'])
    ind_kraken = np.array([True if x in gc_kraken else False for x in tgc['name']])
    
    ind_smooth = ~ind_sequoia & ~ind_ge & ~ind_helmi & ~ind_sgr & ~ind_kraken
    
    labels = ['Sequoia', 'Gaia - Enceladus', 'Helmi', 'Sagittarius', 'Kraken', 'Unassociated']
    indices = [ind_sequoia, ind_ge, ind_helmi, ind_sgr, ind_kraken, ind_smooth]
    markers = ['^', 'H', 'p', 'D', 's', 'o']
    sizes = [60, 80, 80, 40, 60, 20]
    logm = np.array([7.9, 8.43, 7.96, 8.44, 8.28, -1])
    masses = 10**logm
    
    labels_masses = ['{:s} ({:.0f}$\cdot$10$^{:.0f}$ M$_\odot$)'.format(labels[i], masses[i]*10**-int(logm[i]), int(logm[i])) for i in range(5)]
    #masses = [7.9e7, 2.7e8, 9.1e7, 2.7e8, 1.9e8, 0]
    labels_masses += ['Unassociated']
    uscale = 1.

    np.random.seed(182)

    # plotting
    
    plt.close()
    plt.figure(figsize=(9.5,9))
    
    for i in range(6):
        plt.scatter(tgc['lx'][indices[i]], tgc['ly'][indices[i]]*uscale, c=0*tgc['lz'][indices[i]], vmin=-2, vmax=2, zorder=0, cmap='PuOr', s=1.5*sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label='')
    
        # legend entries
        plt.scatter(-100, -100, facecolors='w', s=sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label=labels[i])
    
    ## label globular clusters
    #if labeled:
        #for i in range(len(tgc)):
            #if (tgc['lz'][i]>2) & (tgc['lz'][i]<4):
                #plt.text(tgc['lz'][i], tgc['etot'][i]*uscale, '{:s}'.format(tgc['name'][i]), fontsize='x-small')

    plt.legend(loc=1, frameon=False, fontsize='x-small', handlelength=0.6, title='Globular clusters', title_fontsize='x-small')

    
    for name in names[:]:
        t = Table.read('../data/output/orbit_props_{:s}.fits'.format(name))
        
        lperp = np.nanmedian(np.sqrt(t['lx']**2 + t['ly']**2))
        color = lperp/6.
        plt.plot(t['lx'], t['ly'], '.', alpha=0.2, ms=3, mew=0, label=name, color=mpl.cm.magma(color))
        
        lx = np.nanmedian(t['lx']) + 0.01
        ly = np.nanmedian(t['ly']) + 0.01
        plt.text(lx, ly, name, fontsize='x-small')
    
    #lz_ = np.linspace(-10,10,100)
    #ly_ = -2.5 - 0.3*lz_
    #plt.plot(lz_, ly_, 'k:', lw=1)
    
    plt.xlim(-5,9)
    plt.ylim(-7,7)
    
    plt.xlabel('$L_x$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$L_y$ [kpc$^2$ Myr$^{-1}$]')
    
    # custom colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.magma, norm=plt.Normalize(vmin=0, vmax=6))
    # fake up the array of the scalar mappable. Urgh…
    sm._A = []
    plt.colorbar(sm, label='$L_\perp$ [kpc$^2$ Myr$^{-1}$]', pad=0.02, aspect=30)

    plt.tight_layout()
    plt.savefig('../plots/lylx_streams_dr3.png')


def orbits_sky(T=1000):
    """"""
    t = Table.read('../data/fit_summary_fiducial.fits')
    N = len(t)
    
    dt = 0.1*u.Myr
    T_fwd = T*u.Myr
    nstep_fwd = int((T_fwd/dt).decompose())
    
    T_rr = T*u.Myr
    nstep_rr = int((T_rr/dt).decompose())
    
    wangle = 180*u.deg
    
    plt.close()
    fig = plt.figure(figsize=(12,5.2))
    ax = fig.add_subplot(111, projection='mollweide')
    
    for i in range(N):
        t_ = t[i]
        c = coord.SkyCoord(ra=t_['ra']*u.deg, dec=t_['dec'][0]*u.deg, distance=t_['dist'][0]*u.kpc, pm_ra_cosdec=t_['pmra'][0]*u.mas/u.yr, pm_dec=t_['pmdec'][0]*u.mas/u.yr, radial_velocity=t_['vr'][0]*u.km/u.s, frame='icrs')
        w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)

        orbit_fwd = ham.integrate_orbit(w0, dt=dt, n_steps=nstep_fwd)
        orbiteq_fwd = orbit_fwd.to_coord_frame(coord.ICRS, galactocentric_frame=gc_frame)
        orbitgal_fwd = orbit_fwd.to_coord_frame(coord.Galactic, galactocentric_frame=gc_frame)
        
        orbit_rr = ham.integrate_orbit(w0, dt=-dt, n_steps=nstep_rr)
        orbiteq_rr = orbit_rr.to_coord_frame(coord.ICRS, galactocentric_frame=gc_frame)
        orbitgal_rr = orbit_rr.to_coord_frame(coord.Galactic, galactocentric_frame=gc_frame)
    
        ts = Table.read('../data/output/orbit_props_{:s}.fits'.format(t_['name']))
        lperp = np.nanmedian(np.sqrt(ts['lx']**2 + ts['ly']**2))
        color = lperp/6.
    
        # plot
        #plt.plot(orbiteq_fwd.ra.wrap_at(wangle).rad, orbiteq_fwd.dec.rad, 'o', ms=1, mew=0, color=mpl.cm.magma(color))
        #plt.plot(orbiteq_rr.ra.wrap_at(wangle).rad, orbiteq_rr.dec.rad, 'o', ms=1, mew=0, color=mpl.cm.magma(color))
        
        plt.plot(orbitgal_fwd.l.wrap_at(wangle).rad, orbitgal_fwd.b.rad, 'o', ms=1, mew=0, color=mpl.cm.magma(color))
        plt.plot(orbitgal_rr.l.wrap_at(wangle).rad, orbitgal_rr.b.rad, 'o', ms=1, mew=0, color=mpl.cm.magma(color))
    
    
    plt.xlabel('l [deg]')
    plt.ylabel('b [deg]')
    
    sm = plt.cm.ScalarMappable(cmap=plt.cm.magma, norm=plt.Normalize(vmin=0, vmax=6))
    sm._A = []
    cb = fig.colorbar(sm, ax=ax, pad=0.04, aspect=20)
    cb.set_label('$L_\perp$ [kpc$^2$ Myr$^{-1}$]')
    
    
    plt.tight_layout()
    plt.savefig('../plots/orbits_sky_{:04.0f}.png'.format(T_fwd.value))


def phase_space():
    """"""
    # streams
    ts = Table.read('../../disrupted_gc/data/overall_summary.fits')
    names = get_names()
    
    # add globular clusters
    tgc = Table.read('../../disrupted_gc/data/gc_orbits.fits')

    labels = ['Streams (GC)', 'Streams (dwarf)', 'Streams (?)']
    ind = [ts['progenitor']=='gc', ts['progenitor']=='dwarf', ts['progenitor']=='n/a']
    marker = ['*', 'o', 's']
    ms = [12, 8, 7]
    color_gc = '0.4'
    
    plt.close()
    fig, ax = plt.subplots(1,3,figsize=(15,5))
    
    plt.sca(ax[0])
    plt.plot(tgc['lz'], tgc['etot'], '+', color=color_gc, ms=5, mew=0.5, label='Globular clusters')
    
    plt.xlim(-4.5,4.5)
    plt.ylim(-0.18, -0.03)
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
    
    plt.sca(ax[1])
    plt.plot(tgc['lz'], tgc['ly'], '+', color=color_gc, ms=5, mew=0.5, label='Globular clusters')
    
    plt.xlim(-4.5, 4.5)
    plt.ylim(-7,7)
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$L_y$ [kpc$^2$ Myr$^{-1}$]')
    
    plt.sca(ax[2])
    plt.plot(tgc['lx'], tgc['ly'], '+', color=color_gc, ms=5, mew=0.5, label='Globular clusters')
    
    plt.xlim(-5.5, 5.5)
    plt.ylim(-7,7)
    plt.xlabel('$L_x$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$L_y$ [kpc$^2$ Myr$^{-1}$]')
    
    
    for i in range(3):
        plt.sca(ax[0])
        plt.plot(ts['lz'][ind[i]], ts['etot'][ind[i]], marker=marker[i], ms=ms[i], mec='k', color='none', label=labels[i])
    
        plt.sca(ax[1])
        plt.plot(ts['lz'][ind[i]], ts['ly'][ind[i]], marker=marker[i], ms=ms[i], mec='k', color='none', label=labels[i])
        
        plt.sca(ax[2])
        plt.plot(ts['lx'][ind[i]], ts['ly'][ind[i]], marker=marker[i], ms=ms[i], mec='k', color='none', label=labels[i])
    
    plt.sca(ax[1])
    plt.legend(fontsize='x-small', loc=4, handlelength=0.6)
    
    plt.tight_layout()
    plt.savefig('../plots/phase_space.png')

def phase_space_h3():
    """"""
    # streams
    ts = Table.read('../../disrupted_gc/data/overall_summary.fits')
    names = get_names()
    
    # h3 giants
    th = Table.read('../data/rcat_substructure.fits')
    
    plt.close()
    fig, ax = plt.subplots(1,3,figsize=(18,6))
    
    plt.sca(ax[0])
    
    plt.xlim(-4.5,4.5)
    plt.ylim(-0.18, -0.03)
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
    
    plt.sca(ax[1])
    
    plt.xlim(-4.5, 4.5)
    plt.ylim(-7,7)
    plt.ylim(0, 7)
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$L_y$ [kpc$^2$ Myr$^{-1}$]')
    
    plt.sca(ax[2])
    
    plt.xlim(-5.5, 5.5)
    plt.ylim(-7,7)
    plt.xlabel('$L_x$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$L_y$ [kpc$^2$ Myr$^{-1}$]')
    
    prog_ids = ['GSE', 'HS', 'Iitoi', 'Sequo', 'Arjun', 'Sgr', 'Thamn', 'Wuk']
    ms_ = 1
    
    for pid in prog_ids:
        ind = th['Substructure_ID']==pid
        t_ = th[ind]
        
        plt.sca(ax[0])
        plt.plot(t_['Lz'], t_['E_tot_pot1'], 'o', ms=ms_, label=pid)
        
        plt.sca(ax[1])
        plt.plot(t_['Lz'], np.sqrt(t_['Ly']**2 + t_['Lx']**2), 'o', ms=ms_, label=pid)
        
        plt.sca(ax[2])
        plt.plot(t_['Lx'], t_['Ly'], 'o', ms=ms_, label=pid)
    
    labels = ['GC', 'dwarf', '?']
    ind = [ts['progenitor']=='gc', ts['progenitor']=='dwarf', ts['progenitor']=='n/a']
    marker = ['*', 'o', 's']
    ms = [12, 8, 7]
    
    for i in range(3):
        plt.sca(ax[0])
        plt.plot(ts['lz'][ind[i]], ts['etot'][ind[i]], marker=marker[i], ms=ms[i], mec='k', color='none', label=labels[i])
    
        plt.sca(ax[1])
        plt.plot(ts['lz'][ind[i]], np.sqrt(ts['ly'][ind[i]]**2 + ts['lx'][ind[i]]**2), marker=marker[i], ms=ms[i], mec='k', color='none', label=labels[i])
        
        plt.sca(ax[2])
        plt.plot(ts['lx'][ind[i]], ts['ly'][ind[i]], marker=marker[i], ms=ms[i], mec='k', color='none', label=labels[i])
    
    for name in names[:]:
        t = Table.read('../data/output/orbit_props_{:s}.fits'.format(name))
        
        #lperp = np.nanmedian(np.sqrt(t['lx']**2 + t['ly']**2))
        #color = lperp/6.
        #plt.plot(t['lx'], t['ly'], '.', alpha=0.2, ms=3, mew=0, label=name, color=mpl.cm.magma(color))
        
        lz = np.nanmedian(t['lz']) + 0.03
        etot = np.nanmedian(t['etot']) - 0.003
        
        plt.sca(ax[0])
        plt.text(lz, etot, name, fontsize='x-small')
        
        lx = np.nanmedian(t['lx']) + 0.01
        ly = np.nanmedian(np.sqrt(t['ly']**2 + t['lx']**2)) + 0.01
        lz = np.nanmedian(t['lz']) + 0.01
        
        plt.sca(ax[1])
        plt.text(lz, ly, name, fontsize='x-small')
    
        #plt.sca(ax[2])
        #plt.text(lx, ly, name, fontsize='x-small')
    
    #plt.sca(ax[1])
    #plt.legend(fontsize='x-small', loc=2, ncol=3, handlelength=0.6)
    
    plt.tight_layout()
    #plt.savefig('../plots/phase_space_progenitors.png')


def cetus():
    """"""
    giants = Table.read('../data/rcat_giants.fits')
    cetus_mask = (giants['Ly']>0.7e3*u.kpc*u.km/u.s) & (giants['Lx']>2e3*u.kpc*u.km/u.s) & (giants['FeH']<-1.5) & (giants['Lz']<-1.4e3*u.kpc*u.km/u.s) & (giants['Lz']>-3.4e3*u.kpc*u.km/u.s) & (giants['E_tot_pot1']>-1e5*u.km**2*u.s**-2) & (giants['E_tot_pot1']<-0.6e5*u.km**2*u.s**-2) & (giants['R_gal']>20) & (giants['eccen_pot1']<0.5) & (giants['RA']<150)
    
    print(giants['Lz'])
    print(np.sum(cetus_mask))
    
    plt.close()
    plt.figure()
    
    plt.plot(giants['Lz'], giants['E_tot_pot1'], 'k.', ms=1)
    plt.plot(giants['Lz'][cetus_mask], giants['E_tot_pot1'][cetus_mask], 'ro')
    
    plt.xlim(-5,5)
    plt.ylim(-0.18, -0.03)
    
    plt.tight_layout()
    

def origin_kde():
    """Save smooth distributions for in-situ and ex-situ components"""
    
    # h3 giants
    th = Table.read('../data/rcat_substructure.fits')
    print(np.array(np.unique(th['Substructure_ID'])))
    
    #prog_ids = ['GSE', 'HS', 'Iitoi', 'Sequo', 'Arjun', 'Sgr', 'Thamn', 'Wuk']
    #labels = ['gse', 'helmi', 'iitoi', 'sequoia', 'arjuna', 'sgr', 'thamnos', 'wukong']
    
    ind_in = (th['Substructure_ID']=='HiAlp') | (th['Substructure_ID']=='Aleph') #| (th['Substructure_ID']=='MWTD')
    ind_ex = ~(th['Substructure_ID']=='uncla') & ~((th['Substructure_ID']=='Aleph') | (th['Substructure_ID']=='HiAlp') | (th['Substructure_ID']=='MWTD'))
    ind = [ind_in, ind_ex]
    labels = ['insitu', 'exsitu']
    
    cx = ['Lz', 'Lz', 'Lz', 'Lx']
    cy = ['E_tot_pot1', 'Lperp', 'Ly', 'Ly']
    
    # set up boundaries of the 2D output array
    lz_min, lz_max = -4.5, 4.5
    ly_min, ly_max = -7, 7
    lx_min, lx_max = -5.5, 5.5
    lperp_min, lperp_max = 0, 7
    e_min, e_max = -0.18, -0.03
    
    xmin = np.array([lz_min, lz_min, lz_min, lx_min])
    xmax = np.array([lz_max, lz_max, lz_max, lx_max])
    ymin = np.array([e_min, lperp_min, ly_min, ly_min])
    ymax = np.array([e_max, lperp_max, ly_max, ly_max])
    
    # set up output arrays
    X1, Y1 = np.mgrid[lz_min:lz_max:200j, e_min:e_max:200j]
    X2, Y2 = np.mgrid[lz_min:lz_max:200j, lperp_min:lperp_max:200j]
    X3, Y3 = np.mgrid[lz_min:lz_max:200j, ly_min:ly_max:200j]
    X4, Y4 = np.mgrid[lx_min:lx_max:200j, ly_min:ly_max:200j]
    X = np.array([X1, X2, X3, X4])
    Y = np.array([Y1, Y2, Y3, Y4])
    
    positions1 = np.vstack([X1.ravel(), Y1.ravel()])
    positions2 = np.vstack([X2.ravel(), Y2.ravel()])
    positions3 = np.vstack([X3.ravel(), Y3.ravel()])
    positions4 = np.vstack([X4.ravel(), Y4.ravel()])
    positions = np.array([positions1, positions2, positions3, positions4])
    
    # figure sizing
    nrow = 4
    ncol = len(labels)
    da = 2.5
    
    plt.close()
    fig, ax = plt.subplots(nrow, ncol, figsize=(ncol*da, nrow*da), sharex='row', sharey='row')
    
    for i in range(ncol):
        #ind = th['Substructure_ID']==prog_ids[i]
        t_ = th[ind[i]]
        Z = np.zeros_like(X)
        
        for j in range(nrow):
            values = np.vstack([t_[cx[j]], t_[cy[j]]])
            kernel = gaussian_kde(values)
            Z[j] = np.reshape(kernel(positions[j]).T, X[j].shape)
        
            plt.sca(ax[j][i])
            
            plt.imshow(np.rot90(Z[j]), extent=[xmin[j], xmax[j], ymin[j], ymax[j]])
            plt.plot(values[0], values[1], 'k.', ms=1, alpha=0.1)
            plt.xlim(xmin[j], xmax[j])
            plt.ylim(ymin[j], ymax[j])
            plt.gca().set_aspect('auto')
        
        # save
        np.savez('../data/kde_{:s}'.format(labels[i]), X=X, Y=Y, Z=Z, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    
    plt.tight_layout()

def progenitors_kde():
    """Save smooth distributions for each progenitors"""
    
    # h3 giants
    th = Table.read('../data/rcat_substructure.fits')
    print(np.array(np.unique(th['Substructure_ID'])))
    
    prog_ids = ['GSE', 'HS', 'Iitoi', 'Sequo', 'Arjun', 'Sgr', 'Thamn', 'Wuk']
    labels = ['gse', 'helmi', 'iitoi', 'sequoia', 'arjuna', 'sgr', 'thamnos', 'wukong', 'cetus']
    
    cx = ['Lz', 'Lz', 'Lz', 'Lx']
    cy = ['E_tot_pot1', 'Lperp', 'Ly', 'Ly']
    
    # set up boundaries of the 2D output array
    lz_min, lz_max = -4.5, 4.5
    ly_min, ly_max = -7, 7
    lx_min, lx_max = -5.5, 5.5
    lperp_min, lperp_max = 0, 7
    e_min, e_max = -0.18, -0.03
    
    xmin = np.array([lz_min, lz_min, lz_min, lx_min])
    xmax = np.array([lz_max, lz_max, lz_max, lx_max])
    ymin = np.array([e_min, lperp_min, ly_min, ly_min])
    ymax = np.array([e_max, lperp_max, ly_max, ly_max])
    
    # set up output arrays
    X1, Y1 = np.mgrid[lz_min:lz_max:200j, e_min:e_max:200j]
    X2, Y2 = np.mgrid[lz_min:lz_max:200j, lperp_min:lperp_max:200j]
    X3, Y3 = np.mgrid[lz_min:lz_max:200j, ly_min:ly_max:200j]
    X4, Y4 = np.mgrid[lx_min:lx_max:200j, ly_min:ly_max:200j]
    X = np.array([X1, X2, X3, X4])
    Y = np.array([Y1, Y2, Y3, Y4])
    
    positions1 = np.vstack([X1.ravel(), Y1.ravel()])
    positions2 = np.vstack([X2.ravel(), Y2.ravel()])
    positions3 = np.vstack([X3.ravel(), Y3.ravel()])
    positions4 = np.vstack([X4.ravel(), Y4.ravel()])
    positions = np.array([positions1, positions2, positions3, positions4])
    
    # figure sizing
    nrow = 4
    ncol = len(labels)
    da = 2.5
    
    plt.close()
    fig, ax = plt.subplots(nrow, ncol, figsize=(ncol*da, nrow*da), sharex='row', sharey='row')
    
    for i in range(ncol-1,ncol):
        if labels[i]=='cetus':
            giants = Table.read('../data/rcat_giants.fits')
            ind = (giants['Ly']>0.7e3*u.kpc*u.km/u.s) & (giants['Lx']>2e3*u.kpc*u.km/u.s) & (giants['FeH']<-1.5) & (giants['Lz']<-1.4e3*u.kpc*u.km/u.s) & (giants['Lz']>-3.4e3*u.kpc*u.km/u.s) & (giants['E_tot_pot1']>-1e5*u.km**2*u.s**-2) &   (giants['E_tot_pot1']<-0.6e5*u.km**2*u.s**-2) & (giants['R_gal']>20) & (giants['eccen_pot1']<0.5) & (giants['RA']<150)
            t_ = giants[ind]
        else:
            ind = th['Substructure_ID']==prog_ids[i]
            t_ = th[ind]
        Z = np.zeros_like(X)
        
        for j in range(nrow):
            values = np.vstack([t_[cx[j]], t_[cy[j]]])
            kernel = gaussian_kde(values)
            Z[j] = np.reshape(kernel(positions[j]).T, X[j].shape)
        
            plt.sca(ax[j][i])
            
            plt.imshow(np.rot90(Z[j]), extent=[xmin[j], xmax[j], ymin[j], ymax[j]])
            plt.plot(values[0], values[1], 'k.', ms=1, alpha=0.1)
            plt.xlim(xmin[j], xmax[j])
            plt.ylim(ymin[j], ymax[j])
            plt.gca().set_aspect('auto')
        
        # save
        np.savez('../data/kde_{:s}'.format(labels[i]), X=X, Y=Y, Z=Z, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    
    plt.tight_layout()




def stream_gc_connections():
    """"""
    
    t = Table.read('/home/ana/data/baumgardt_positions.fits')
    
    #Omega Cen NGC 5139---Fimbulthul, NGC~3201---Gj\" oll, NGC~4590---Fjorm, NGC~5024---Sylgr and Ravi, NGC~5272---Sv\" ol, NGC~5824---Triangulum
    clusters = ['NGC 3201', 'NGC 4590', 'NGC 5024', 'NGC 5139', 'NGC 5272', 'NGC 5824']
    streams = ['gjoll', 'fjorm', 'sylgr', 'fimbulthul', 'svol', 'triangulum']
    direction = [-1, 1, -1, -1, 1, 1]
    nstep = [50, 100, 300, 70, 70, 600]
    N = len(clusters)
    
    dt = 0.5*u.Myr
    wangle = 180*u.deg
    l_off = 0*u.deg
    
    plt.close()
    fig = plt.figure(figsize=(12,5.2))
    ax = fig.add_subplot(111, projection='mollweide')
    
    for i in range(6):
        ind = t['Name']== clusters[i]
        t_ = t[ind]
        
        c = coord.SkyCoord(ra=t_['RAJ2000'], dec=t_['DEJ2000'], distance=t_['Rsun'], pm_ra_cosdec=t_['pmRA_'], pm_dec=t_['pmDE'], radial_velocity=t_['RV'], frame='icrs')
        cgal = c.transform_to(coord.Galactic)
        w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)[0]
        
        color = next(ax._get_lines.prop_cycler)['color']
        plt.plot((cgal.l+l_off).wrap_at(wangle).rad, cgal.b.rad, '+', color=color, mew=2, ms=10, label=t_['Name'][0])
        plt.text((cgal.l+l_off).wrap_at(wangle).rad, (cgal.b+3*u.deg).rad, t_['Name'][0], fontsize='xx-small', ha='center')

        orbit = ham.integrate_orbit(w0, dt=dt*direction[i], n_steps=nstep[i])
        orbit_gal = orbit.to_coord_frame(coord.Galactic, galactocentric_frame=gc_frame)
        plt.plot((orbit_gal.l+l_off).wrap_at(wangle).rad, orbit_gal.b.rad, '-', color=color, lw=1, label='')
        
        pkl = pickle.load(open('../data/streams/data_{:s}.pkl'.format(streams[i]), 'rb'))
        cs = coord.SkyCoord(ra=pkl['dec'][0], dec=pkl['dec'][1], frame='icrs')
        cs_gal = cs.transform_to(coord.Galactic)
        plt.plot((cs_gal.l+l_off).wrap_at(wangle).rad, cs_gal.b.rad, 'o', color=color, ms=8, label=streams[i])
        plt.text((cs_gal.l+l_off).wrap_at(wangle).rad[-1], cs_gal.b.rad[-1], streams[i], fontsize='xx-small')
    
    #plt.legend()
    plt.grid(ls=':')
    plt.xlabel('l [deg]')
    plt.ylabel('b [deg]')
    
    plt.tight_layout()

def stream_gc_connections_eq():
    """"""
    
    t = Table.read('/home/ana/data/baumgardt_positions.fits')
    
    #Omega Cen NGC 5139---Fimbulthul, NGC~3201---Gj\" oll, NGC~4590---Fjorm, NGC~5024---Sylgr and Ravi, NGC~5272---Sv\" ol, NGC~5824---Triangulum
    clusters = ['NGC 3201', 'NGC 4590', 'NGC 5024', 'NGC 5139', 'NGC 5272', 'NGC 5824']
    streams = [['gjoll'], ['fjorm'], ['sylgr', 'ravi'], ['fimbulthul'], ['svol'], ['triangulum', 'turbio']]
    direction = [[-1], [1], [-1, 1], [-1], [1], [1, 1]]
    nstep = [[50], [100], [300, 500], [70], [70], [600, 1]]
    N = len(clusters)
    
    dt = 0.5*u.Myr
    wangle = 180*u.deg
    ra_off = 142*u.deg
    #ra_off = 22*u.deg
    
    plt.close()
    fig = plt.figure(figsize=(12,5.2))
    ax = fig.add_subplot(111, projection='mollweide')
    
    for i in range(N):
        ind = t['Name']== clusters[i]
        t_ = t[ind]
        
        c = coord.SkyCoord(ra=t_['RAJ2000'], dec=t_['DEJ2000'], distance=t_['Rsun'], pm_ra_cosdec=t_['pmRA_'], pm_dec=t_['pmDE'], radial_velocity=t_['RV'], frame='icrs')
        w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)[0]
        
        color = next(ax._get_lines.prop_cycler)['color']
        plt.plot((c.ra+ra_off).wrap_at(wangle).rad, c.dec.rad, '+', color=color, mew=2, ms=10, label=t_['Name'][0])
        plt.text((c.ra+ra_off).wrap_at(wangle).rad, (c.dec+3*u.deg).rad, t_['Name'][0], fontsize='xx-small', ha='center')

        for j in range(len(direction[i])):
            orbit = ham.integrate_orbit(w0, dt=dt*direction[i][j], n_steps=nstep[i][j])
            orbit_eq = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=gc_frame)
            plt.plot((orbit_eq.ra+ra_off).wrap_at(wangle).rad, orbit_eq.dec.rad, '-', color=color, lw=1, label='')
            
            pkl = pickle.load(open('../data/streams/data_{:s}.pkl'.format(streams[i][j]), 'rb'))
            plt.plot((pkl['dec'][0]+ra_off).wrap_at(wangle).rad, pkl['dec'][1].rad, 'o', color=color, ms=8, label=streams[i][j])
            plt.text((pkl['dec'][0]+ra_off).wrap_at(wangle).rad[-1], pkl['dec'][1].rad[-1], streams[i][j], fontsize='xx-small')
    
    #plt.legend()
    plt.grid(ls=':')
    plt.xlabel('R.A. + {:.0f} [deg]'.format(ra_off.value))
    plt.ylabel('Dec [deg]')
    
    plt.tight_layout()
    plt.savefig('../plots/gc_streams_sky.png')


def stream_planes():
    """"""
    streams = [['gjoll', 'leiptr', 'phlegethon'], ['gd1', 'wambelong','ylgr']]
    N = len(streams)
    colors = ['r', 'b']
    
    wangle = 180*u.deg
    ra_off = -90*u.deg
    
    plt.close()
    fig = plt.figure(figsize=(12,5.2))
    ax = fig.add_subplot(111, projection='mollweide')
    
    for i in range(N):
        for j in range(len(streams[i])):
            pkl = pickle.load(open('../data/streams/data_{:s}.pkl'.format(streams[i][j]), 'rb'))
            plt.plot((pkl['dec'][0]+ra_off).wrap_at(wangle).rad, pkl['dec'][1].rad, 'o', color=colors[i], ms=8, label=streams[i][j])
            plt.text((pkl['dec'][0]+ra_off).wrap_at(wangle).rad[-1], pkl['dec'][1].rad[-1], streams[i][j], fontsize='xx-small')
    
    plt.grid(ls=':')
    plt.xlabel('R.A. + {:.0f} [deg]'.format(ra_off.value))
    plt.ylabel('Dec [deg]')
    
    plt.tight_layout()
    plt.savefig('../plots/retrograde_streams.png')
