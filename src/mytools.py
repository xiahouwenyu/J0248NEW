from astropy.io import fits
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
try:
    import tevcat as TeVCat
    haveTeVCat = True
except ImportError as e:
    haveTeVCat = False
    print(e)
import MapPalette
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.coordinates import FK5
import threeML
import gammapy as gp
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.simbad import Simbad
from astropy.coordinates import Angle

try:
    tevcat = TeVCat.TeVCat()
except IOError as e:
    print(e)
    print("Downloading data from tevcat.uchicago.edu")
    tevcat = TeVCat.TeVCat()
except:
    print("Why caught here?")
    print("Downloading data from tevcat.uchicago.edu")
    tevcat = TeVCat.TeVCat()

def Drawgascontour():
        from matplotlib.colors import Normalize
        with fits.open('./J0248_co_-55--30_all.fits') as hdul:
                # 输出文件信息
                qtsj = hdul[0].data
                hdul[0].header
        hdul[0].header
        # 获取坐标信息
        crval1 = hdul[0].header['CRVAL1']
        crval2 = hdul[0].header['CRVAL2']
        cdelt1 = hdul[0].header['CDELT1']
        cdelt2 = hdul[0].header['CDELT2']
        crpix1 = hdul[0].header['CRPIX1']
        crpix2 = hdul[0].header['CRPIX2']

        # 计算x轴和y轴的坐标范围
        naxis1, naxis2 = qtsj.shape
        xmin = (1 - crpix1) * cdelt1 + crval1
        xmax = (naxis1 - crpix1) * cdelt1 + crval1
        ymin = (1 - crpix2) * cdelt2 + crval2
        ymax = (naxis2 - crpix2) * cdelt2 + crval2
        glon = np.linspace(xmin,xmax,naxis1)
        glat = np.linspace(ymin,ymax,naxis2)
        Glon,Glat = np.meshgrid(glon,glat)
        galactic_coord = SkyCoord(Glon* u.degree, Glat* u.degree, frame='galactic')
        j2000_coords = galactic_coord.transform_to('fk5')
        Glon,Glat = j2000_coords.ra,j2000_coords.dec
        plt.contour(Glon,Glat,qtsj,5,cmap="Greys",levels=np.array([0.2,0.3,0.5,0.7,1,1.5,2,3,4])*1e22,norm=Normalize(vmin=0.2e22,vmax=1e22),alpha=0.7)

def GetTeVcat(xmin,xmax,ymin,ymax):
    xa=[]
    ya=[]
    assoca=[]
    sources = tevcat.getSources()
    for i in range(0,len(sources)):
        sourceFK5 = sources[i].getFK5()
        ras=sourceFK5.ra.degree
        decs=sourceFK5.dec.degree
        assoc = sources[i].getCanonicalName()
        if (ras < xmin) or (ras > xmax) or (decs < ymin) or (decs>ymax):
            continue
        if assoc in assoca:
            continue
        xa.append(ras)			
        ya.append(decs)			
        assoca.append(assoc)
    sources_tmp = list(zip(xa,ya,assoca))
    sources_tmp.sort(key=lambda source: source[0])
    return sources_tmp

def GetFermicat(xmin,xmax,ymin,ymax):
    # 打开FITS文件
    header = []
    data = []
    with fits.open('./gll_psch_v13.fit') as hdul:

        # 输出文件信息
        hdul.info()

        # 输出每个HDU的信息
        for i, hdu in enumerate(hdul):
            print(f'HDU {i}:')
            header.append(hdu.header)
            data.append(hdu.data)
    
    xa=[]
    ya=[]
    assoca=[]
    sources = tevcat.getSources()
    for i in range(0,len(data[1])):
        ras=data[1][i][1]
        decs=data[1][i][2]
        assoc = data[1][i][0]
        if (ras < xmin) or (ras > xmax) or (decs < ymin) or (decs>ymax):
            continue
        if assoc in assoca:
            continue
        xa.append(ras)			
        ya.append(decs)			
        assoca.append(assoc)
    sources_tmp = list(zip(xa,ya,assoca))
    sources_tmp.sort(key=lambda source: source[0])
    return sources_tmp

def GetPSRcat(xmin,xmax,ymin,ymax):
    # 设置max_catalogs参数为100
    Vizier.ROW_LIMIT = -1  # 无限制
    Vizier.MAX_RESULTS = -1 # 无限制
    Vizier.TIMEOUT = 180 # 设置超时时间
    xa=[]
    ya=[]
    assoca=[]
    # 获取ATNF pulsar目录的数据
    atnf_catalog = Vizier.get_catalogs('B/psr')[0]
    for i in range(0,len(atnf_catalog)):
        try:
            ras=float(Angle(atnf_catalog[i][2], unit='hourangle').degree)
            decs=float(Angle(atnf_catalog[i][3], unit='degree').value)
        except:
            continue
        assoc = atnf_catalog[i][0]
        if (ras < xmin) or (ras > xmax) or (decs < ymin) or (decs>ymax):
            continue
        if assoc in assoca:
            continue
        xa.append(ras)			
        ya.append(decs)			
        assoca.append(assoc)
    sources_tmp = list(zip(xa,ya,assoca))
    sources_tmp.sort(key=lambda source: source[0])
    return sources_tmp

def GetSNRcat(xmin,xmax,ymin,ymax):
    # 设置max_catalogs参数为100
    Vizier.ROW_LIMIT = -1  # 无限制
    Vizier.MAX_RESULTS = -1 # 无限制
    Vizier.TIMEOUT = 180 # 设置超时时间
    xa=[]
    ya=[]
    assoca=[]
    # 获取ATNF pulsar目录的数据
    green_catalog = Vizier.get_catalogs('VII/278')[0]
    for i in range(0,len(green_catalog)):
        try:
            ras=float(Angle(green_catalog[i][1], unit='hourangle').degree)
            decs=float(Angle(green_catalog[i][2], unit='degree').value)
        except:
            continue
        assoc = green_catalog[i][0]
        if (ras < xmin) or (ras > xmax) or (decs < ymin) or (decs>ymax):
            continue
        if assoc in assoca:
            continue
        xa.append(ras)			
        ya.append(decs)			
        assoca.append(assoc)
    sources_tmp = list(zip(xa,ya,assoca))
    sources_tmp.sort(key=lambda source: source[0])
    return sources_tmp

def Drawcat(xmin,xmax,ymin,ymax,cat="TeV",mark="s",c="black",angle=45):
    if cat=="TeV":
        sources_tmp = GetTeVcat(xmin,xmax,ymin,ymax)
    elif cat=="fermi":
        sources_tmp = GetFermicat(xmin,xmax,ymin,ymax)
    elif cat=="psr":
        sources_tmp = GetPSRcat(xmin,xmax,ymin,ymax)
    elif cat=="snr":
        sources_tmp = GetSNRcat(xmin,xmax,ymin,ymax)
    ymid = np.mean([ymin,ymax])
    i=0
    dt=0
    rt=0
    pre_rt1=0
    pre_rt2=0
    dr=(xmax-xmin)/2/(len(sources_tmp)/2.+1)/2
    for r, d, s in sources_tmp:
            print(r, d, s)
            
            if d>ymid:
    #                        if np.abs(r-pre_rt1) <dr:
                if np.abs(r)-np.abs(pre_rt1) <dr:
                    rt=pre_rt1+dr
                else :
                    rt=r
                pre_rt1=rt	
                dt=(ymax+ymid)/2.
                Rotation=angle
                Va='bottom'
            else:
    #                        if np.abs(r-pre_rt2) <dr:
                if np.abs(r)-np.abs(pre_rt2) <dr:
                    rt=pre_rt2+dr
                else :
                    rt=r
                pre_rt2=rt	
                dt=(ymid+ymin)/2.
                Rotation=360-angle
                Va='top'
            i+=1
            plt.text(rt,dt, s+'', color=c,
                    rotation=Rotation, #catLabelsAngle,
                    va=Va,
                    fontdict={'family': 'sans-serif',
                                'size': 7,
                                'weight': 'bold'})
            plt.plot([r,rt],[d,dt],'k--',c=c)
            plt.scatter(r,d, color=c, facecolors="none", marker=mark)

def interpimg(hp_map,xmin,xmax,ymin,ymax,xsize):
    faspect = abs(xmax - xmin)/abs(ymax-ymin)
    phi   = np.linspace(xmin, xmax, xsize)
    theta = np.linspace(ymin, ymax, int(xsize/faspect))
    Phi, Theta = np.meshgrid(phi, theta)
    rotimg = hp.get_interp_val(hp_map, Phi,Theta,lonlat=True) #,nest=True
    # plt.contourf(Phi,Theta,rotimg)
    # plt.imshow(rotimg, origin="lower",extent=[xmin,xmax,ymin,ymax])
    # plt.colorbar()
    return rotimg

def hpDraw(map, ra, dec, rad=5, radx=5,rady=2.5,contours=[3,5],colorlabel="Excess",color="Fermi"):
    ra_crab, dec_crab =  83.63,22.02
    ra_j0248, dec_j0248 = 42.19,60.35
    # ra,dec = ra_crab,dec_crab
    # ra,dec = ra_j0248, dec_j0248
    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_subplot(111, projection='mollweide')
    # rad = 5
    # radx = 5
    # rady = 2.5


    ymax = dec+rady/2
    ymin = dec-rady/2
    xmin = ra-radx/2
    xmax = ra+radx/2

    tfig   = plt.figure(num=2)
    rot = (0, 0, 0)
    coord = 'C'
    xsize = 2048
    # img = hp.cartview(hp_map,fig=2,lonra=[ra-rad,ra+rad],latra=[dec-rad,dec+rad],return_projected_map=True, rot=rot, coord=coord, xsize=xsize)
    img = interpimg(map, xmin,xmax,ymin,ymax,xsize)
    plt.close(tfig)


    plt.figure(dpi=200)
    dMin = -5
    dMax = 15
    dMin = np.min(img) if np.min(img) != None else -5
    dMax = np.max(img) if np.max(img) != None else 15
    if color == "Milagro":
        textcolor, colormap = MapPalette.setupMilagroColormap(dMin-1, dMax+1, 3, 100)
    elif color == "Fermi":
        textcolor, colormap = MapPalette.setupGammaColormap(10000)
    plt.imshow(img, origin="lower",extent=[xmin,xmax,ymin,ymax], cmap=colormap) #


    plt.grid(linestyle="--")
    cbar = plt.colorbar(format='%.2f',orientation="horizontal",shrink=0.6,
                            fraction=0.1,
                            #aspect=25,
                            pad=0.15)

    cbar.set_label(colorlabel)
    cbar.set_ticks(np.concatenate(([np.min(img)],[np.mean(img)],[np.max(img)]))) #,cbar.get_ticks()

    contp = plt.contour(img,levels=np.sort(contours),colors='g',linestyles = '-',linewidths = 2,origin='upper',extent=[xmin, xmax, ymax, ymin])
    fmt = {}
    strs=[]
    for i in range(len(contours)):
        strs.append('%d$\sigma$'%(contours[i]))
    for l, s in zip(contp.levels, strs):
        fmt[l] = s

    CLabel = plt.clabel(contp, contp.levels, use_clabeltext=True, rightside_up=True, inline=1, fmt=fmt, fontsize=10)

    for l in CLabel:
        l.set_rotation(180)


    plt.xlabel(r"$\alpha$ [$^\circ$]")
    plt.ylabel(r"$\delta$ [$^\circ$]")

    plt.xlim(ra+radx/2,ra-radx/2)
    plt.ylim(dec-rady/2,dec+rady/2)
    Drawgascontour()

    plt.gca().set_aspect(1./np.cos((ymax+ymin)/2*np.pi/180))
    plt.scatter(ra, dec, s=20**2,marker="+", facecolor="#000000", color="#000000")
    Drawcat(xmin,xmax,ymin,ymax,"TeV","s","black",60)
    Drawcat(xmin,xmax,ymin,ymax,"psr","*","black",60)
    Drawcat(xmin,xmax,ymin,ymax,"snr","o","black",60)
    plt.savefig("J0248_sig_llh.png",dpi=300)
