import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
from scipy.optimize import curve_fit
import sys
import os

rcParams.update({'font.size': 26})
rcParams['text.usetex'] = True
#rc('text', usetex=True)

def func1(x,a,b):
    return  a*(1-b*x)  #

def func2(x,a,b,c):
    ep = x[0]
    es = x[1]
    return  a*(1-b*ep-c*es)

def Combine_files(fn_out):
    if os.path.exists(fn_out):
        os.remove(fn_out)
    filelist = [f for f in os.listdir(output_fldr) if f.startswith('Stab_') and f.endswith('.txt') and f != 'Stab_moon.txt']

    with open(output_fldr + fn_out, "w") as outfile:
        for fn in filelist:
            with open(output_fldr+fn,"r") as infile:
                outfile.write(infile.read())

def colorbar(mappable):
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax,orientation='horizontal')
    cbar.ax.xaxis.set_ticks_position('top')
    plt.sca(last_axes)
    return cbar

def plot_ep(ax,data_dir,data_file,lbl,idx,fname):
    pos = ax.figbox.get_points()
    cmap = cm.gnuplot_r
    vmin = 0.
    vmax = 1.
    my_cmap=cm.get_cmap(cmap)
    norm = colors.Normalize(vmin,vmax)
    cmmapable =cm.ScalarMappable(norm,my_cmap)
    cmmapable.set_array(range(0,1))
    cmap.set_under('w')

    data = np.genfromtxt(data_dir+data_file,delimiter=',',comments='#')
    stab = np.where(np.abs(data[:,-1]-1e5)<1e-6)[0]
    print(len(stab),len(data))

    xi = np.arange(0,0.51,0.01)
    yi = np.arange(0.25,0.55,0.01)
    X = []
    Y = []
    Z = []
    out = open(fname,'w')
    out.write("#e_p,R_H,f_stab\n")
    for xp in xi:
        for yp in yi:
            X.append(xp)
            Y.append(yp)
            xy_cut = np.where(np.logical_and(np.abs(data[:,2]-xp)<1e-6,np.abs(data[:,0]-yp)<1e-6))[0]
            stab_frac = 0.
            if len(xy_cut)>0:
                for cnt in range(0,len(xy_cut)):
                    if np.abs(data[xy_cut[cnt],-1]-1e5)<1e-6:
                        stab_frac += 1.
                    #else:
                    #    stab_frac -= 1.
                if stab_frac <= 0.05:
                    stab_frac -= 1.
                stab_frac /= float(len(xy_cut))
                #print(float(len(xy_cut)))
            else:
                stab_frac = -1.
            Z.append(stab_frac)
            out.write("%1.3f,%1.3f,%1.5f\n" % (xp,yp,stab_frac))
    out.close()
    X = np.asarray(X)
    Y = np.asarray(Y)
    Z = np.asarray(Z)

    Z[Z>=0.95] = 1.

    zi = griddata((X,Y),Z,(xi[None,:],yi[:,None]),method = 'nearest') #,fill_value=0
    CS = ax.pcolormesh(xi,yi,zi,cmap = cmap,vmin=vmin,vmax=vmax)
    x_data = []
    y_data = []
    prob_cut = 1.#float(sys.argv[1])
    for ecc in np.arange(0,0.5,0.01):
        e_cut = np.where(np.logical_and(np.abs(X-ecc)<1e-6,Z>=prob_cut))[0]
        if len(e_cut)>=1:
            x_data.append(ecc)
            #print(cut_y)
            y_data.append(np.max(Y[e_cut]))
            '''cut_y = Y[e_cut]
            cut_y[::-1].sort()
            if len(cut_y)>=2:
                y_data.append(np.mean(cut_y[:2]))
            else:
                y_data.append(np.mean(cut_y))'''
    popt,pcov = curve_fit(func1,x_data,y_data,method='lm') #bounds=([0.35,0.5],[0.55,1.5])
    print(popt,np.sqrt(np.diag(pcov)))
    sigma = np.sqrt(np.diag(pcov))
    ecc = np.arange(0,0.5,0.01)
    #ax.plot(ecc,func(ecc,0.4895,1.0305),'k--',lw=lw,label='Domingos et al. 2006')
    #ax.plot(ecc,func(ecc,*popt),'-',color='gray',lw=lw,label='This work')
    for i in range(0,300):
        a_i = np.random.normal(popt[0],sigma[0])
        b_i = np.random.normal(popt[1],sigma[1])
        y_i = func1(ecc,a_i,b_i)
        ax.plot(ecc,y_i,'-',color='gray',lw=3,alpha=0.1)
    if idx==1:
        ax.plot(ecc,func1(ecc,0.4895,1.0305),'k--',lw=lw,label='DWY06')
        ax.plot(ecc,func1(ecc,*popt),'-',color='r',lw=lw,label='This Work')
        ax.legend(loc='upper right',fontsize='x-large')
    else:
        ax.plot(ecc,func1(ecc,*popt),'-',color='r',lw=lw)
    #ax.legend(loc='best',fontsize=fs)
    if idx ==3:
        ax.set_xlabel(r"$e_{p}$",fontsize=fs)
    if idx ==1:
        ax.set_ylabel(r"$a_{%s}\;(R_{H,p})$" % lbl,fontsize=fs)
    else:
        ax.set_ylabel(r"$a_{%s}\;(R_{H,sat})$" % lbl,fontsize=fs)
    ax.minorticks_on()
    ax.tick_params(which='major',axis='both', direction='out',length = 10.0, width = 8.0,labelsize=fs)
    ax.tick_params(which='minor',axis='both', direction='out',length = 6.0, width = 8.0)
    ax.text(0.02,0.92,sublbl[idx-1], color='k',fontsize='xx-large',weight='bold',horizontalalignment='left',transform=ax.transAxes)

    ax.set_ylim(0.25,0.55)
    ax.set_xlim(0.,0.5)
    ax.set_yticks(np.arange(0.25,0.6,0.05))
    ax.set_xticks(np.arange(0,0.6,0.1))
    if idx ==1:
        color_label=r'$f_{stab}$'
        #cax = fig.add_axes([0.92,0.11,0.015,0.77])
        #cax=plt.axes([pos[1,0]-0.03,pos[0,1]-0.005,0.01,pos[1,1]-pos[0,1]+0.015])
        #cbar=plt.colorbar(cmmapable,cax=cax,orientation='horizontal')
        #cbar=fig.colorbar(CS,ax=ax,orientation='horizontal')
        cbar=colorbar(CS)
        #cbar.set_label(color_label,fontsize=fs)
        cbar.set_ticks(np.arange(0.,1.25,0.25))
        cbar.ax.tick_params(axis='both', direction='out',length = 8.0, width = 8.0,labelsize=fs)
        fig.text(0.3,0.94,color_label, color='black',fontsize=fs,ha='center', va='center')

def plot_esat(ax,data_dir,idx,fname):
    cmap = cm.jet #cm.nipy_spectral_r
    pos = ax.figbox.get_points()
    vmin = 0.1
    vmax = 0.4
    my_cmap=cm.get_cmap(cmap)
    norm = colors.Normalize(vmin,vmax)
    cmmapable =cm.ScalarMappable(norm,my_cmap)
    cmmapable.set_array(range(0,1))
    cmap.set_under('gray')

    xi = np.arange(0.0,0.51,0.01)
    yi = np.arange(0.0,0.51,0.01)
    X = []
    Y = []
    Z = []
    stab_num = 19
    if idx == 4:
        stab_num = 19
    out = open(fname,'w')
    out.write("#e_p,e_sat,a_crit\n")
    for xp in xi:
        for yp in yi:
            a_sat = -1.
            fname = "MaxEcc_[%1.3f,%1.3f].txt" % (xp,yp)
             
            if os.path.exists(data_dir+fname):
                with open(data_dir+fname,'r') as f:
                    lines = f.readlines()
                if len(lines)>1:
                    data = np.genfromtxt(data_dir+fname,delimiter=',',comments='#')
                    if len(data.shape)>1:
                        a_cut = np.where(np.abs(data[:,0]-np.min(data[:,0]))<1e-6)[0]
                        if len(a_cut) >= stab_num:
                            a_sat = np.min(data[:,0])
                            #if idx == 4:
                            #    a_sat = func2(np.array([xp,yp]),0.31662589, 0.23965448, 1.02990124)
                                #print(a_sat)
            if a_sat > 0:
                X.append(xp)
                Y.append(yp)
                Z.append(a_sat)
                out.write("%1.3f,%1.3f,%1.5f\n" % (xp,yp,a_sat))
    out.close()
    X = np.asarray(X)
    Y = np.asarray(Y)
    Z = np.asarray(Z)

    zi = griddata((X,Y),Z,(xi[None,:],yi[:,None]),method = 'linear',fill_value=-1)
    stab = np.where(Z>0)[0]
    x_data = [X[stab],Y[stab]]
    Z = Z[stab]
    print(np.min(Z))
    popt,pcov = curve_fit(func2,x_data,Z,method='lm') #bounds=([0.35,0.5],[0.55,1.5])
    print(popt,np.sqrt(np.diag(pcov)))
    CS = ax.pcolormesh(xi,yi,zi,cmap = cmap,vmin=vmin,vmax=vmax)
    
    if idx==4:
        ax.set_xlabel(r"$e_{p}$",fontsize=fs)
    ax.set_ylabel(r"$e_{sat}$",fontsize=fs)
    ax.minorticks_on()
    ax.tick_params(which='major',axis='both', direction='out',length = 10.0, width = 8.0,labelsize=fs)
    ax.tick_params(which='minor',axis='both', direction='out',length = 6.0, width = 8.0)
    ax.set_ylim(0.0,0.5)
    ax.set_xlim(0.0,0.5)
    ax.set_yticks(np.arange(0,0.6,0.1))
    ax.set_xticks(np.arange(0,0.6,0.1))
    ax.text(0.02,0.92,sublbl[idx-1], color='w',fontsize='xx-large',weight='bold',horizontalalignment='left',transform=ax.transAxes)
    #ax.text(0.5,1.025,sub_label[i],color='k',fontsize=fs,horizontalalignment='center',transform=ax_list[i].transAxes)
    if idx == 2:
        color_label=r'$a_{\rm crit}\;(R_H)$'
        #cax = fig.add_axes([0.92,0.11,0.015,0.77])
        #cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]])
        #cbar= plt.colorbar(cmmapable,cax=cax,orientation='vertical')
        cbar=colorbar(CS)
        #cbar.set_label(color_label,fontsize=fs)
        cbar.set_ticks(np.arange(0.1,0.45,0.05))
        cbar.ax.tick_params(axis='both', direction='out',length = 8.0, width = 8.0,labelsize=fs)
        fig.text(0.7,0.94,color_label, color='black',fontsize=fs,ha='center', va='center')


home = os.getcwd() + "/"
output_fldr = home + "Submoon_circ/"

Combine_files('Stab_moon.txt')
fs = 'xx-large'
width = 10.
aspect = 2.1
ms = 6.5
lw=5
sublbl = [r'$\textbf{a}$',r'$\textbf{b}$',r'$\textbf{c}$',r'$\textbf{d}$']


fig = plt.figure(figsize=(aspect*width,2.*width),dpi=300)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)


plot_ep(ax1,home+"Domingos06_circ/","Stab_moon.txt","sat",1,"Stab_frac_moon.txt")
plot_esat(ax2,home+"GenRuns_out/",2,"Contour_moon.txt")
plot_ep(ax3,home+"Submoon_circ/","Stab_moon.txt","sub",3,"Stab_frac_submoon.txt")
plot_esat(ax4,home+"Submoon_out/",4,"Contour_submoon.txt")

fig.subplots_adjust(wspace=0.25,hspace=0.15)
fig.savefig("Stability_Fig.png",bbox_inches='tight',dpi=300)
plt.close()
