#%%



def machplot(ylim = [0.1,10],fontsize = 8):
    # plt.plot(ma, color = label = 'No Avg')
    fig,ax = plt.subplots(1,2,figsize = (10,5))
    # ax[0].plot(ma, label = 'Very Short Avg',color='k')
    ax[0].plot(mameans[0], label = 'Short Avg',color='darkblue')
    ax[0].plot(mameans[1],label = 'Medium Avg',color='limegreen')
    ax[0].plot(mameans[2],label = 'Long Avg',color = 'orange')
    ax[0].plot(mameans[3],label = 'Very Long Avg',color = 'r')
    ax[0].axhline(y=1,color='k', label = '$R_A$')
    ax[0].axhline(y=1.6,color='k',linestyle='dashed', label = 'CL limits')
    ax[0].axhline(y=0.6,color='k',linestyle='dashed')
    ax[0].set_ylim(ylim[0],ylim[1])
    ax[0].set_title('Alfvén Mach number ($M_A$)')
    ax[0].legend(fontsize = fontsize)
    ax[0].set_yscale('log')

    ax[1].plot(abs(dva[0]/vameans[0]), label = 'Veru Short Scales (12 s)',color='darkblue')
    ax[1].plot(abs(dva[1]/vameans[1]),label = 'Short Scales (2 min)',color='limegreen')
    ax[1].plot(abs(dva[2]/vameans[2]),label = 'Medium Scales (20 min)',color = 'r')
    ax[1].set_ylim(0,1)
    ax[1].axhline(y=0,color='k')
    ax[1].set_title('$\delta v_A/v_A$')
    ax[1].legend(fontsize = fontsize)
    # ax[1].set_yscale('log')
machplot(ylim=[0.1,10])
#%%

def vaplot(ylim = [0,500],xlim = [0,len(vameans[0])], fontsize=10):
    # plt.plot(ma, color = label = 'No Avg')
    fig,ax = plt.subplots(4,1,figsize = (8,15))
    ax[0].plot(mameans[0],label = 'Small Scales (~34 $Mm$)',color='darkblue')
    ax[0].plot(mameans[1],label = 'Medium Scales (~0.5 $R_s$)',color='limegreen')
    ax[0].plot(mameans[2],label = 'Large Scales (~5 $R_s$ )',color = 'orange')
    ax[0].plot(mameans[3],label = 'Very Large Scales (~50 $R_s$)',color ='r')
    ax[0].axhline(y=1,color='k', label = '$R_A$')
    ax[0].axhline(y=1.6,color='k',linestyle='dashed', label = 'CL limits')
    ax[0].axhline(y=0.6,color='k',linestyle='dashed')
    ax[0].set_ylim(0.1,10)
    ax[0].set_xlim(xlim[0],xlim[1])
    ax[0].set_ylabel('Alfvén Mach number ($M_A$)')
    ax[0].legend(fontsize = fontsize,loc = 'upper left')
    ax[0].set_yscale('log')

    # ax[1].plot(1e-3*abs(dva[0]), label = 'Very Short Scales (12 s)',color='grey')
    # ax[1].plot(1e-3*abs(dva[1]),label = 'Short Scales (2 min)',color='darkblue')
    # ax[1].plot(1e-3*abs(dva[2]),label = 'Medium Scales (20 min)',color = 'limegreen')
    # ax[1].plot(1e-3*abs(dva[3]),label = 'Large Scales (3.3 hr)',color = 'orange')
    # ax[1].set_ylim(ylim[0],200)
    # ax[1].set_ylabel('$|\delta v_A|$')
    # ax[1].legend(fontsize = fontsize)
    
    ax[1].plot(abs(dva[0]/vameans[0]),label = 'Very Small Scales (~3.5 $Mm$) ',color='grey')
    ax[1].plot(abs(dva[1]/vameans[1]),label = 'Small Scales (~34 $Mm$)',color='darkblue')
    ax[1].plot(abs(dva[2]/vameans[2]),label = 'Medium Scales ( ~0.5 $R_s$)',color='limegreen')
    ax[1].plot(abs(dva[3]/vameans[3]),label = 'Large Scales (~5 $R_s$ )',color = 'orange')
    ax[1].set_ylim(0,1)
    ax[1].set_xlim(xlim[0],xlim[1])
    ax[1].set_ylabel('$|\delta v_A/v_A|$')
    ax[1].legend(fontsize = 10,loc='upper left')

    ax[2].plot(abs(dvth[0]/vthmeans[0]),label = 'Very Small Scales (~3.5 $Mm$) ',color='grey')
    ax[2].plot(abs(dvth[1]/vthmeans[1]),label = 'Small Scales (~34 $Mm$)',color='darkblue')
    ax[2].plot(abs(dvth[2]/vthmeans[2]),label = 'Medium Scales ( ~0.5 $R_s$)',color='limegreen')
    ax[2].plot(abs(dvth[3]/vthmeans[3]),label = 'Large Scales (~5 $R_s$ )',color = 'orange')
    ax[2].set_ylim(0,0.5)
    ax[2].set_xlim(xlim[0],xlim[1])
    ax[2].set_ylabel(r'$|\delta v_{th}/v_{th}|$')
    ax[2].legend(fontsize = 10,loc = 'upper left')

    ax[3].plot(lperp_over_lpar[0], label = 'Very Short Scales (~ 3.5 Mm)',color='grey')
    ax[3].plot(lperp_over_lpar[1], label = 'Short Scales (~ 34 Mm)',color='darkblue')
    ax[3].plot(lperp_over_lpar[2],label = 'Medium Scales (~ 0.5 $R_s$ m)',color='limegreen')
    ax[3].plot(lperp_over_lpar[3],label = 'Large Scales (~ 5 $R_s$)',color = 'orange')
    ax[3].plot(lperp_over_lpar[4],label = 'Very Large Scales (~ 5 $R_s$)',color = 'r')
    # ax[3].plot(lperp_over_lpar[4],label = 'Very Large Scales (~50 $R_s$)',color = 'r')
    ax[3].axhline(y=1,color = 'k')
    ax[3].set_ylim(1e-2,1e2)
    ax[3].set_ylabel(r'$\ell_\perp / \ell_\parallel$',fontsize=fontsize)
    # [3]0].set_ylabel(r'$\nabla_\perp$',loc = 'top',rotation =90,fontsize=fontsize)
    ax[3].set_yscale('log')
    # ax[3].set_title(r'$\ell_\perp / \ell_\parallel$')
    ax[3].legend(fontsize = 10)

vaplot(fontsize=5)



#%%


def lplot(ylim = [0.001,1000],fontsize = 10):
    # plt.plot(ma, color = label = 'No Avg')
    fig,ax = plt.subplots(1,1,figsize = (6,5))
    # ax.plot(1/np.array(lperp_over_lpar), label = 'Very Short Avg',color='k')
    ax.plot(lpar_over_l[0]**2, label = 'Very Short Scales (~ 3.5 Mm)',color='grey')
    ax.plot(lpar_over_l[1]**2, label = 'Short Scales (~ 34 Mm)',color='darkblue')
    ax.plot(lpar_over_l[2]**2,label = 'Medium Scales (~ 0.5 $R_s$ m)',color='limegreen')
    ax.plot(lpar_over_l[3]**2,label = 'Large Scales (~ 5 $R_s$)',color = 'orange')
    ax.plot(lpar_over_l[4]**2,label = 'Very Large Scales (~50 $R_s$)',color = 'r')
    ax.axhline(y=1,color = 'k')
    ax.set_ylim(ylim[0],ylim[1])
    # 0].set_ylabel(r'$\nabla_\perp$',loc = 'bottom',rotation ='horizontal',fontsize=fontsize)
    # 0].set_ylabel(r'$\nabla_\perp$',loc = 'top',rotation =90,fontsize=fontsize)
    ax.set_yscale('linear')
    ax.set_title(r'$\ell_\perp^2/ \ell^2$')
    ax.legend(fontsize = fontsize)

    # ax[1].set_yscale('log')




# def lplot(ylim = [0.001,1000],fontsize = 10):
#     # plt.plot(ma, color = label = 'No Avg')
#     fig,ax = plt.subplots(1,1,figsize = (6,5))
#     # ax.plot(1/np.array(lperp_over_lpar), label = 'Very Short Avg',color='k')
#     ax.plot(lperp_over_l[0], label = 'Very Short Scales (~ 3.5 Mm)',color='grey')
#     ax.plot(lperp_over_l[1], label = 'Short Scales (~ 34 Mm)',color='darkblue')
#     ax.plot(lperp_over_l[2],label = 'Medium Scales (~ 0.5 $R_s$ m)',color='limegreen')
#     ax.plot(lperp_over_l[3],label = 'Large Scales (~ 5 $R_s$)',color = 'orange')
#     ax.plot(lperp_over_l[4],label = 'Very Large Scales (~50 $R_s$)',color = 'r')
#     ax.axhline(y=1,color = 'k')
#     ax.set_ylim(ylim[0],ylim[1])
#     ax.set_ylabel(r'$\ell_\perp/\ell$',fontsize=fontsize)
#     # 0].set_ylabel(r'$\nabla_\perp$',loc = 'top',rotation =90,fontsize=fontsize)
#     ax.set_yscale('log')
#     # ax.set_title(r'$\nabla_\perp / \nabla_\parallel$')
#     ax.legend(fontsize = fontsize)

#     # ax[1].set_yscale('log')
lplot(ylim = [1e-1,1])

# %%

#%%


def nplot(ylim = [1e-5,1e4],fontsize = 8):
    fig,ax = plt.subplots(1,3,figsize = (12,4))
    ax[0].plot(1e-6*nmeans[0], label = 'Short Avg',color='darkblue')
    ax[0].plot(1e-6*nmeans[1],label = 'Medium Avg',color='limegreen')
    ax[0].plot(1e-6*nmeans[2],label = 'Long Avg',color = 'r')
    ax[0].set_ylim(ylim[0],ylim[1])
    ax[0].set_yscale('log')
    ax[0].set_title('Density ($n_i$)')
    ax[0].legend(fontsize = fontsize)

    ax[1].plot(1e-6*abs(dn[0]), label = 'Short Avg',color='darkblue')
    ax[1].plot(1e-6*abs(dn[1]),label = 'Medium Avg',color='limegreen')
    ax[1].plot(1e-6*abs(dn[2]),label = 'Long Avg',color = 'r')
    ax[1].set_ylim(ylim[0],5000)
    ax[1].set_title('$|\delta n_i|$')
    ax[1].set_yscale('log')
    ax[1].legend(fontsize = fontsize)
    
    ax[2].plot(abs(dn[0]/nmeans[0]), label = 'Short Avg',color='darkblue')
    ax[2].plot(abs(dn[1]/nmeans[1]),label = 'Medium Avg',color='limegreen')
    ax[2].plot(abs(dn[2]/nmeans[2]),label = 'Long Avg',color = 'r')
    ax[2].set_ylim(0.1,5)
    ax[2].set_yscale('log')
    ax[2].set_title('$|\delta n_i/n_i|$')
    ax[2].legend(fontsize = fontsize)
    # ax[1].set_yscale('log')
    # ax[1].set_yscale('log')

nplot()
# %%

plt.plot(1.602e19*Tperp[0], label = 'Short Avg',color='darkblue')
plt.plot(1.602e19*Tperp[1],label = 'Medium Avg',color='limegreen')
plt.plot(1.602e19*Tperp[2],label = 'Long Avg',color = 'r')
# plt.ylim(ylim[0],ylim[1])
plt.yscale('log')
plt.title('Density ($n_i$)')
plt.legend(fontsize =10)
#%%

plt.scatter(abs(dva[0]/vameans[0]),Tparperp[0])
# plt.scatter(abs(dva[1]/vameans[1]),Tparperp[1])
# plt.scatter(abs(dva[2]/vameans[2]),Tparperp[2])
# plt.scatter(abs(dva[3]/vameans[3]),Tparperp[3])
plt.yscale('log')
plt.axhline(y=1,color='k',linestyle='dashed')
plt.ylabel('$T_\perp/T_\parallel$')
plt.xlabel('$|\delta v_A/v_A|$')
#%%
# %%