"""Regenerate CF ternary plots with optimized parameters."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os, csv

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'output', 'figures')
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'output', 'data')

def to_cartesian(phi1, phi2, phi3):
    return phi1 + phi2*0.5, phi2*np.sqrt(3)/2

def draw_frame(ax):
    s3=np.sqrt(3); c=[(0,0),(1,0),(0.5,s3/2)]
    for i in range(3): ax.plot(*zip(c[i],c[(i+1)%3]),'-',color='#333',lw=2)
    for lv in np.arange(0.1,1.0,0.1):
        y=lv*s3/2; ax.plot([lv/2,1-lv/2],[y,y],'-',color='#ccc',lw=0.4,alpha=0.6)
        ax.plot([lv,lv/2],[0,y],'-',color='#ccc',lw=0.4,alpha=0.6)
        ax.plot([1-lv,1-lv/2],[0,y],'-',color='#ccc',lw=0.4,alpha=0.6)

def vlabels(ax,t,bl,br):
    s3=np.sqrt(3)
    ax.text(0.5,s3/2+0.06,t,ha='center',va='bottom',fontsize=13,fontweight='bold')
    ax.text(-0.06,-0.02,bl,ha='right',va='top',fontsize=13,fontweight='bold')
    ax.text(1.06,-0.02,br,ha='left',va='top',fontsize=13,fontweight='bold')

def load(tag):
    pairs=[]; spin=[]
    with open(os.path.join(DATA_DIR,f'binodal_{tag}.csv')) as f:
        for r in csv.DictReader(f):
            pairs.append((float(r['phi1_L8Bo']),float(r['phi2_Chloroform']),float(r['phi3_Donor'])))
    pl=[]; [pl.append((pairs[i],pairs[i+1])) for i in range(0,len(pairs)-1,2)]
    with open(os.path.join(DATA_DIR,f'spinodal_{tag}.csv')) as f:
        for r in csv.DictReader(f):
            spin.append((float(r['phi1']),float(r['phi2']),float(r['phi3'])))
    with open(os.path.join(DATA_DIR,f'critical_{tag}.csv')) as f:
        r=next(csv.DictReader(f))
        crit=(float(r['phi1']),float(r['phi2']),float(r['phi3']))
    return pl,spin,crit

SPECS={
    'PM6_CF':('PM6 / L8-Bo / CF','#4f46e5',('Chloroform','PM6','L8-Bo')),
    'D18_CF':('D18 / L8-Bo / CF','#f59e0b',('Chloroform','D18','L8-Bo')),
}

# 2x1 grid
fig,axes=plt.subplots(1,2,figsize=(16,7.5))
fig.suptitle('Chloroform Systems — Optimized Flory-Huggins Ternary Phase Diagrams',fontsize=15,fontweight='bold',y=0.98)
for idx,(tag,(title,color,vl)) in enumerate(SPECS.items()):
    ax=axes[idx]; pairs,spin,crit=load(tag)
    ax.set_aspect('equal'); ax.set_xlim(-0.1,1.1); ax.set_ylim(-0.08,np.sqrt(3)/2+0.14); ax.axis('off')
    draw_frame(ax)
    conc=sorted([p[0] for p in pairs],key=lambda x:x[2],reverse=True)
    xs,ys=zip(*[to_cartesian(p[0],p[1],p[2]) for p in conc])
    ax.plot(xs,ys,'-',color=color,lw=2.5,zorder=10)
    dil=sorted([p[1] for p in pairs],key=lambda x:x[2])
    xs2,ys2=zip(*[to_cartesian(p[0],p[1],p[2]) for p in dil])
    ax.plot(xs2,ys2,'-',color=color,lw=2.5,zorder=10)
    for i,(c,d) in enumerate(pairs):
        if i%6!=0: continue
        x1,y1=to_cartesian(c[0],c[1],c[2]); x2,y2=to_cartesian(d[0],d[1],d[2])
        ax.plot([x1,x2],[y1,y2],'-',color='#999',lw=0.5,alpha=0.35,zorder=5)
    valid=[(p[0],p[1],p[2]) for p in spin if 0<=p[0]<=1 and 0<=p[1]<=1 and 0<=p[2]<=1]
    valid.sort(key=lambda x:(x[1],x[0]))
    if len(valid)>2:
        xss,yss=zip(*[to_cartesian(p[0],p[1],p[2]) for p in valid])
        ax.plot(xss,yss,'--',color=color,lw=1.2,alpha=0.5,zorder=8)
    xc,yc=to_cartesian(crit[0],crit[1],crit[2])
    ax.scatter(xc,yc,marker='D',s=140,c=color,edgecolors='black',lw=1.8,zorder=20)
    ax.annotate(f'CP $\phi_2$={crit[1]:.3f}',(xc,yc),textcoords="offset points",xytext=(8,8),fontsize=9,color=color,fontweight='bold',bbox=dict(boxstyle='round,pad=0.3',fc='white',ec=color,alpha=0.85))
    vlabels(ax,vl[0],vl[1],vl[2])
    ax.set_title(f"{title}\n{len(pairs)} tie-lines | CP $\phi_2$={crit[1]:.3f}",fontsize=13,fontweight='bold',color=color)
plt.subplots_adjust(wspace=0.12,left=0.05,right=0.95)
fig.savefig(os.path.join(OUT_DIR,'ternary_real_CF_2x1.png'),dpi=200,bbox_inches='tight',facecolor='white')
print(f"Saved: ternary_real_CF_2x1.png")

# Single plots
for tag,(title,color,vl) in SPECS.items():
    fig,ax=plt.subplots(figsize=(10,10))
    pairs,spin,crit=load(tag)
    ax.set_aspect('equal'); ax.set_xlim(-0.12,1.12); ax.set_ylim(-0.1,np.sqrt(3)/2+0.15); ax.axis('off')
    draw_frame(ax)
    conc=sorted([p[0] for p in pairs],key=lambda x:x[2],reverse=True)
    xs,ys=zip(*[to_cartesian(p[0],p[1],p[2]) for p in conc])
    ax.plot(xs,ys,'-',color=color,lw=2.5,zorder=10)
    dil=sorted([p[1] for p in pairs],key=lambda x:x[2])
    xs2,ys2=zip(*[to_cartesian(p[0],p[1],p[2]) for p in dil])
    ax.plot(xs2,ys2,'-',color=color,lw=2.5,zorder=10)
    for i,(c,d) in enumerate(pairs):
        if i%4!=0: continue
        x1,y1=to_cartesian(c[0],c[1],c[2]); x2,y2=to_cartesian(d[0],d[1],d[2])
        ax.plot([x1,x2],[y1,y2],'-',color='#999',lw=0.5,alpha=0.35,zorder=5)
    valid=[(p[0],p[1],p[2]) for p in spin if 0<=p[0]<=1 and 0<=p[1]<=1 and 0<=p[2]<=1]
    valid.sort(key=lambda x:(x[1],x[0]))
    if len(valid)>2:
        xss,yss=zip(*[to_cartesian(p[0],p[1],p[2]) for p in valid])
        ax.plot(xss,yss,'--',color=color,lw=1.2,alpha=0.5,zorder=8)
    xc,yc=to_cartesian(crit[0],crit[1],crit[2])
    ax.scatter(xc,yc,marker='D',s=140,c=color,edgecolors='black',lw=1.8,zorder=20)
    vlabels(ax,vl[0],vl[1],vl[2])
    ax.set_title(f"{title}\nBinodal + Spinodal + {len(pairs)} Tie-lines + Critical Point",fontsize=15,fontweight='bold',color=color)
    leg=[Line2D([0],[0],color=color,lw=2.5,label='Binodal'),Line2D([0],[0],color=color,lw=1.2,ls='--',alpha=0.5,label='Spinodal'),Line2D([0],[0],color='#999',lw=0.6,alpha=0.5,label='Tie-lines'),Line2D([0],[0],marker='D',color='w',markerfacecolor=color,markersize=8,markeredgecolor='black',markeredgewidth=1.5,label='Critical Point')]
    ax.legend(handles=leg,loc='upper left',fontsize=11,framealpha=0.9)
    fig.savefig(os.path.join(OUT_DIR,f'ternary_real_{tag}.png'),dpi=200,bbox_inches='tight',facecolor='white')
    print(f"Saved: ternary_real_{tag}.png")
    plt.close(fig)

print("Done!")

