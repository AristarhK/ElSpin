"""Plot transient air plasma maps and histories. Requires numpy and matplotlib."""
from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', default='.')
    parser.add_argument('--prefix', default='')
    parser.add_argument('--out', default='plasma_results.png')
    parser.add_argument('--final-only',action='store_true',help='Show positive ions instead of the initial map.')
    args = parser.parse_args()
    root = Path(args.directory)
    status = (root/(args.prefix+'run_status.txt')).read_text()
    if 'PLASMA_TRANSIENT_COMPLETED' not in status:
        raise RuntimeError('The plasma transient did not finish; inspect run_status.txt.')
    final = np.genfromtxt(root/(args.prefix+'plasma_final.csv'),names=True,delimiter=',')
    initial = final if args.final_only else np.genfromtxt(root/(args.prefix+'plasma_initial.csv'),names=True,delimiter=',')
    history = np.genfromtxt(root/(args.prefix+'plasma_history.csv'),names=True,delimiter=',')
    r,z=np.unique(final['r_m']),np.unique(final['z_m'])
    shape=(len(z),len(r))
    keep_r,keep_z=r<=100e-6,z<=500e-6
    gas=final['gas_cell'].reshape(shape)[np.ix_(keep_z,keep_r)].T.astype(bool)
    def local(a,key):
        return np.ma.array(a[key].reshape(shape)[np.ix_(keep_z,keep_r)].T,mask=~gas)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(14,9),layout='constrained')
    fig.set_constrained_layout_pads(h_pad=.08,w_pad=.04,hspace=.12,wspace=.05)
    vmax=max(initial['ne_m3'].max(),final['ne_m3'].max())
    norm=LogNorm(vmin=max(1,vmax*1e-8),vmax=max(vmax,10))
    titles=['Концентрация положительных ионов' if args.final_only else 'Начальная концентрация электронов',
            f"Концентрация при t = {final['time_s'][0]*1e9:.3g} нс"]
    for ax,data,title,key,label in zip(axes[0,:2],[initial,final],titles,
                                    ['n_positive_m3' if args.final_only else 'ne_m3','ne_m3'],
                                    ['n₊, м⁻³' if args.final_only else 'nₑ, м⁻³','nₑ, м⁻³']):
        im=ax.pcolormesh(z[keep_z]*1e6,r[keep_r]*1e6,local(data,key),shading='auto',norm=norm,cmap='viridis')
        ax.set(title=title,xlabel='z, мкм',ylabel='r, мкм')
        fig.colorbar(im,ax=ax,label=label)
    im=axes[0,2].pcolormesh(z[keep_z]*1e6,r[keep_r]*1e6,local(final,'Te_eV'),shading='auto',cmap='inferno')
    axes[0,2].set(title='Эффективная температура электронов',xlabel='z, мкм',ylabel='r, мкм')
    fig.colorbar(im,ax=axes[0,2],label='Tₑ, эВ')
    time=history['time_s']*1e9
    axes[1,0].plot(time,history['ne_peak_m3'],lw=2,color='#176f9e')
    axes[1,0].set(xlabel='Время, нс',ylabel='Максимальная nₑ, м⁻³',title='Развитие заданной затравки')
    axes[1,1].plot(time,history['Te_peak_eV'],lw=2,color='#b25d20')
    axes[1,1].set(xlabel='Время, нс',ylabel='Максимальная Tₑ, эВ',title='Изменение энергии электронов')
    axes[1,2].plot(time,history['phi_feedback_relative']*100,lw=2,color='#4b8764')
    axes[1,2].set(xlabel='Время, нс',ylabel='Изменение потенциала / |U|, %',title='Обратное влияние заряда газа')
    for ax in axes[1,:]:
        ax.grid(alpha=.2)
    fig.suptitle(f'Воздух при 1 атм · сетка {len(r)}×{len(z)} · пример с заданной затравкой',fontsize=16,fontweight='bold')
    out=Path(args.out); out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out,dpi=170,facecolor='white')
    print(out.resolve())


if __name__=='__main__':
    main()
