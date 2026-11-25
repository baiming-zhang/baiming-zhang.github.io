"""Render the stored full-grid VOF fields used by the operator benchmark."""
from pathlib import Path
import json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from large_grid_operator_fixtures import R,save

def main():
    folder=R/'figures';folder.mkdir(exist_ok=True)
    if (folder/'field_figure_manifest.json').exists() and (folder/'field_figure_manifest.json').stat().st_mtime>=Path(__file__).stat().st_mtime:return
    started=time.perf_counter()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.linewidth':.65,'pdf.fonttype':42,'xtick.major.width':.6,'ytick.major.width':.6})
    palette=LinearSegmentedColormap.from_list('water_fraction',['#ffffff','#287abd'])
    fig,axes=plt.subplots(1,2,figsize=(6.6,3.15))
    fig.subplots_adjust(left=.075,right=.99,bottom=.17,top=.91,wspace=.28)
    for i,(ax,scene) in enumerate(zip(axes,['single_drop','many_drops'])):
        data=np.load(R/'fixtures'/f'{scene}_N4096.npz');field=data['field']
        ax.imshow(field,origin='lower',extent=[0,1,0,1],vmin=0,vmax=1,cmap=palette,interpolation='nearest',aspect='equal',rasterized=True)
        ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1]);ax.set_xlabel('x / L');ax.set_ylabel('y / L')
        ax.text(-.15,1.025,'ab'[i],transform=ax.transAxes,weight='bold',fontsize=9)
        for spine in ax.spines.values():spine.set_visible(True)
    fig.savefig(folder/'actual_vof_fields.png',dpi=300)
    fig.savefig(folder/'actual_vof_fields.pdf')
    plt.close(fig)
    caption='实际用于本轮单步算子测试的4096×4096体积分数场。a，单圆滴，6842个参与预测的界面格；b，256个不相交圆滴，127090个参与预测的界面格。蓝色表示液体，白色表示气体；两图保持真实物理比例并显示完整边界。这些是静态几何测试场，用于测量通量/曲率局部预测的精度与时间，未进行完整流体推进。'
    (folder/'actual_vof_fields-caption.txt').write_text(caption,encoding='utf-8')
    save(folder/'field_figure_manifest.json',dict(state='rendered',source='stored full4096x4096 VOF arrays,not reconstructed illustrations',physical_aspect='equal',inside_titles=False,full_boundaries=True,files=['actual_vof_fields.png','actual_vof_fields.pdf'],elapsed_seconds=time.perf_counter()-started))
    print('Actual VOF field figure rendered',time.perf_counter()-started)

if __name__=='__main__':main()
