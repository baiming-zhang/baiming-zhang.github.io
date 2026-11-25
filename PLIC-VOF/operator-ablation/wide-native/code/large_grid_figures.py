"""Plot completed deployment measurements without selecting successful seeds."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from large_grid_operator_fixtures import R,save

COLORS={3:'#b6d7ec',5:'#6baed6',7:'#3182bd',9:'#08519c',11:'#08306b'}
CPU={'aphros':'#7b7164','simliquid':'#4b5563'}

def rows(name):
    with (R/'analysis'/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def number(row,key):
    value=row.get(key,'')
    return float(value) if value not in ('',None) else None

def old(row):return str(row['historical']).lower()=='true'
def adjusted(row):return row['policy'] not in ('control','baseline','historical')
def color(row):return '#8066a5' if old(row) else COLORS[int(row['stencil'])]
def marker(row):return '*' if old(row) else ('^' if adjusted(row) else 'o')

def main():
    status=json.loads((R/'analysis/status.json').read_text())
    if status['state']!='complete':raise RuntimeError('Completed formal measurements required; no placeholder figures')
    groups=rows('gpu_three_seed_groups.csv');records=rows('all_measured_records.csv')
    pinned=json.loads((R/'analysis/summary.json').read_text()).get('cpu_measurements')=='fixed_pool_shared_host_retest'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,
        'axes.linewidth':.65,'pdf.fonttype':42,'xtick.major.width':.6,'ytick.major.width':.6})
    folder=R/'figures';folder.mkdir(exist_ok=True)
    fig,axes=plt.subplots(2,2,figsize=(7.2,4.7))
    fig.subplots_adjust(left=.10,right=.76,bottom=.11,top=.95,wspace=.33,hspace=.40)
    for ri,op in enumerate(['flux','curvature']):
        for ci,mode in enumerate(['cpu_gpu_cpu','resident']):
            ax=axes[ri,ci]
            subset=[g for g in groups if g['scene']=='many_drops' and int(g['N'])==4096 and g['operator']==op and g['mode']==mode]
            for g in subset:
                x=number(g,'time_mean_ms');y=number(g,'full_task_mae')
                if y is None or y<=0 or number(g,'coverage_min')!=1:continue
                ax.errorbar(x,y,xerr=number(g,'time_seed_sd_ms'),yerr=number(g,'full_task_seed_sd'),
                    marker=marker(g),color=color(g),linestyle='none',markersize=5 if old(g) else 4,
                    capsize=1.5,elinewidth=.65,markeredgewidth=.7,
                    markerfacecolor='white' if adjusted(g) and not old(g) else color(g))
            for backend in CPU:
                refs=[r for r in records if r['scene']=='many_drops' and int(r['N'])==4096 and r['operator']==op and r['device']=='CPU' and r['implementation']==backend]
                refs.sort(key=lambda r:int(r['workers']))
                good=[r for r in refs if number(r,'full_task_mae') is not None and number(r,'full_task_coverage')==1]
                ax.plot([number(r,'time_median_ms') for r in good],[number(r,'full_task_mae') for r in good],
                    color=CPU[backend],marker='s' if backend=='aphros' else 'D',markersize=3,
                    linewidth=.65,markerfacecolor='none')
            ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Inference time (ms)')
            ax.set_ylabel('Flux-density MAE' if op=='flux' else r'Curvature MAE ($L^{-1}$)')
            ax.text(-.24,1.05,'abcd'[ri*2+ci],transform=ax.transAxes,weight='bold',fontsize=9)
            ax.grid(True,which='major',linewidth=.35,color='#e3e8ef');ax.set_axisbelow(True)
            for spine in ax.spines.values():spine.set_visible(True)
    handles=[Line2D([],[],color=c,marker='o',linestyle='none',label=f'{k}×{k} NN',markersize=4) for k,c in COLORS.items()]
    handles.extend([Line2D([],[],color='#8066a5',marker='*',linestyle='none',label='Earlier joint NN (3×3)',markersize=5),
        Line2D([],[],color='#3182bd',marker='o',linestyle='none',label='Baseline data',markersize=4),
        Line2D([],[],color='#3182bd',marker='^',markerfacecolor='white',linestyle='none',label='Adjusted data',markersize=4)])
    handles.extend(Line2D([],[],color=c,marker='s' if b=='aphros' else 'D',markerfacecolor='none',label=b.upper() if b=='aphros' else 'SimLiquid',linewidth=.65,markersize=3) for b,c in CPU.items())
    fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.78,.52),frameon=False,fontsize=6.5,handlelength=1.1,labelspacing=.8)
    fig.savefig(folder/'large_grid_accuracy_time.png',dpi=300);fig.savefig(folder/'large_grid_accuracy_time.pdf');plt.close(fig)
    caption='4096×4096、多液滴真实体积分数场上的单步精度与时间，参与计算的界面格为127090个。a、b为通量，c、d为曲率；左列包含CPU到GPU传输、推理与结果回传，右列为输入和输出均驻留GPU。新网络展示三个种子的均值与标准差，此前3×3联合网络为单个固定参考；传统算法分别显示1、8、24线程的实测点，连线只连接同一算法的线程配置。8、24线程使用外部持久线程池适配器。误差越小、时间越短越好；GPU采用单精度输入及原模型算术，传统数值算子采用双精度。仅计局部算子，不含邻域提取和流体推进。'
    if pinned:caption=caption.replace('1、8、24','1、4、8').replace('8、24线程使用外部持久线程池适配器。','4、8线程使用外部持久线程池适配器。CPU计算线程固定在0–7号逻辑核内，整机共享使用。')
    (folder/'large_grid_accuracy_time-caption.txt').write_text(caption,encoding='utf-8')
    fig,axes=plt.subplots(2,2,figsize=(7.2,4.7))
    fig.subplots_adjust(left=.10,right=.76,bottom=.11,top=.95,wspace=.33,hspace=.40)
    for ri,op in enumerate(['flux','curvature']):
        for ci,policy in enumerate(['control','enhanced']):
            ax=axes[ri,ci]
            for k in [3,5,7,9,11]:
                gs=[g for g in groups if g['scene']=='many_drops' and g['operator']==op and g['mode']=='cpu_gpu_cpu' and g['policy']==policy and int(g['stencil'])==k and not old(g)]
                gs.sort(key=lambda g:int(g['N']))
                assert len(gs)==4
                ax.errorbar([int(g['N']) for g in gs],[number(g,'time_mean_ms') for g in gs],
                    yerr=[number(g,'time_seed_sd_ms') for g in gs],color=COLORS[k],marker='o',markersize=3,linewidth=.8,capsize=1.5,elinewidth=.65)
            gs=[g for g in groups if g['scene']=='many_drops' and g['operator']==op and g['mode']=='cpu_gpu_cpu' and old(g)]
            gs.sort(key=lambda g:int(g['N']))
            ax.plot([int(g['N']) for g in gs],[number(g,'time_mean_ms') for g in gs],color='#8066a5',marker='*',markersize=4,linewidth=.8)
            for backend in CPU:
                fastest=[]
                for n in [512,1024,2048,4096]:
                    refs=[r for r in records if r['scene']=='many_drops' and int(r['N'])==n and r['operator']==op and r['device']=='CPU' and r['implementation']==backend and number(r,'full_task_coverage')==1]
                    fastest.append(min(refs,key=lambda r:number(r,'time_median_ms')) if refs else None)
                ax.plot([512,1024,2048,4096],[number(r,'time_median_ms') if r else np.nan for r in fastest],color=CPU[backend],marker='s' if backend=='aphros' else 'D',markerfacecolor='none',markersize=3,linewidth=.8)
            ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xticks([512,1024,2048,4096],labels=['512','1024','2048','4096'])
            ax.set_xlabel('Grid cells per side, N');ax.set_ylabel('Inference time (ms)')
            ax.text(-.24,1.05,'abcd'[ri*2+ci],transform=ax.transAxes,weight='bold',fontsize=9)
            ax.grid(True,which='major',linewidth=.35,color='#e3e8ef');ax.set_axisbelow(True)
            for spine in ax.spines.values():spine.set_visible(True)
    fig.legend(handles=handles[:6]+handles[8:],loc='center left',bbox_to_anchor=(.78,.52),frameon=False,fontsize=6.5,handlelength=1.1,labelspacing=.8)
    fig.savefig(folder/'large_grid_time_scaling.png',dpi=300);fig.savefig(folder/'large_grid_time_scaling.pdf');plt.close(fig)
    caption='多液滴场景随网格加密的单步推理时间，GPU时间包含传入、推理和传回。a、b为通量，c、d为曲率；左列使用基准训练数据，右列使用调整后的训练数据。每条新网络曲线固定窗口和数据策略，展示全部三个种子的均值与标准差；此前3×3联合网络单列。CPU曲线在各网格分别取1、8、24线程中实测最快的完整覆盖配置。网格由512²增至4096²时，实际计算界面格由15894增至127090，未计算纯相格来人为增加批量。'
    if pinned:caption=caption.replace('1、8、24','1、4、8')+'CPU计算线程使用固定8个逻辑核，宿主机可同时运行其他任务。'
    (folder/'large_grid_time_scaling-caption.txt').write_text(caption,encoding='utf-8')
    save(folder/'deployment_figure_manifest.json',dict(state='rendered_visual_check_pending',source='complete formal measured CSVs',files=[s+e for s in ['large_grid_accuracy_time','large_grid_time_scaling'] for e in ['.png','.pdf']],new_models='all three seeds',scene='many_drops',N=[512,1024,2048,4096]))
    print('Completed-data deployment figure rendered')

if __name__=='__main__':main()
