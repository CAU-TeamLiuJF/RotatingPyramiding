import matplotlib.pyplot as plt
import numpy as np

target_gene_num = 5
# scenarios = [3]
scenarios = [1, 2, 3, 4]
scenarios_file_name = [
    'Cascading',
    'NonCascading',
    'Rotating',
    'Baseline'
]
scenarios_name = {
    1: 'Cascading',
    2: 'Noncascading',
    3: 'Rotating(mRMC)',
    4: 'Baseline'
}
repeat_num = 100
generation = 15

data_path = f'../output/20250202_100'

data = {i: {j: {} for j in range(1, repeat_num + 1)} for i in scenarios}

for s_idx in scenarios:
    for r_idx in range(1, repeat_num + 1):

        with open(f'{data_path}/result_strategy_{scenarios_file_name[s_idx - 1]}_repeat_{r_idx}.csv', encoding='utf-8') as f:
            lines = f.readlines()

            # 属性
            header = lines[0].strip().split(',')
            line_data = lines[1:]
            if r_idx == 1:
                line_data = lines[6:]
            for key in header:
                data[s_idx][r_idx][key] = []
            for idx, line in enumerate(line_data):
                line_split = line.strip().split(',')[:-1]
                for i in range(len(header)):
                    if i != len(header) - 1:
                        # if i >= 8:
                        #     data[s_idx][r_idx][header[i]].append(float(line_split[i + 1]))
                        # else:
                        data[s_idx][r_idx][header[i]].append(float(line_split[i]))
                    # 最后五列是 不同阳性基因个数的个体分布
                    else:
                        data[s_idx][r_idx][header[i]].append(line_split[-5:])

            #

draw_titles = [
   # 'avg_target_percent',
    'avg_inb', 'avg_rel', 'avg_MAF', 'mean_rst'
   # 'mean_gv_JZRL', 'mean_gv_JZBBH', 'mean_gv_CZS',
   # 'var_gv_CZS', 'var_gv_JZRL', 'var_gv_JZBBH',
   # 'mean_pheno_CZS', 'mean_pheno_JZRL', 'mean_pheno_JZBBH',
   # 'var_pheno_CZS', 'var_pheno_JZRL', 'var_pheno_JZBBH',
   # 'heritability_CZS', 'heritability_JZRL', 'heritability_JZBBH',
   #  "target_gene_str",
]

color = ['#376795', '#FFD06F', '#E76254', '#72BCD5']

fig = plt.figure()
plt_idx = 0
for draw_title in draw_titles:
    print(draw_title)
    if draw_title == 'target_gene_str':
        continue
    plt_idx += 1
    plt.rcParams['font.sans-serif'] = ['Times New Roman']
    plt.rcParams['xtick.direction'] = 'in'  # 将x周的刻度线方向设置向内
    plt.rcParams['ytick.direction'] = 'in'  # 将y轴的刻度方向设置向内
    plt.subplot(2, 2, plt_idx)

    # sanlian
    if plt_idx == 1:
        plt.text(-1.8, 0.021, 'b)', fontsize=20)
    if plt_idx == 2:
        plt.text(-1.8, 0.028, 'd)', fontsize=20)
    if plt_idx == 3:
        plt.text(-1.8, 0.183, 'a)', fontsize=20)
    if plt_idx == 4:
        plt.text(-1.8, -1.38, 'c)', fontsize=20)
    x_list = [f'G{i}' for i in range(generation - 4)]
    with open(f'../results/{draw_title}.csv', 'w', encoding='utf-8') as f:
        for s_idx in scenarios:
            y_list = []
            for g_idx in range(11):
                tmp_y_list = []
                for r_idx in range(1, repeat_num + 1):
                    if draw_title.startswith('heritability'):
                        trait = draw_title.split('_')[-1]
                        var_gv = data[s_idx][r_idx][f'var_gv_{trait}'][g_idx]
                        var_pheno = data[s_idx][r_idx][f'var_pheno_{trait}'][g_idx]
                        tmp_y_list.append(var_gv / var_pheno)
                    else:
                        tmp_y_list.append(data[s_idx][r_idx][draw_title][g_idx])
                y_list.append(sum(tmp_y_list) / len(tmp_y_list))
            f.writelines(f'S{s_idx},' + ','.join(list(map(str, y_list))) + '\n')
            plt.plot(x_list[:-1], y_list[:-1], label=f'{scenarios_name[s_idx]}', color=color[s_idx - 1], marker='o', ms=7)


    plt.xticks(x_list[:-1])
    # plt.title(draw_title.replace('_', ' '))
    if plt_idx == 1:
        plt.legend(fontsize=12)
    # plt.xlabel(f'trait {(plt_idx - 1) % 3 + 1}', fontsize=15)
    # draw_title = draw_title.replace('rst', 'index')
    # plt.xlabel(draw_title.replace('_', ' '), fontsize=15)
    plt.grid(True)
plt.show()
pic_path = f'../images'
fig.savefig(f'{pic_path}/silian.svg', bbox_inches='tight')
fig.savefig(f'{pic_path}/silian.png', bbox_inches='tight')
plt.clf()

# target_gene_percent >= 100% 的世代数

for s_idx in scenarios[:-1]:
    generation_list = []
    for r_idx in range(1, repeat_num + 1):
        target_100_generation = 15
        for g_idx in range(10):
            target_gene_percent = data[s_idx][r_idx]['avg_target_percent'][g_idx]
            if target_gene_percent >= 1:
                target_100_generation = (g_idx + 6)
                break
        generation_list.append(target_100_generation)
    mean_generation = np.mean(generation_list)
    var_generation = np.var(generation_list)
    print(f'S{s_idx} mean {mean_generation} var {var_generation}')

