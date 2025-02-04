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
repeat_num = 10
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
   'avg_target_percent',
   #  'avg_MAF', 'avg_inb', 'mean_rst',
   # 'mean_gv_JZRL', 'mean_gv_JZBBH', 'mean_gv_CZS',
   # 'var_gv_CZS', 'var_gv_JZRL', 'var_gv_JZBBH',
   # 'mean_pheno_CZS', 'mean_pheno_JZRL', 'mean_pheno_JZBBH',
   # 'var_pheno_CZS', 'var_pheno_JZRL', 'var_pheno_JZBBH',
   # 'heritability_CZS', 'heritability_JZRL', 'heritability_JZBBH',
    "target_gene_str",
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
    plt.subplot(1, 1, plt_idx)
    # liulian
    # if plt_idx == 1:
    #     plt.text(-1.8, 158, 'a)', fontsize=20)
    # if plt_idx == 4:
    #     plt.text(-1.8, 3.21, 'b)', fontsize=20)

    # sanlian
    # if plt_idx == 1:
    #     plt.text(-2.3, 0.1813, 'a)', fontsize=20)
    # if plt_idx == 2:
    #     plt.text(-2.3, 0.0245, 'b)', fontsize=20)
    # if plt_idx == 3:
    #     plt.text(-2.3, -1.02, 'c)', fontsize=20)
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
fig.savefig(f'{pic_path}/yilian.svg', bbox_inches='tight')
fig.savefig(f'{pic_path}/yilian.png', bbox_inches='tight')
plt.clf()

if "target_gene_str" in draw_titles:

    gene_dict = {
        1: 'One gene',
        2: 'Two genes',
        3: 'Three genes',
        4: 'Four genes',
        5: 'Five genes',
    }

    colors = ['white', '#376795', '#72BCD5', '#FFE6B7', '#F7AA58', '#E76254']
    fig = plt.gcf()
    ax = plt.axes()

    total_x_list = list(range(11 * len(scenarios)))
    total_x_ticks_tmp = ['S1', 'S2', 'S3', ''] * 11
    total_x_ticks = []
    for idx, x_tick in enumerate(total_x_ticks_tmp):
        if x_tick == 'S2':
            g_idx = idx // 4 + 6
            total_x_ticks.append(f'{x_tick}\nG{g_idx - 6}')
        else:
            total_x_ticks.append(x_tick)

    for s_idx in scenarios[:-1]:

        target_gene_distribution = {i: {j: 0.0 for j in range(6)} for i in range(6, 17)}
        with open(f'../results/S{s_idx}_target_distribution.csv', 'w', encoding='utf-8') as f:
            for g_idx in range(11):
                dis_list = []
                for gene_num in range(5):
                    tmp_y_list = []
                    for r_idx in range(1, repeat_num + 1):
                        # print(data[s_idx][r_idx]["target_gene_str"])
                        dis = data[s_idx][r_idx]["target_gene_str"][g_idx][gene_num]
                        tmp_y_list.append(int(dis))
                    tmp_y = sum(tmp_y_list) / len(tmp_y_list)
                    dis_list.append(tmp_y)
                    target_gene_distribution[g_idx + 6][gene_num + 1] = tmp_y
                f.writelines(f'G{g_idx},' + ','.join(list(map(str, dis_list))) + '\n')


        # print(target_gene_distribution)
        for g_idx in range(6, 17):
            #
            # target_gene_distribution[g_idx][0] += 5000
            # for i in range(1, 6):
            #     target_gene_distribution[g_idx][0] = (
            #             target_gene_distribution[g_idx][0] - target_gene_distribution[g_idx][i])

            for i in range(6):
                target_gene_distribution[g_idx][i] = (target_gene_distribution[g_idx][i] / 5000)
        # print(total_x_ticks)
        # print(target_gene_distribution)

        # if s_idx == 3:
        #     x_list = list(range(s_idx - 2, 10 * (len(scenarios)) + 1, len(scenarios)))
        #     # print(s_idx, x_list)
        # else:
        x_list = list(range(s_idx - 1, 11 * (len(scenarios)), len(scenarios)))
            # print(s_idx, x_list)

        # print(x_list)
        for gene_num in range(1, 6):
            values = []
            bottoms = []
            for generation in range(6, 17):
                tmp_value = 0
                tmp_bottom = 0
                for i in range(gene_num + 1):
                    if i == gene_num:
                        tmp_value += target_gene_distribution[generation][i]
                    else:
                        tmp_bottom += target_gene_distribution[generation][i]
                values.append(tmp_value)
                bottoms.append(tmp_bottom)
            if sum(values) == 0:
                continue
            # print(values)
            # print(bottoms)
            # print(gene_num)
            tmp_x_list = x_list[:-1]
            values = values[:-1]
            bottoms = bottoms[:-1]
            if s_idx == 3:
                plt.bar(tmp_x_list, values, bottom=bottoms, color=colors[gene_num], alpha=1,
                        label=f'{gene_dict[gene_num]}', width=0.8, linewidth=2, edgecolor='black')
            else:
                plt.bar(tmp_x_list, values, bottom=bottoms, color=colors[gene_num], alpha=1,
                        width=0.8, linewidth=2, edgecolor='black')

    plt.rcParams['font.sans-serif'] = ['Times New Roman']
    ax.tick_params(direction='in')
    # plt.title("target gene distribution")
    plt.legend(loc='best')
    # plt.xlabel('generation', fontsize=15)
    plt.xticks(total_x_list[:-6], labels=total_x_ticks[:-6])
    plt.ylim(0, 1)
    plt.show()
    pic_path = f'../images/target_gene_distribution_ruanzhu'
    fig.savefig(f'{pic_path}.svg', bbox_inches='tight')
    fig.savefig(f'{pic_path}.png', bbox_inches='tight')
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

