import matplotlib.pyplot as plt
import numpy as np

target_gene_num = 5
scenarios = [1, 2, 3]
scenarios_name = {
    1: 'Cascading',
    2: 'Noncascading',
    3: 'Looping',
    4: 'Baseline'
}
repeat_num = 10
generation = 16
family_num = 5

data_path = f'G{target_gene_num}/mao_edit_5m5f_first_100_0816/gene_dis'

gene_dict = {
    1: 'One gene',
    2: 'Two genes',
    3: 'Three genes',
    4: 'Four genes',
    5: 'Five genes',
}
figure_list = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

total_x_list = list(range(11 * (family_num + 1)))
total_x_ticks_tmp = ['L1', 'L2', 'L3', 'L4', 'L5', ''] * 11
total_x_ticks = []
for idx, x_tick in enumerate(total_x_ticks_tmp):
    if x_tick == 'F3':
        g_idx = idx // 6 + 6
        total_x_ticks.append(f'{x_tick}\nG{g_idx - 6}')
    else:
        total_x_ticks.append(x_tick)
print(total_x_ticks)

fig = plt.figure()
# ax = plt.axes()
for s_idx in scenarios:

    colors = ['white', '#376795', '#72BCD5', '#FFE6B7', '#F7AA58', '#E76254']
    plt.subplot(1, 3, s_idx)

    target_gene_distribution = {
        i: {j: {k: 0.0 for k in range(6)} for j in range(1, family_num + 1)} for i in range(generation - 5)}
    # print(target_gene_distribution[0])

    for r_idx in range(1, repeat_num + 1):

        with open(data_path + f'/target_gene_dis_scenario_{s_idx}_repeat_{r_idx}.csv', encoding='utf-8') as f:

            lines = f.readlines()[1:]

            for g_idx, line in enumerate(lines):

                line_list = line.strip().split(',')[1:]

                for f_idx, family_dis in enumerate(line_list):

                    family_dis_list = family_dis.strip().split(':')[:-1]

                    for gene_num, sample_num in enumerate(family_dis_list):

                        target_gene_distribution[g_idx][f_idx + 1][gene_num + 1] += (float(sample_num) / 1000)

    for f_idx in range(1, family_num + 1):
        x_list = list(range(f_idx - 1, 11 * (family_num + 1), family_num + 1))

        for gene_num in range(1, 6):
            values = []
            bottoms = []
            for g_idx in range(generation - 7, -1, -1):
                tmp_value = 0
                tmp_bottom = 0
                for i in range(gene_num + 1):
                    if i == gene_num:
                        tmp_value += (target_gene_distribution[g_idx][f_idx][i] / 10)
                    else:
                        tmp_bottom += (target_gene_distribution[g_idx][f_idx][i] / 10)
                values.append(tmp_value)
                bottoms.append(tmp_bottom)
            if sum(values) == 0:
                continue
            tmp_x_list = x_list[:-2]
            values = values[:-1]
            bottoms = bottoms[:-1]
            print(len(tmp_x_list), len(values), len(bottoms))
            # print(gene_num)
            if f_idx == 3:
                plt.barh(tmp_x_list, values, left=bottoms, color=colors[gene_num], alpha=1,
                         label=f'{gene_dict[gene_num]}', height=0.8)
            else:
                plt.barh(tmp_x_list, values, left=bottoms, color=colors[gene_num], alpha=1,
                         height=0.8)

    plt.text(-0.15, 59, figure_list[s_idx - 1], fontsize=20)
    plt.rcParams['font.sans-serif'] = ['Times New Roman']
    # ax.tick_params(direction='in')
    if s_idx == 2:
        plt.legend(loc='best')
    plt.yticks(total_x_list[:-6], labels=total_x_ticks[:-6])
    for g_idx in range(generation - 7, -1, -1):
        plt.text(-0.115, g_idx * 6 + 1.52, f'G{9 - g_idx}')
    plt.xlim(0, 1)
    plt_title = f'target_gene_distribution_S{s_idx}'
    # plt.title(plt_title)
plt.show()
pic_path = f'images/gene_h'
fig.savefig(f'{pic_path}.svg', bbox_inches='tight')
fig.savefig(f'{pic_path}.png', bbox_inches='tight')
plt.clf()



