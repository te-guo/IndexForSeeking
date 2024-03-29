###################
# $ python3 plot.py [output_name] [txt_1] [txt_2] [txt_3]
###################
import sys
import re
import matplotlib.pyplot as plt

n = 0
result_name = sys.argv[1]
colors = ['red', 'orange', 'green', 'blue', 'purple', 'pink', 'gray', 'black']
names = []
data = []

for filter in sys.argv[2:]:
    n += 1
    lines = []
    f = open("./log/" + filter + ".txt", 'r')
    names.append(filter)
    data.append([])
    while True:
        line = f.readline().split()
        if(line[0] == 'Evaluation:'):
            break
    while True:
        line = re.split("[\s\n=,][\s\n=,]*", f.readline())
        if len(line) <= 1:
            break
        else:
            lines.append(line)
    f.close()
    for line in lines:
        status = {}
        for i in range(1, len(line)-1, 2):
            status[line[i-1]] = float(line[i])
        data[-1].append(status)

fig = plt.figure()

subfig = fig.add_subplot(2, 3, 1)
subfig.set_title('Insert Throughput', fontsize=8)
subfig.set_xlabel('BPK', fontsize=5)
subfig.set_ylabel('Insert Throughput (M/s)', fontsize=5)
subfig.tick_params(axis='both', labelsize=5)
for i in range(n):
    x = [status['BPK'] for status in data[i]]
    y = [status['InsertTP'] for status in data[i]]
    subfig.plot(x, y, color = colors[i], linewidth = 0.8, linestyle='-', label=names[i])
    subfig.scatter(x, y, color = colors[i], s = 1, marker = 'o')
subfig.set_xlim(xmin=10)
subfig.set_ylim(ymin=0)


subfig = fig.add_subplot(2, 3, 2)
subfig.set_title('Query Throughput', fontsize=8)
subfig.set_xlabel('BPK', fontsize=5)
subfig.set_ylabel('Query Throughput (M/s)', fontsize=5)
subfig.tick_params(axis='both', labelsize=5)
for i in range(n):
    x = [status['BPK'] for status in data[i]]
    y = [status['QueryTP'] for status in data[i]]
    subfig.plot(x, y, color = colors[i], linewidth = 0.8, linestyle='-', label=names[i])
    subfig.scatter(x, y, color = colors[i], s = 1, marker = 'o')
subfig.set_xlim(xmin=10)
subfig.set_ylim(ymin=0)


subfig = fig.add_subplot(2, 3, 3)
subfig.set_title('I/O Cost', fontsize=8)
subfig.set_xlabel('BPK', fontsize=5)
subfig.set_ylabel('I/O cost', fontsize=5)
subfig.tick_params(axis='both', labelsize=5)
for i in range(n):
    x = [status['BPK'] for status in data[i]]
    y = [status['IO'] for status in data[i]]
    subfig.plot(x, y, color = colors[i], linewidth = 0.8, linestyle='-', label=names[i])
    subfig.scatter(x, y, color = colors[i], s = 1, marker = 'o')
subfig.set_xlim(xmin=10)
subfig.set_ylim(ymin=1-0.005)
# subfig.set_yscale('log')
subfig.legend(fontsize=4.5)


# subfig = fig.add_subplot(1, 2, 2)
# subfig.set_title('Standard Deviation of I/O Cost', fontsize=8)
# subfig.set_xlabel('BPK', fontsize=5)
# subfig.set_ylabel('Standard Deviation', fontsize=5)
# subfig.tick_params(axis='both', labelsize=5)
# for i in range(n):
#     x = [status['BPK'] for status in data[i]]
#     y = [status['SD_IO'] for status in data[i]]
#     subfig.plot(x, y, color = colors[i], linewidth = 0.8, linestyle='-', label=names[i])
#     subfig.scatter(x, y, color = colors[i], s = 1, marker = 'o')
# subfig.set_xlim(xmin=10)
# subfig.set_ylim(ymin=0)
# subfig.legend(fontsize=4.5)


subfig = fig.add_subplot(2, 3, 4)
subfig.set_title('Extra I/O Cost', fontsize=8)
subfig.set_xlabel('BPK', fontsize=5)
subfig.set_ylabel('I/O cost', fontsize=5)
subfig.tick_params(axis='both', labelsize=5)
for i in range(n):
    x = [status['BPK'] for status in data[i]]
    y = [status['IO'] - 1 for status in data[i]]
    subfig.plot(x, y, color = colors[i], linewidth = 0.8, linestyle='-', label=names[i])
    subfig.scatter(x, y, color = colors[i], s = 1, marker = 'o')
subfig.set_xlim(xmin=10)
subfig.set_yscale('log')


subfig = fig.add_subplot(2, 3, 5)
subfig.set_title('I/O Cost (by Hash Conflicts)', fontsize=8)
subfig.set_xlabel('BPK', fontsize=5)
subfig.set_ylabel('I/O cost', fontsize=5)
subfig.tick_params(axis='both', labelsize=5)
for i in range(n):
    x = [status['BPK'] for status in data[i]]
    y = [status['IO_hash'] for status in data[i]]
    subfig.plot(x, y, color = colors[i], linewidth = 0.8, linestyle='-', label=names[i])
    subfig.scatter(x, y, color = colors[i], s = 1, marker = 'o')
subfig.set_xlim(xmin=10)
subfig.set_ylim(ymin=-0.005)
# subfig.set_yscale('log')


subfig = fig.add_subplot(2, 3, 6)
subfig.set_title('I/O Cost (by Keys on the Left)', fontsize=8)
subfig.set_xlabel('BPK', fontsize=5)
subfig.set_ylabel('I/O cost', fontsize=5)
subfig.tick_params(axis='both', labelsize=5)
for i in range(n):
    x = [status['BPK'] for status in data[i]]
    y = [status['IO'] - 1 - status['IO_hash'] for status in data[i]]
    subfig.plot(x, y, color = colors[i], linewidth = 0.8, linestyle='-', label=names[i])
    subfig.scatter(x, y, color = colors[i], s = 1, marker = 'o')
subfig.set_xlim(xmin=10)
subfig.set_ylim(ymin=-0.0001)
# subfig.set_yscale('log')

#fig.suptitle('Result', fontsize=9)
fig.tight_layout(pad=0.7, w_pad=0.7, h_pad=0.7)
fig.savefig("./log/" + result_name + '.png', dpi=1000)