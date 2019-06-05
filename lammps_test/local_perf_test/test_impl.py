from subprocess import Popen, PIPE
import numpy as np
import sys


def make_test_set(command, nruns):
    average_run_times = None
    have_first = False
    it = 0
    while it < nruns:
        it += 1

        proc = Popen(command, shell=True, stdout=PIPE, stderr=None)
        proc.wait()
        ack = proc.communicate()

        times_str = ack[0].decode().strip()
        times_vec = []
        for single_launch_str in times_str.split(" | "):
            times_vec.append([int(v) for v in single_launch_str.split(" ")])

        if not have_first:
            average_run_times = np.array(times_vec)
            have_first = True
        else:
            average_run_times += np.array(times_vec)

    average_run_times = average_run_times / nruns
    return average_run_times

test_cpu_set = [1, 2, 4]
interaction_map = {0: 'None',
                   1: 'ii',
                   2: 'ei',
                   3: 'ee',
                   4: 'all'}

test_times_map = {'None' : {'raw': []},
                  'ii': {'raw': []},
                  'ei': {'raw': []},
                  'ee': {'raw': []},
                  'all': {'raw': []}}

for cpu_num in test_cpu_set:
    command = "./exec_test_set.sh 'mpirun -np {0} {1}'".format(cpu_num, sys.argv[1])
    exec_times = make_test_set(command, 1)

    for i, v in enumerate(exec_times):
        test_times_map[interaction_map[i]]['raw'].append(v)

for key in test_times_map.keys():
    test_times_map[key]['interaction'] = [test_times_map[key]['raw'][i] - v for i, v in enumerate(test_times_map['None']['raw'])]
    test_times_map[key]['raw_acc'] = test_times_map[key]['raw'][0][0] / np.array(test_times_map[key]['raw'])
    if key != 'None':
        test_times_map[key]['interaction_acc'] = test_times_map[key]['interaction'][0][0] / np.array(test_times_map[key]['interaction'])
    else:
        test_times_map[key]['interaction_acc'] = np.array(test_times_map[key]['interaction']) / (test_times_map[key]['interaction'][0][0]+0.1)

print("CPU_NUM,{0},\t,\t,\t,\t".format(",\t,\t,\t,".join(interaction_map[k] for k in [0, 1, 2, 3, 4])))
subcolumn=["raw","raw_acc","interaction","interaction_acc"]
print("\t,{0}".format(",".join(subcolumn*5)))

for i, cpu_num in enumerate(test_cpu_set):
    for entry in range(cpu_num):
        line = []
        for k in range(5):
            for subcol in subcolumn:
                line.append(test_times_map[interaction_map[k]][subcol][i][entry])

        print ("{0},{1}".format(cpu_num, ",".join(["{0:.3f}".format(v) for v in line])))
