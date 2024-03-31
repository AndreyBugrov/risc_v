import csv
import subprocess
import argparse
import shlex
import sys
import os


def error_message(msg: str):
    print(f"Error: {msg}")
    sys.exit(-1)


def compile_source(source_file_list: list[str], bin_path: str, optimization_flag: str):
    args = 'g++ ' + ' '.join(source_file_list) + ' -o ' + bin_path + ' ' + optimization_flag + ' -fopenmp -I open_blas/ -lopenblas'
    cmd = shlex.split(args)
    subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()


def get_available_cores():
    with open('/proc/cpuinfo', 'r') as f:
        return [int(line[:-1].split(': ')[1]) for line in f.readlines() if line.startswith('processor')]


def get_min_max_frequencies():
    """
    Returns:
        min and max frequencies
    """
    with open(f'/sys/devices/system/cpu/cpu0/cpufreq/scaling_available_frequencies', 'r') as f:
        avaliable_frequencies = [int(num) for num in f.readline()[:-1].split()]
    return min(avaliable_frequencies), max(avaliable_frequencies)


def set_min_core_frequency_limit(frequency, core_num):
    line = f"sudo sh -c 'echo {frequency} > /sys/devices/system/cpu/cpu{core_num}/cpufreq/scaling_min_freq'"
    subprocess.Popen(line, shell=True).communicate()


def run_matrix_exp(bin_path: str, function_name: str, matrix_sizes: list[int], exp_num: int, device_type: str, frequency: int):
    min_n = matrix_sizes[0]
    max_n = matrix_sizes[1] + 1
    step = matrix_sizes[2]
    csv_file_name = function_name + '_' + device_type + '_' + str(frequency / 1000*1000) + 'GHz'
    with open(os.path.join('csv_results', csv_file_name), 'w', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['Row size','OpenBLAS', 'Current', 'Ratio', 'Inaccuracy'])
        for i in range(min_n, max_n, step):
            num = str(i)
            args = f"{bin_path} a {function_name} {num} {exp_num}"
            cmd = shlex.split(args)
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            out = proc.communicate()[0].decode('utf-8')
            result = ';'.join(out[:-1].split('\n'))
            row = result.split(';')
            row.insert(0, num)
            writer.writerow(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Matrix multiplication experiments run automatization")
    parser.add_argument('-f', '--function-name', choices=['all', 'omp', 'single_thread', 'base', 'base_omp', 
                                                          'row', 'row_omp', 'tr', 'tr_omp', 'tr_omp_simd', 'strassen', 
                                                          'strassen_omp', 'strassen_rec_omp'], 
                        help="\tShort name for matrix multiplication function."
                        "\n'all': all functions, omp: omp functions only, single_thread: single-thread functions only",
                        required=True)
    parser.add_argument('-s', "--matrix-sizes", help="Matrix sizes: 1) min n 2) max n 3) step", type=int, nargs=3)
    parser.add_argument('-l', '--opt-level', help="Optimization level in the execution file",
                        choices=['release', 'opt', 'fast'], default='opt')
    parser.add_argument('-n', '--exp-num', help="Number of experiments with equal parameters", type=int, required=True)
    parser.add_argument('-d', '--device-name', help="RISC-V device name", choices=["sf2", "lichee", "mango", "kendryte", "x86"])
    args = parser.parse_args()
    function_name = args.function_name
    matrix_sizes = args.matrix_sizes
    opt_level = args.opt_level
    exp_num = args.exp_num
    device_name = args.device_name
    if int(exp_num) < 1:
        error_message("choose at least one experiment!")
    if matrix_sizes[2] <= 0:
        error_message("step should be more than 0!")
    if matrix_sizes[0] > matrix_sizes[1]:
        error_message("min_n should be less or equal max_n!")

    type_handlings = {'release': '-O2', 'opt': '-O3', 'fast': '-Ofast'}
    
    root_source_dir = 'mat_opt'
    root_bin_dir = 'bin'

    if function_name == 'all':
        function_name_list = ['base', 'base_omp', 'row', 'row_omp', 'tr', 'tr_omp', 'tr_omp_simd', 'strassen', 'strassen_omp', 'strassen_rec_omp']
    elif function_name == 'omp':
        function_name_list = ['base_omp', 'row_omp', 'tr_omp', 'tr_omp_simd', 'strassen_omp', 'strassen_rec_omp']
    elif function_name == 'single-thread':
        function_name_list = ['base', 'row', 'tr', 'strassen']
    else:
        function_name_list = [function_name]

    core_nums = get_available_cores()
    frequencies = get_min_max_frequencies()
    for core in core_nums:
        set_min_core_frequency_limit(frequencies[1], core)
    
    for function_item in function_name_list:
        bin_file_name = function_item + '_' + opt_level
        bin_path = os.path.join(root_bin_dir, bin_file_name)
        source_file_list = [os.path.join(root_source_dir, 'common.cpp'),
                            os.path.join(root_source_dir, 'multiplication.cpp'),
                            os.path.join(root_source_dir, 'experiment.cpp'),
                            os.path.join(root_source_dir, 'main_experiment.cpp')]
        compile_source(source_file_list, bin_path, type_handlings[opt_level])
        run_matrix_exp(bin_path, function_item, matrix_sizes, exp_num, device_name, frequencies[1])
        print(f"\"{function_item}\" function experiment is conducted")
    
    for core in core_nums:
        set_min_core_frequency_limit(frequencies[0], core)
