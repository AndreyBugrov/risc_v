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
    compiler_errors = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[1]
    if compiler_errors:
        error_message('\n'+compiler_errors.decode('utf-8'))


def create_cache_list_file(file_path: str):
    args = 'getconf -a | grep CACHE'
    cache_info = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True).communicate()[0].decode('utf-8').split('\n')[:-1]
    with open(file_path, 'w') as f:
        for line in cache_info:
            cache_level_info = line.split()
            if(len(cache_level_info)>1):
                f.write(f'{cache_level_info[1]}\n')
            else:
                f.write('0\n')


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


def run_matrix_exp(bin_path: str, function_name: str, matrix_sizes: list[int], exp_num: int, device_type: str, frequency: int, is_temporary: bool):
    min_n = matrix_sizes[0]
    max_n = matrix_sizes[1] + 1
    step = matrix_sizes[2]
    csv_file_name = function_name + '_' + device_type + '_' + str(frequency / (1000*1000)) + 'GHz'
    if is_temporary:
        csv_file_name = os.path.join('csv_results_tmp', csv_file_name)
    else:
        csv_file_name = os.path.join('csv_results', csv_file_name)
    with open(csv_file_name, 'w', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['Row size','OpenBLAS', 'Current', 'Ratio', 'Inaccuracy'])
        for i in range(min_n, max_n, step):
            num = str(i)
            args = f"{bin_path} a {function_name} {num} {exp_num}"
            cmd = shlex.split(args)
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            full_out = proc.communicate()
            out = full_out[0].decode('utf-8')
            if full_out[1]:
                print(f"Unhandled C++ exception:\n{full_out[1].decode('utf-8')}") #proc.communicate()[0].decode('utf-8')
            result = ';'.join(out[:-1].split('\n'))
            row = result.split(';')
            row.insert(0, num)
            writer.writerow(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Matrix multiplication experiments run automatization", formatter_class=argparse.RawTextHelpFormatter)
    simd_function_set = {'tr_omp_simd', 'row_simd', 'row_opt_simd','row_omp_simd', 'row_opt_omp_simd'}

    omp_function_set = {'base_omp', 'tr_omp', 'tr_omp_simd', 'row_omp', 'row_omp_simd', 'row_opt_omp', 'row_opt_omp_simd', 'strassen_omp', 'strassen_rec_omp'}
    omp_no_simd_function_set = omp_function_set.difference(simd_function_set)
    single_thread_function_set = {'base', 'row', 'row_simd', 'row_opt', 'row_opt_simd', 'tr', 'strassen'}
    single_thread_no_simd_function_set = single_thread_function_set.difference(simd_function_set)
    all_function_list = omp_function_set.union(single_thread_function_set)
    all_no_simd_function_set = all_function_list.difference(simd_function_set)
    all_function_list = list(all_function_list)
    all_function_list.sort()
    function_choices = all_function_list
    function_choices.extend(['omp', 'omp_no_simd' 'all', 'all_no_simd', 'single_thread', 'single_thread_no_simd'])
    parser.add_argument('-f', '--function-names', choices=function_choices, nargs='+', metavar='FUNC',
                        help="Matrix multiplication functions short names"
                        f"\nChoices:\n{function_choices}\n'all': all functions, 'omp': omp functions only, "
                        "single_thread': single-thread functions only. You can add '_no_simd' to 'all', 'omp' or "
                        "'single_thread' to avoid simd version of functions",
                        required=True)
    parser.add_argument('-s', "--matrix-sizes", help="Matrix sizes", metavar=('MIN_SIZE', 'MAX_SIZE', 'STEP'),type=int, nargs=3, required=True)
    parser.add_argument('-l', '--opt-level', help="Optimization level in the executable file",
                        choices=['release', 'opt', 'fast'], default='opt')
    parser.add_argument('-n', '--exp-num', help="Number of experiments with equal parameters", type=int, required=True)
    parser.add_argument('-d', '--device-name', help="RISC-V device name", choices=["sf2", "lichee", "mango", "kendryte", "x86"], required=True)
    parser.add_argument('--is-temporary', help="Should the results be saved to temporary directory or not", required=True, choices=['true', 'false'])
    args = parser.parse_args()
    function_names = args.function_name
    matrix_sizes = args.matrix_sizes
    opt_level = args.opt_level
    exp_num = args.exp_num
    device_name = args.device_name
    is_temporary = True if args.is_temporary == 'true' else False 
    if int(exp_num) < 1:
        error_message("choose at least one experiment!")
    if matrix_sizes[2] <= 0:
        error_message("step should be more than 0!")
    if matrix_sizes[0] > matrix_sizes[1]:
        error_message("min_n should be less or equal max_n!")

    type_handlings = {'release': '-O2', 'opt': '-O3', 'fast': '-Ofast'}
    
    root_source_dir = 'mat_opt'
    root_bin_dir = 'bin'

    function_name_list = {}
    if 'all' in function_names:
        function_name_list = all_function_list
    else:
        function_name_list = {item for item in function_names if item != 'omp' and item != 'omp_no_simd' and item != 'single_thread' and item != 'single_thread_no_simd'}
        if 'omp' in function_names:
            function_name_list = function_name_list.union(omp_function_set)
        if 'single_thread' in function_names:
            function_name_list = function_name_list.union(single_thread_function_set)
        if 'all_no_simd' in function_names:
            unction_name_list = function_name_list.union(all_no_simd_function_set)
        else:
            if 'omp_no_simd' in function_names:
                function_name_list = function_name_list.union(omp_no_simd_function_set)
            if 'single_thread_no_simd' in function_names:
                function_name_list = function_name_list.union(single_thread_no_simd_function_set)

    function_name_list = list(function_name_list)
    function_name_list.sort()

    bin_path = os.path.join(root_bin_dir, device_name + '_exp')
    source_file_list = [os.path.join(root_source_dir, 'common.cpp'),
                            os.path.join(root_source_dir, 'multiplication.cpp'),
                            os.path.join(root_source_dir, 'experiment.cpp'),
                            os.path.join(root_source_dir, 'main_experiment.cpp')]
    print("Preprocessing...")
    compile_source(source_file_list, bin_path, type_handlings[opt_level]) # recompilation does not garantees cold start (garanteed only for first function)

    create_cache_list_file('cache.txt')

    core_nums = get_available_cores()
    frequencies = get_min_max_frequencies()
    for core in core_nums:
        set_min_core_frequency_limit(frequencies[1], core)

    print("Conducting an experiment...")    
    for function_item in function_name_list:
        print(f'Process \"{function_item}\" function...')
        run_matrix_exp(bin_path, function_item, matrix_sizes, exp_num, device_name, frequencies[1], is_temporary)
    
    for core in core_nums:
        set_min_core_frequency_limit(frequencies[0], core)
    print("Done!")
