import csv
import subprocess
import argparse
import shlex
import sys
import os
import logging


logger = logging.getLogger(__name__)


def create_function_dict() -> dict[str, set[str]]:
    all_function_set = {'base', 'base_omp','tr', 'tr_omp', 'tr_omp_simd', 'row', 'row_omp', 'row_opt', 'row_opt_omp',
                        'row_simd', 'row_omp_simd', 'row_opt_omp_simd', 'row_opt_simd', 'strassen', 'strassen_omp',
                        'strassen_rec_omp'}
    omp_function_set = {item for item in all_function_set if 'omp' in item}
    single_thread_function_set = all_function_set.difference(omp_function_set)
    simd_function_set = {item for item in all_function_set if 'simd' in item}
    omp_function_set_no_simd = omp_function_set.difference(simd_function_set)
    single_thread_function_set_no_simd = single_thread_function_set.difference(simd_function_set)
    function_names_dict = {'all': all_function_set, 'omp': omp_function_set,
                        'single_thread': single_thread_function_set, 'omp_no_simd': omp_function_set_no_simd,
                        'single_thread_no_simd': single_thread_function_set_no_simd}
    for item in all_function_set:
        function_names_dict[item] = {item}
    return function_names_dict


def critical_message(msg: str):
    logger.critical(f"{msg}")
    sys.exit(-1)


def compile_source(source_file_list: list[str], bin_path: str, optimization_flag: str):
    args = 'g++ ' + ' '.join(source_file_list) + ' -o ' + bin_path + ' ' + optimization_flag + ' -fopenmp -I open_blas/ -lopenblas'
    cmd = shlex.split(args)
    compiler_errors = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[1]
    if compiler_errors:
        critical_message('Compilation errors:\n'+compiler_errors.decode('utf-8'))


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
    csv_file_name = function_name + '_' + device_type + '_' + str(frequency / (1000*1000)) + 'GHz' +'.csv'
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
        out = proc.communicate()[0].decode('utf-8')  # we can't catch termination message (it is neither in stderr nor in stdout)
        if proc.returncode:
            critical_message(f'Process has completed with non-zero return code: {proc.returncode}')
        result = ';'.join(out[:-1].split('\n'))
        row = result.split(';')
        row.insert(0, num)
        with open(csv_file_name, 'a', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter=';')
            writer.writerow(row)
        logger.debug(f'Experiment with {i}x{i} matrix was carried out')
    logger.info(f'Results saved to {csv_file_name}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Matrix multiplication experiments run automatization", formatter_class=argparse.RawTextHelpFormatter)
    function_names_dict = create_function_dict()
    function_names_list = list(function_names_dict.keys())
    function_names_list.sort()
    parser.add_argument('-f', '--function-names', choices=function_names_dict.keys(), nargs='+', metavar='FUNC',
                        help="Matrix multiplication functions short names"
                        f"\nChoices:\n{function_names_list}\n'all': all functions, 'omp': omp functions only, "
                        "single_thread': single-thread functions only. You can add '_no_simd' to 'all', 'omp' or "
                        "'single_thread' to avoid simd version of functions",
                        required=True)
    parser.add_argument('-s', "--matrix-sizes", help="Matrix sizes", metavar=('MIN_SIZE', 'MAX_SIZE', 'STEP'),type=int, nargs=3, required=True)
    parser.add_argument('-l', '--opt-level', help="Optimization level in the executable file",
                        choices=['release', 'opt', 'fast'], default='opt')
    parser.add_argument('-n', '--exp-num', help="Number of experiments with equal parameters", type=int, required=True)
    parser.add_argument('-d', '--device-name', help="RISC-V device name", choices=["sf2", "lichee", "mango", "kendryte", "x86"], required=True)
    parser.add_argument('--is-temporary', help="Should the results be saved to temporary directory or not", required=True, choices=['true', 'false'])
    parser.add_argument('-r', '--recompile',help='Should source files be recompiled', choices=['true', 'false'], default='true')
    args = parser.parse_args()
    function_names = args.function_names
    matrix_sizes = args.matrix_sizes
    opt_level = args.opt_level
    exp_num = args.exp_num
    device_name = args.device_name
    is_temporary = True if args.is_temporary == 'true' else False
    recompile = True if args.recompile == 'true' else False

    logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s [%(name)s:%(lineno)s] %(message)s')
    
    if is_temporary:
        logger.warning('Results will be saved to temporary directory')
    else:
        logger.warning('Results will be saved to the base directory')
    if int(exp_num) < 1:
        critical_message("Choose at least one experiment!")
    if matrix_sizes[2] <= 0:
        critical_message("Step should be more than 0!")
    if matrix_sizes[0] > matrix_sizes[1]:
        critical_message("\"min_n\" should be less or equal \"max_n\"!")

    type_handlings = {'release': '-O2', 'opt': '-O3', 'fast': '-Ofast'}
    
    root_source_dir = 'mat_opt'
    root_bin_dir = 'bin'

    function_names_set = set()
    for function_item in function_names:
        function_names_set = function_names_set.union(function_names_dict[function_item])
    logger.debug(f'Chosen functions: {function_names_set}')

    bin_path = os.path.join(root_bin_dir, device_name + '_exp')
    source_file_list = [os.path.join(root_source_dir, 'common.cpp'),
                            os.path.join(root_source_dir, 'multiplication.cpp'),
                            os.path.join(root_source_dir, 'experiment.cpp'),
                            os.path.join(root_source_dir, 'main_experiment.cpp')]
    logger.info("Start of preprocessing phase")
    if recompile:
        logger.info('Compilation')
        compile_source(source_file_list, bin_path, type_handlings[opt_level]) # recompilation does not garantees cold start (garanteed only for first function)
    else:
        logger.warning('Compilation is skipped')

    create_cache_list_file('cache.txt')

    core_nums = get_available_cores()
    frequencies = get_min_max_frequencies()
    try:
        logger.info('Frequency setting')
        for core in core_nums:
            set_min_core_frequency_limit(frequencies[1], core)

        logger.info("Start of the experiment execution phase")
        for function_item in function_names_set:
            logger.info(f'Process \"{function_item}\" function')
            run_matrix_exp(bin_path, function_item, matrix_sizes, exp_num, device_name, frequencies[1], is_temporary)
    except KeyboardInterrupt:
        critical_message('Program has been interrupted')
    finally:
        for core in core_nums:
            set_min_core_frequency_limit(frequencies[0], core)
        print("Done!")
