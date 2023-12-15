import csv
import subprocess
import argparse
import shlex
import sys


def error_message(msg: str):
    print(f"Error: {msg}")
    sys.exit(-1)


def compile_prog(source_path: str, exe_path: str, optimization_flag: str):
    args = 'g++ ' + source_path + ' -o ' + exe_path + ' ' + optimization_flag
    if exe_path.endswith('omp'):
        args+= ' -fopenmp'
    cmd = shlex.split(args)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    proc.communicate()



def run_pi_exp(exe_path: str, pi_args: list[str], exp_num: str, output_fn: str):
    max_deg = int(pi_args[1]) if len(pi_args) >=2 else 9
    rectangle_type = '2'
    with open(output_fn, 'w', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['Число шагов','Погрешность','Время'])
        for i in range(2, max_deg+1):
            num = str(pow(10, i))
            args = exe_path + ' a ' + num + ' ' + rectangle_type + ' ' + exp_num
            cmd = shlex.split(args)
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            out = proc.communicate()[0].decode('utf-8')
            row = [num, out.split(' ')[0], out.split(' ')[1]]
            writer.writerow(row)


def run_matrix_exp(exe_path: str, matrix_args: list[str], exp_num: str, output_fn: str):
    if len(matrix_args) >= 4:
        min_n = int(matrix_args[1])
        max_n = int(matrix_args[2]) + 1
        step = int(matrix_args[3])
        with open(output_fn, 'w', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter=';')
            writer.writerow(['Число элементов в строке','Время'])
            for i in range(min_n, max_n, step):
                num = str(i)
                args = exe_path + ' a ' + num + ' ' + exp_num
                cmd = shlex.split(args)
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                out = proc.communicate()[0].decode('utf-8')
                row = [num, out]
                writer.writerow(row)
    else:
        error_message("Too few matrix arguments")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run experiments")
    parser.add_argument('-p', "--pi-args", help="Arguments for pi calculating: "
                        "1) double (d) or long double (ld) version 2) max deg (better for 9-th) - not necessary", 
                        nargs="*")
    parser.add_argument('-m', "--matrix-args", help="Arguments for matrix multiplication: "
                        "1) double (d) or long double (ld) version 2) min n 3) max n 4) step",
                        nargs="*")
    parser.add_argument('-t', '--exe-type', help="Type of execution file", choices=['release', 'normal', 'fast', 'omp', 'optimized'], required=True)
    parser.add_argument('-n', '--exp-num', help="Number of experiments with equal parameters", required=True)
    parser.add_argument("-o", "--output-fn", help="Path to csv file", required=True)
    args = parser.parse_args()
    pi_args = args.pi_args
    matrix_args = args.matrix_args
    exe_type = args.exe_type
    exp_num = args.exp_num
    output_fn = args.output_fn

    if pi_args and matrix_args:
        error_message("choose only one mode")
    if not pi_args and not matrix_args:
        error_message("choose any mode")
    if int(exp_num) < 1:
        error_message("choose at least one experiment")
    
    type_handlings = {'release': ('release', '-O2'), 'normal': ('normal', '-O3'), 'fast': ('fast', '-Ofast'), 
                      'omp': ('omp','-O2'), 'optimized': ('optimized', '-O3')}   
    source_path = ''
    exe_path = ''

    if pi_args:
        exe_path = 'pi/pi'
        if pi_args[0] != 'd' and pi_args[0] != 'ld':
            error_message("wrong type of pi calculating")
        exe_path += '_' + pi_args[0]
    elif matrix_args:
        exe_path = 'matrix_mult/matrix_mult'
        if matrix_args[0] != 'd' and matrix_args[0] != 'ld':
            error_message("wrong type of matrix multiplication")
        exe_path += '_' + matrix_args[0]
    
    if exe_type == 'optimized':
        source_path = exe_path + '_optimized.cpp'
    elif exe_type == 'omp':
        source_path = exe_path + '_omp.cpp'
    else:
        source_path = exe_path + '.cpp'

    exe_path += '_' + type_handlings[exe_type][0]

    compile_prog(source_path, exe_path, type_handlings[exe_type][1])
    if pi_args:
        run_pi_exp(exe_path, pi_args, exp_num, output_fn)
    if matrix_args:
        run_matrix_exp(exe_path, matrix_args, exp_num, output_fn)
