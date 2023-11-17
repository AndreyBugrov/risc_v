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
    cmd = shlex.split(args)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    proc.communicate()



def run_pi_exp(exe_path: str, pi_args: list[str], output_fn: str):
    max_deg = int(pi_args[0])
    rectangle_type = '2' if len(pi_args) < 2 else pi_args[1]
    with open(output_fn, 'w') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['Число шагов','Погрешность','Время'])
        for i in range(max_deg+1):
            num = str(pow(10, i))
            args = exe_path + ' a ' + num + ' ' + rectangle_type
            cmd = shlex.split(args)
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            out = proc.communicate()[0].decode('utf-8')
            row = [num, out.split(' ')[0], out.split(' ')[1]]
            writer.writerow(row)


def run_matrix_exp(exe_path: str, matrix_args: list[str], output_fn: str):
    if len(matrix_args) < 4:
        min_n = int(matrix_args[0])
        step = int(matrix_args[1])
        max_n = int(matrix_args[2]) + 1
        with open(output_fn, 'w') as f:
            writer = csv.writer(f, delimiter=';')
            writer.writerow(['Число элементов в строке','Время'])
            for i in range(min_n, max_n, step):
                num = str(i)
                args = exe_path + ' a ' + num
                cmd = shlex.split(args)
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                out = proc.communicate()[0].decode('utf-8')
                row = [num, out]
                writer.writerow(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run experiments")
    parser.add_argument('-p', "--pi-args", help="Arguments for pi calculating: 1) max deg 2) rectangle type (unnecessary)",
                        nargs="*")
    parser.add_argument('-m', "--matrix-args", help="Arguments for matrix multiplication: 1) min n 2)step 3) max n",
                        nargs="*")
    parser.add_argument('-t', '--exe-type', help="Type of execution file", choices=['normal', 'fast', 'optimized'], required=True)
    parser.add_argument("-o", "--output-fn", help="Path to csv file", required=True)
    args = parser.parse_args()
    pi_args = args.pi_args
    matrix_args = args.matrix_args
    exe_type = args.exe_type
    output_fn = args.output_fn

    if pi_args and matrix_args:
        error_message("choose only one mode")
    if not pi_args and not matrix_args:
        error_message("choose any mode")
    
    type_handlings = {'normal': ('normal', '-O3'), 'fast': ('fast', '-Ofast'), 'optimized': ('optimized', '-O3')}   
    source_path = ''
    exe_path = ''

    if pi_args:
        exe_path = 'pi/pi'
    elif matrix_args:
        exe_path = 'matrix_mult/matrix_mult'
    
    if exe_type == 'optimized':
        source_path = exe_path + '_optimized.cpp'
    else:
        source_path = exe_path + '.cpp'

    exe_path += '_' + type_handlings[exe_type][0]

    compile_prog(source_path, exe_path, type_handlings[exe_type][1])
    if pi_args:
        run_pi_exp(exe_path, pi_args, output_fn)
    if matrix_args:
        run_matrix_exp(exe_path, matrix_args, output_fn)
