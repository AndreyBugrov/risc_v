import csv
import subprocess
import argparse
import shlex

import time


def compile_prog(source_path: str, prog_path: str, type: str):  
    args = 'g++ ' + source_path + ' -o'
    if type == 'O3': # can't compile with any '-O' option (even with O0)
        args += type 
    args+= ' ' + prog_path
    cmd = shlex.split(args)
    print(cmd)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    proc.communicate()


def run_pi_exp(prog_path: str, pi_args: list[str], output_fn: str):
    max_deg = int(pi_args[0])
    rectangle_type = '2' if len(pi_args) < 2 else pi_args[1]
    with open(output_fn, 'w') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['Число шагов','Погрешность','Время'])
        for i in range(max_deg+1):
            num = str(pow(10, i))
            args = prog_path + ' a ' + num + ' ' + rectangle_type
            cmd = shlex.split(args)
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            out = proc.communicate()[0].decode('utf-8')
            row = [num, out.split(' ')[0], out.split(' ')[1]]
            writer.writerow(row)


def run_matrix_exp(prog_path: str, matrix_args: list[str], output_fn: str):
    if len(matrix_args) < 4:
        min_n = int(matrix_args[0])
        step = int(matrix_args[1])
        max_n = int(matrix_args[2]) + 1
        with open(output_fn, 'w') as f:
            writer = csv.writer(f, delimiter=';')
            writer.writerow(['Число элементов в строке','Время'])
            for i in range(min_n, max_n, step):
                num = str(i)
                args = prog_path + ' a ' + num
                cmd = shlex.split(args)
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                out = proc.communicate()[0].decode('utf-8')
                row = [num, out]
                writer.writerow(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run experiments")
    parser.add_argument("--pi-args", help="Arguments for pi calculating: 1) max deg 2) rectangle type (unnecessary)",
                        nargs="*")
    parser.add_argument("--pi-type", help="Type of pi calculating execution file", choices=['O0', 'O3', 'optimized'])
    parser.add_argument("--matrix-args", help="Arguments for matrix multiplication: 1) min n 2)step 3) max n",
                        nargs="*")
    parser.add_argument("--matrix-type", help="Type of matrix multiplication execution file",
                        choices=['O0', 'O3', 'optimized'])
    parser.add_argument("-d", "--device", help="Device name", choices=["sf2", "laptop", "licheepi"], required=True)
    parser.add_argument("-o", "--output-fn", help="Path to csv file", required=True)
    args = parser.parse_args()
    pi_args = args.pi_args
    pi_type = args.pi_type
    matrix_args = args.matrix_args
    matrix_type = args.matrix_type
    device = args.device
    output_fn = args.output_fn
    if pi_args:
        pi_source = 'pi/pi.cpp'
        pi_path = "pi/pi"
        if pi_type == 'O3':
            pi_path += '_O3'
        elif pi_type == 'optimized':
            pi_source = 'pi/pi_optimized.cpp'
            pi_path += '_o'
        compile_prog(pi_source, pi_path, pi_type)
        run_pi_exp(pi_path, pi_args, output_fn)
    if matrix_args:
        matrix_source = 'matrix_mult/matrix_mult.cpp'
        matrix_path = 'matrix_mult/matrix_mult'
        if matrix_type == 'O3':
            matrix_path += '_O3'
        elif matrix_type == 'optimized':
            matrix_source = 'matrix_mult/matrix_mult_optimized.cpp'
            matrix_path += '_o'
        compile_prog(matrix_source, matrix_path, matrix_type)
        run_matrix_exp(matrix_path,matrix_args,output_fn)
