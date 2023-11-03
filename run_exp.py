import csv
import subprocess
import argparse
import shlex


def compile_prog(source_path: str, prog_path: str, type: str):
    
    args = 'g++ ' + source_path + ' -o'
    if type == 'O3':
        args += type 
    args+= ' ' + prog_path
    cmd = shlex.split(args)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)


def run_pi_exp(prog_path: str, pi_args: list[str], output_fn: str):
    max_deg = int(pi_args[0])
    rectangle_type = 2 if len(pi_args) < 2 else pi_args[1]
    with open(output_fn, 'w') as f:
        writer = csv.writer(f, delimiter=';')
        for i in range(max_deg+1):
            num = str(pow(10, i))
            args = prog_path + ' ' + num + ' ' + rectangle_type
            cmd = shlex.split(args)
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            out = proc.communicate()[0]
            row = [num, out.decode('utf-8')]
            writer.writerow(row)


def run_matrix_exp(prog_path: str, matrix_args: list[str], output_fn: str):
    if len(matrix_args) < 2:
        min_n = int(matrix_args[0])
        max_n = int(matrix_args[1]) + 1
        step = int(matrix_args[2])
        with open(output_fn, 'w') as f:
            writer = csv.writer(f, delimiter=';')
            for i in range(min_n, max_n, step):
                num = str(pow(10, i))
                args = prog_path + ' ' + num
                cmd = shlex.split(args)
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                out = proc.communicate()[0]
                row = [num, out.decode('utf-8')]
                writer.writerow(row)
    ...


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
        pi_path = "./pi/pi"
        if pi_type == 'O3':
            pi_path += 'O3'
        elif pi_type == 'optimized':
            pi_path += 'o'
        compile_prog('./pi/pi.cpp', pi_path, pi_type)
        run_pi_exp(pi_path, pi_args, output_fn)
    if matrix_args:
        matrix_path = './matrix_mult/matrix_mult'
        if matrix_type == 'O3':
            matrix_path += 'O3'
        elif matrix_type == 'optimized':
            matrix_path += 'o'
        compile_prog('./matrix_mult/matrix_mult.cpp', matrix_path, matrix_type)
        run_matrix_exp(matrix_path,matrix_args,output_fn)
