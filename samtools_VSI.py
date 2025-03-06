#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
from datetime import datetime

def run_command(cmd, description, cmd_log_file):
    print(f"Executing: {description}")
    # コマンド内容をCMD.txtに追記
    with open(cmd_log_file, "a") as f:
        f.write(f"{cmd}\n")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        sys.exit(f"Error: {description} failed.")
    
def main():
    parser = argparse.ArgumentParser(description="A script for SAM to BAM conversion, sorting, and indexing using samtools")
    parser.add_argument("-i", required=True, help="Input SAM file path")
    parser.add_argument("-o", required=True, help="Output directory path (must exist)")
    parser.add_argument("-@", type=int, default=1, help="Number of threads to use (default: 1)")
    
    args = parser.parse_args()
    
    input_sam = os.path.abspath(args.i)
    output_dir = os.path.abspath(args.o)
    threads = args.__dict__['@']  # -@ の値
    
    # 入力ファイルの存在確認
    if not os.path.isfile(input_sam):
        sys.exit(f"Error: Input file {input_sam} does not exist.")
    
    # 出力ディレクトリの存在確認
    if not os.path.isdir(output_dir):
        sys.exit(f"Error: Output directory {output_dir} does not exist.")
    
    # サブディレクトリ名の生成（YYMMDD_samtools_VSI_results）
    date_str = datetime.now().strftime("%y%m%d")
    result_dir = os.path.join(output_dir, f"{date_str}_samtools_VSI_results")
    
    if os.path.exists(result_dir):
        sys.exit(f"Error: Output subdirectory {result_dir} already exists. Overwriting is not allowed.")
    
    # サブディレクトリ作成
    os.makedirs(result_dir)
    print(f"Created result directory: {result_dir}")
    
    # CMD.txtのパスを設定
    cmd_log_file = os.path.join(result_dir, "CMD.txt")
    
    # 作業ディレクトリをサブディレクトリに変更
    os.chdir(result_dir)
    
    # 1. SAM -> BAM 変換
    bam_file = "alignment.bam"
    cmd_view = f"samtools view -bS {input_sam} -@ {threads} -o {bam_file}"
    run_command(cmd_view, "Conversion from SAM to BAM", cmd_log_file)
    
    # 2. BAMファイルのソート
    sorted_bam = "alignment.sorted.bam"
    cmd_sort = f"samtools sort {bam_file} -@ {threads} -o {sorted_bam}"
    run_command(cmd_sort, "Sorting BAM file", cmd_log_file)
    
    # 3. BAMファイルのインデックス作成
    cmd_index = f"samtools index {sorted_bam}"
    run_command(cmd_index, "Indexing BAM file", cmd_log_file)
    
    print("All processes completed successfully.")
    
if __name__ == '__main__':
    main()

