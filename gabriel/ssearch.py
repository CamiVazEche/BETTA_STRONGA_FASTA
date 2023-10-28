#! usr/bin/env pythong3

import os
import sys
import math
import re

#functions

def process_text_file(file_path):
    data = {}
    with open(file_path, 'r') as file:
        matrix_name = None
        for line in file:
            if line.startswith('#'):
                continue
            if matrix_name is None:
                matrix_name = line.strip()
                data[matrix_name] = {
                    'percid': [],
                    'alen': [],
                    'evalue': [],
                }
            else:
                parts = line.strip().split('\t')
                if len(parts) == 12:
                    percid, alen, evalue = float(parts[2]), int(parts[3]), float(parts[10])
                    data[matrix_name]['percid'].append(percid)
                    data[matrix_name]['alen'].append(alen)
                    data[matrix_name]['evalue'].append(evalue)
    return data

def process_directory(directory_path):
    results = {}
    for filename in os.listdir(directory_path):
        if filename.endswith('.txt'):
            file_path = os.path.join(directory_path, filename)
            unique_key = filename.rsplit('_', 1)[1][:-4]
            results.update({unique_key: process_text_file(file_path)})
    return results

#main

if __name__ == '__main__':
    directory_path = input("Enter the directory path containing text files: ")
    if os.path.isdir(directory_path):
        results = process_directory(directory_path)
        for unique_key, data in results.items():
            for matrix_name, metrics in data.items():
                print(f"Matrix: {unique_key}")
                percid_values = metrics['percid']
                alen_values = metrics['alen']
                evalue_values = metrics['evalue']
                print(f"    identity: {max(percid_values)}")
                print(f"    alignment length: {max(alen_values)}")
                print(f"    evalue: {min(evalue_values)}")
    else:
        print("Invalid directory path.")
