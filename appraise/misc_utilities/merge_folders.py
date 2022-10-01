"""
A small tool to merge unzipped large folders downloaded from Google Drive
Author: Xiaozhe Ding
"""
import os
import sys
folder_path = sys.argv[1]
if folder_path[-1] != '/':
    folder_path += '/'
for i in range(19):
    merge_command = 'ditto {} {}'.format(folder_path[:-1] + '\ ' + str(i+2) + '/', folder_path)
    print(merge_command)
    os.system(merge_command)
