import os, sys
import MSA_Convex

args = sys.argv[1:]
root_path = args[0]

files_list = os.listdir(root_path)
files_list = [f for f in files_list if f.endswith(".msa")]
print(files_list)
for f in files_list:
    print("Running", f)
    full_path = os.path.join(root_path, f)
    try:
        MSA_Convex.solve(full_path)
    except:
        print("Fail")
