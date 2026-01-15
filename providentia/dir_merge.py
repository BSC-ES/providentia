import os 
import shutil
from tqdm import tqdm

def merge_and_remove(exp, mod):
    if os.path.exists(exp): 
        if os.listdir(exp):
            for root, dirs, files in os.walk(exp):
                # compute relative path from exp
                rel_path = os.path.relpath(root, exp)
                dest_dir = os.path.join(mod, rel_path)

                # create destination directory if needed and keep permissions
                if not os.path.exists(dest_dir):
                    src_mode = os.stat(root).st_mode
                    os.makedirs(dest_dir)
                    os.chmod(dest_dir, src_mode)

                # copy files (overwrite if exists)
                if files:
                    for file in tqdm(files, bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"Copying files to {dest_dir}"):
                        src_file = os.path.join(root, file)
                        dst_file = os.path.join(dest_dir, file)
                        shutil.copy2(src_file, dst_file)

        # remove old directory
        print(f"\nRemoving source directory: {exp}\n")
        shutil.rmtree(exp)

mod_root = os.path.expanduser("~/data/providentia/mod")
exp_root = os.path.expanduser("~/data/providentia/exp")
mod_to_interp_root = os.path.expanduser("~/data/providentia/mod_to_interp")
exp_to_interp_root = os.path.expanduser("~/data/providentia/exp_to_interp")

os.makedirs(mod_root, exist_ok=True)
os.makedirs(mod_to_interp_root, exist_ok=True)

merge_and_remove(exp_root, mod_root)
merge_and_remove(exp_to_interp_root, mod_to_interp_root)