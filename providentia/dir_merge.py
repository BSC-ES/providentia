import os
import shutil
from tqdm import tqdm
import yaml

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

# get current path and providentia root path
PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

data_paths = yaml.safe_load(
    open(os.path.join(PROVIDENTIA_ROOT, "settings", "data_paths.yaml"))
)


def merge_and_remove(exp, mod):
    if os.path.exists(exp):
        # rename the directory if there is no mod directory
        if not os.path.exists(mod):
            os.rename(exp, mod)

            print(
                f"\nRename of '{os.path.basename(exp)}' to '{os.path.basename(mod)}' finished!"
            )

        # copy from exp to mod if both exist
        else:
            if os.listdir(exp):
                for root, _, files in os.walk(exp):
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
                        for file in tqdm(
                            files,
                            total=len(files),
                            desc=f"Copying to {dest_dir}",
                            bar_format="{l_bar}{bar}|{n_fmt}/{total_fmt}",
                        ):
                            src_file = os.path.join(root, file)
                            dst_file = os.path.join(dest_dir, file)
                            shutil.copy2(src_file, dst_file)

                print("\nCopy finished!")

            # remove old directory
            print(f"\nRemoving source directory: {exp}\n")
            shutil.rmtree(exp)
    else:
        print(f"\nNo '{os.path.basename(exp)}' directory already, all good! Exiting...")


def dir_handler(mod_dir, exp_dir):
    if f"{mod_dir}_root" in data_paths["local"]:
        mod_root = os.path.expanduser(data_paths["local"][f"{mod_dir}_root"])

        if mod_root.endswith(mod_dir):
            exp_root = os.path.join(os.path.dirname(mod_root), exp_dir)
            merge_and_remove(exp_root, mod_root)
        else:
            print(
                "\nThe directory name is not {mod_dir} due to the user changing it in 'data_paths.yaml'. Exiting..."
            )
    else:
        print(
            f"\nOld 'data_paths.yaml'. Execute this again after running 'git pull'. Exiting..."
        )


dir_handler("mod", "exp")
dir_handler("mod_to_interp", "exp_to_interp")
