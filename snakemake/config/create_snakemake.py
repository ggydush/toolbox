#!/usr/bin/env python
import shutil
import os


def check_if_file_in_dir(filename):
    """ Return boolean if file in directory """
    return os.path.exists(os.path.basename(filename))


def copy_file(f):
    """ Copy file to current directory unless already there """
    try:
        shutil.copy(f, os.path.basename(f))
    except Exception:
        print("Files are identical, nothing done.")


def copy_files_to_dir(files):
    """ Copy all relevant files to directory """
    for f in files:
        fname = os.path.basename(f)
        if check_if_file_in_dir(f):
            overwrite = False
            while overwrite not in ["y", "n"]:
                overwrite = input(
                    f"File: {fname} found in directory. Overwrite (y/n)? "
                ).lower()
                if overwrite in ["Y", "y"]:
                    copy_file(f)
        else:
            copy_file(f)


def main():
    """ Copy Snakemake files to current directory """
    base_dir = os.path.abspath(os.path.dirname(__file__))
    qsub = os.path.join(base_dir, "qsub_wrapper.py")
    snake = os.path.join(base_dir, "run_snakemake.sh")
    job = os.path.join(base_dir, "jobscript.sh")
    files = [qsub, snake, job]
    copy_files_to_dir(files)


if __name__ == "__main__":
    main()
