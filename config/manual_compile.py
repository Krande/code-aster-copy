import subprocess


def main():
    # Now call that batch file
    subprocess.run(link_bat_path.name, shell=True, cwd=ROOT_DIR)