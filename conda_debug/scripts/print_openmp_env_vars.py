import os


def main():
    print(f"{os.environ.get('KMP_VERSION')=}")
    print(f"{os.environ.get('KMP_AFFINITY')=}")
    print(f"{os.environ.get('KMP_SETTINGS')=}")
    print(f"{os.environ.get('OMP_DISPLAY_ENV')=}")
    print(f"{os.environ.get('MKL_VERBOSE')=}")
    print(f"{os.environ.get('MKL_DEBUG_CPU_TYPE')=}")
    print(f"{os.environ.get('OMP_NUM_THREADS')=}")
    print(f"{os.environ.get('MKL_NUM_THREADS')=}")
    print(f"{os.environ.get('MKL_DYNAMIC')=}")
    print(f"{os.environ.get('MKL_THREADING_LAYER')=}")


if __name__ == '__main__':
    main()
