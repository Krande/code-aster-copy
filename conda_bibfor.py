import pathlib

THIS_DIR = pathlib.Path(__file__).parent


def main():
    bibfor_list = []
    root_dir = THIS_DIR / 'build/debug'

    for f in (root_dir / 'bibfor').rglob('*.o'):
        bibfor_list.append(f.relative_to(root_dir))

    with open(THIS_DIR / 'config' / 'bibfor.list', 'w') as f:
        f.write('\n'.join(map(str, bibfor_list)))

    print(bibfor_list)


if __name__ == '__main__':
    main()
