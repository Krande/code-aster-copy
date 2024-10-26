import pathlib


def main():
    existing_def = pathlib.Path("../../bibfor/bibfor.def")
    new_file = existing_def.with_name("bibfor_new.def")
    total_lines = 0
    new_lines = 0
    with new_file.open('w') as new_def:
        for line in existing_def.read_text().splitlines():
            total_lines += 1
            if 'module' in line or '.' in line:
                continue
            new_def.write(line+'\n')
            new_lines += 1

    print(f"New .def file created at: {new_file} with {new_lines} symbols out of {total_lines} lines")


if __name__ == "__main__":
    main()