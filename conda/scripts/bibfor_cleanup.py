import pathlib


def main():
    existing_def = pathlib.Path("../../bibfor/bibfor.def")
    missing_symbols = []
    with open('files/bibfor_def_error_1.txt', 'r') as f:
        for line in f:
            if not line.startswith('bibfor.def'):
                continue
            missing_symbol = line.split('unresolved external symbol')[-1].strip()
            missing_symbols.append(missing_symbol)

    print(f"Missing symbols: {missing_symbols}")
    new_file = existing_def.with_name("bibfor_new.def")
    new_lines = 0
    with new_file.open('w') as new_def:
        for line in existing_def.read_text().splitlines():
            symbol = line.strip()
            if symbol in missing_symbols:
                continue
            new_def.write(line+'\n')
            new_lines += 1

    print(f"New .def file created at: {new_file} with {new_lines-2} symbols")


if __name__ == "__main__":
    main()