import os

class File_Handler:
    def __init__(self,path_sv: str,print_flag: bool = False):
        """
        param path_save: 结果保存根目录
        """
        self.path_save = self._check_and_create_dir(path_sv)
        self.print_flag = print_flag
        print(f"[FileHandler] save path: {self.path_save}. Print_save_flag: {self.print_flag}\n")

    def _check_and_create_dir(self, path: str) -> str:
        
        abs_path = os.path.abspath(path)
        if not os.path.exists(abs_path):
            os.makedirs(abs_path, exist_ok=True)
            print(f"[FileHandler] make absPath: {abs_path}\n")
        return abs_path

    def write_content(self, content: str, file_prefix: str, count: int = 0) -> None:

        file_path = os.path.join(self.path_save, f"{file_prefix}_{count}.dat")
        with open(file_path, "w", ) as f:
            f.write(f"\n{'*' * 40}\n")
            f.write(content)
            f.write(f"\n{'*' * 40}\n")
        if self.print_flag:
            print(f"[FileHandler] done save {file_prefix}_{count}.dat\n")

    def write_content_code_out(self, code_out_list: list, file_prefix: str, count: int = 0) -> None:

        file_path = os.path.join(self.path_save, f"{file_prefix}_{count}.dat")
        with open(file_path, "w", ) as f:
            f.write(f"\n{'*' * 40}\n")
            for item in code_out_list:
                f.write(item)
            f.write(f"\n{'*' * 40}\n")
        if self.print_flag:
            print(f"[FileHandler] done save: {file_prefix}_{count}.dat\n")

    def read_content(self, file_prefix: str, count: int = 0) -> str:

        file_path = os.path.join(self.path_save, f"{file_prefix}_{count}.dat")
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"[FileHandler] File not found: {file_path}")
        with open(file_path, "r") as f:
            content = f.read()
        return content