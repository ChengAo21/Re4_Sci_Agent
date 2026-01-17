import sys
import io
import warnings
import traceback
import textwrap
import re
from utils import TextUtils

class CodeProcessor:
    def __init__(self):
        self.exec_env = {}
        self.text_utils = TextUtils
    
    def capture_code_warnings(self, code: str) -> str:
        """
        执行代码并捕获警告、错误、输出
        :param code: 待执行的代码字符串
        :return: 包含执行结果、错误、警告的文本
        """

        code_clean = self.text_utils.remove_thinkmark(code)
        code_no_main = self.text_utils.remove_main_guard(code_clean)
        code_st = self.text_utils.print_plot_to_streamlit(code_no_main, print_flag=False, plot_flag=False)


        old_stdout, old_stderr = sys.stdout, sys.stderr
        new_stdout, new_stderr = io.StringIO(), io.StringIO()
        sys.stdout, sys.stderr = new_stdout, new_stderr


        captured_warnings = []
        seen_warnings = set()
        error_info = ""

        try:
            with warnings.catch_warnings(record=True) as w_list:
                warnings.simplefilter("always")
                exec(code_st, self.exec_env)


                for warn in w_list:
                    warn_key = (warn.filename, warn.lineno, str(warn.message))
                    if warn_key not in seen_warnings:
                        seen_warnings.add(warn_key)
                        captured_warnings.append(f"{warn.filename}:{warn.lineno} - {warn.message}")

  
            if not captured_warnings:
                error_info += textwrap.dedent("""
                    *The code runs without warnings, following is the output text content of the code. If there is no text content, it indicates that no text was printed as output.*
                    """)
            else:
                error_info += "\n**Warnings in code:**\n" + "\n".join(captured_warnings) + "\n"

        except Exception as e:
            # 捕获执行错误
            error_trace = traceback.format_exc().splitlines()
            error_info += "\n**Errors in code:**\n" + "\n".join(error_trace) + "\n"

        # 提取输出并修剪过长内容
        code_output = new_stdout.getvalue().strip()
        code_output_cut = self.text_utils.cut_out_text(code_output, tol=1200)
        code_error = new_stderr.getvalue().strip()

        # 组装最终结果
        error_info += "\n**Code output (stdout):**\n" + code_output_cut if code_output else ""
        error_info += "\n**Code error (stderr):**\n" + code_error if code_error else ""
        whole_output = f"\n**Executed Codes:**\n{code}\n\n{error_info}"

        # 恢复stdout/stderr
        sys.stdout, sys.stderr = old_stdout, old_stderr
        new_stdout.close()
        new_stderr.close()

        return whole_output
    

    def run_codes(self, pycode: str, output_flag: bool = False)->list:
        """
        :param pycode: 包含代码块的文本
        :param output_flag: 是否打印执行过程
        :return: 每个代码块的执行结果列表
        """

        code_pattern = re.compile(r'```python:?.*?\s*(.*)```', re.DOTALL)
        code_matches = code_pattern.findall(pycode)

        if not code_matches:
            return ["No valid Python code bolcks found"]

        valid_codes = []
        for code in code_matches:
            cleaned = re.sub(r'^```python:?\s*', '', code.strip(), flags=re.DOTALL)
            if cleaned:
                valid_codes.append(cleaned)

        execution_results = []
        for code in valid_codes:
            result = self.capture_code_warnings(code)
            execution_results.append(result)
            if output_flag:
                print(result)

        return execution_results        