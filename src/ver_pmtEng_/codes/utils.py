import tiktoken
import re
import textwrap

class TextUtils:
    @staticmethod
    def count_tokens(text: str, model: str = 'gpt-4')->int:
        """
        Safely calculate the number of tokens in a text and handle non-string inputs
        Args:
            text (str): The input text to be tokenized.
            model (str): The model name to determine the encoding.Default is 'gpt-4'.
        Returns:
            int: The number of tokens in the input text.
        """

        if text is None:
            return 0
        
        if not isinstance(text, str):
            try:
                text = str(text)
            except Exception as e:
                print(f'Error converting input to string: {e}')
                return 0
            
        try:
            encoding = tiktoken.encoding_for_model(model)
        except KeyError:
            print(f'Warning: model {model} not found. Using cl100k_base encoding.')
            encoding = tiktoken.get_encoding("cl100k_base")

        try:
            return len(encoding.encode(text))
        except Exception as e:
            print(f'Error encoding text: {e}')
            return 0
        
    @staticmethod
    def cut_out_text(out_str: str, tol: int = 800) -> str:
        """
        Cut down the text if it exceeds the token limit by sampling lines.
        Args:
            out_str (str): The input text to be processed.
            tol (int): The token limit threshold. Default is 800.
        Returns:
            str: The processed text, potentially trimmed.
        """

        token_count = TextUtils.count_tokens(out_str, "gpt-4")  # 固定用gpt-4编码估算
        ratio = 1

        if token_count > tol:
            ratio = round(token_count / tol) + 1
            print(f"\n[Utils] 文本过长（{token_count} Token），超过阈值 {tol}，按比例 {ratio} 修剪")
        
        # 按比例抽取行
        lines = out_str.split('\n')
        if ratio != 1:
            lines = lines[::ratio]  # 每隔ratio行保留1行
            print(f"[Utils] 修剪后文本保留 {len(lines)} 行\n")
        
        return '\n'.join(lines)
    
    @staticmethod
    def remove_thinkmark(text: str) -> str:
        """提取文本中最后一个以 "```python" 开头的代码块（移除标记）"""
        start_marker = "```python"
        last_start_idx = text.rfind(start_marker)
        
        if last_start_idx == -1:
            return text
        
        # 查找代码块结束标记
        end_marker = "```"
        end_idx = text.find(end_marker, last_start_idx + len(start_marker))
        
        if end_idx == -1:
            # 无结束标记时，提取到文本末尾
            code = text[last_start_idx + len(start_marker):].strip()
        else:
            # 提取标记间的代码
            code = text[last_start_idx + len(start_marker):end_idx].strip()
        
        return code


    @staticmethod    
    def print_plot_to_streamlit(code_str: str, print_flag: bool = False, plot_flag: bool = False) -> str:
        """
        转换代码适配Streamlit：print→st.write，plt→st.pyplot
        :param code_str: 原始代码字符串
        :param print_flag: 是否替换print语句
        :param plot_flag: 是否替换绘图语句
        :return: 处理后的代码
        """

        lines = code_str.split('\n')
        processed_lines = []
        figure_var = "fig"
        have_savefig = False
        have_show = False

        for line in lines:
            stripped = line.lstrip()
            indent = len(line) - len(stripped)

            if stripped.startswith(('```python:', '```python')):
                continue

            elif stripped.startswith('```'):
                break

            elif stripped.startswith('print') and print_flag:
                modified = stripped.replace('print', 'st.write', 1)
                processed_lines.append(' ' * indent + modified)

            elif stripped.startswith('plt.figure') and plot_flag:
                modified = stripped.replace('plt.figure', f'{figure_var} = plt.figure', 1)
                processed_lines.append(' ' * indent + modified)

            elif stripped.startswith('plt.savefig') and plot_flag:
                processed_lines.append(' ' * indent + f'st.pyplot({figure_var})')
                have_savefig = True

            elif stripped.startswith('plt.show()') and plot_flag:
                modified = stripped.replace('plt.show()', f'st.pyplot({figure_var})', 1)
                have_show = True
                if not have_savefig:
                    processed_lines.append(' ' * indent + modified)
            
            elif stripped.startswith('```'):
                break
            
            else:
                processed_lines.append(line)

        return 'import streamlit as st\n' + '\n'.join(processed_lines)


    @staticmethod
    def remove_main_guard(code_str: str) -> str:
        """
        移除代码中的 if __name__ == "__main__": 块，并调整缩进
        """

        lines = code_str.split('\n')
        processed_lines= []
        in_main_block = False
        default_space = 4

        for i,line in enumerate(lines):
            stripped = line.strip()

            if stripped.startswith(('if __name__ == "__main__":', "if __name__ == '__main__':")):
                in_main_block = True

                if i + 1 < len(lines):
                    next_stripped = lines[i+1].strip()
                    default_space = len(lines[i+1]) - len(next_stripped)
                continue

            if in_main_block:
                current_indent = len(line) - len(line.strip())
                if current_indent >= default_space:
                    adjusted_indent = current_indent - default_space
                    processed_lines.append(' ' * adjusted_indent + line.strip())
                else:
                    in_main_block = False
                    processed_lines.append(line)
            else:
                processed_lines.append(line)

        return '\n'.join(processed_lines)