import sys
import os
import io
import re
import contextlib
import traceback
import matplotlib.pyplot as plt
import glob
from IPython.display import display, Image, Markdown, HTML
import base64
import time
from config import WORKSPACE, CUT_OUTEXT, MAX_ITERATIONS, MULTI_MODAL_REVIEW, ALTERNATIVE_NUM, MAX_COMMENTS, ConsModule_ON

# display funcs
def collapsible_mardown(title: str, content: str):
    """
    创建一个可折叠的 Markdown 区块，用于展示长文本内容。
    """
    html_content = f"""
    <details>
      <summary><strong>{title}</strong></summary>
      <pre style="white-space: pre-wrap;">{content}</pre>
    </details>
    """
    display(HTML(html_content))


# SAVE func
def save_response_toMD(resp: dict, oput_file: str = f"./{WORKSPACE}/temp_rept.md"):
    try:
        with open(oput_file, "w", encoding="utf-8") as f:
            f.write(f"### 🛌🏻 Re4gent with Multi-Modal Review\n\n")
            f.write(f"#### Configs 🏷️:\nMaxCount: {MAX_ITERATIONS}, MultiModalReview: {MULTI_MODAL_REVIEW}, AlternativeNum: {ALTERNATIVE_NUM}, MaxComments: {MAX_COMMENTS}, ConsModule_ON: {ConsModule_ON}, CUT_OUTEXT: {CUT_OUTEXT}\n\n")
            for key, value in resp.items():
                f.write(f"### {key}\n")
                if isinstance(value,list):
                    for i,item in enumerate(value):
                        if key == "python_codes":
                            f.write(f"\n#### Script block{i+1}:\n")
                            f.write(f"\n{item}\n\n")

                        elif key == "runtime_outputs":
                            f.write(f"\n#### Output block{i+1}\n")
                            f.write(f"\n{item}\n\n")
                        else:
                            f.write(f"\tCurrent Stage [{chr(i+96+1).upper()}/{len(value)}]\n{item}\n\n")
                    f.write("\n\n")
                else:
                    f.write(f"{str(value)}\n\n")
        print(f"✅ 保存成功: {oput_file}")
    except Exception as e:
        print(f"❌ {oput_file}保存失败: {e}")


# text truncate func
def truncate_outext(text: str, max_chars: int = 800) -> str:
    """
    高效截断过长文本：保留开头和结尾，中间替换
    """
    if not text:
        return ""
        
    if len(text) <= max_chars:
        return text
        
    head_len = max_chars // 2
    tail_len = max_chars // 2
    
    # 字符串切片，内存开销最小
    return (
        text[:head_len] 
        + f"\n\n... [Truncated: Content too long ({len(text)} chars). Hidden middle part.] ...\n\n" 
        + text[-tail_len:]
    )

def execute_code_tool(code_input: str) -> str:
    """
    执行 Python 代码并捕获输出（stdout/stderr）和错误。
    会自动切换到 workspace 目录以保存图片。
    """
    # 1. 代码去除 ```python 标记
    pattern = r"```python.*?\s+(.*?)(?:```|$)"
    # re.DOTALL 让 . 可以匹配换行符，re.IGNORECASE 忽略 python 大小写
    match = re.search(pattern, code_input, re.DOTALL | re.IGNORECASE)
    if match:
        code = match.group(1).strip()
    else:
        # 如果没匹配到 python 标记，尝试匹配裸的 ```
        pattern_bare = r"```.*?\s+(.*?)(?:```|$)"
        match_bare = re.search(pattern_bare, code_input, re.DOTALL)
        if match_bare:
            code = match_bare.group(1).strip()
        else:
            # 兜底：如果完全没有 markdown 标记，假设整段文本就是代码
            code = code_input.strip()
    # 如果代码是空的，直接返回
    if not code:
        return "Error: No executable code found."


    # 2. 执行环境
    output_capture = io.StringIO()
    # 限制 global 作用域，但保留必要的 builtins
    local_scope = {}
    global_scope = {"__name__": "__main__", "__builtins__": __builtins__}

    # 保存当前目录，执行完后切回来
    original_cwd = os.getcwd()
    
    try:
        # 切换到工作目录, plt.savefig('fig.png') 就会存在这里
        os.chdir(WORKSPACE)
        
        # 清理旧图片
        for f in glob.glob("*.png"): os.remove(f)

        # 3. 执行代码并重定向输出
        with contextlib.redirect_stdout(output_capture), contextlib.redirect_stderr(output_capture):
            exec(code, global_scope, global_scope)
            
        raw_result = output_capture.getvalue()
        result = truncate_outext(raw_result, max_chars=CUT_OUTEXT)

    except Exception as e:
        # 捕获 Traceback，这对 Reviewer 修复 bug 非常重要
        tb_str = traceback.format_exc()
        raw_result = output_capture.getvalue() + f"\nRuntime Error:\n{tb_str}"
        result = truncate_outext(raw_result, max_chars=CUT_OUTEXT)

    finally:
        # 恢复目录
        os.chdir(original_cwd)
        # 关闭所有 plot，防止内存泄漏
        plt.close('all')

    return result

def load_image_tool():
    """
    从 workspace 读取最新生成的图片，以可折叠 Markdown 展示。
    """
    image_files = glob.glob(os.path.join(WORKSPACE, "*.png")) + \
                  glob.glob(os.path.join(WORKSPACE, "*.jpg"))
    
    image_files.sort(key=os.path.getmtime)
    
    if not image_files:
        return None

    # 构建折叠内容
    img_html = ""
    for img_path in image_files:
        filename = os.path.basename(img_path)

        img_html += f'<p><b>{filename}</b></p>\n'
        img_html += f'<img src="{img_path}" style="max-width:50%; height:auto; margin:10px 0;">\n'

    # 用 details 折叠
    collapsible_html = f"""
    <details>
        <summary><b>Generated Visualization ({len(image_files)} image(s)) 📊</b></summary>
        <div style="padding: 10px; margin-top: 5px; border-left: 2px solid #ddd;">
            {img_html}
        </div>
    </details>
    """
    
    display(HTML(collapsible_html))
    return image_files

def encode_image_tool(image_path: str):
    """
    图片编码为 base64 字符串，作为输入用于Review
    """
    with open(image_path, "rb") as img_file:
        encoded_string = base64.b64encode(img_file.read()).decode('utf-8')
    return encoded_string