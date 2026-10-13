import os

# Config 参数
MAX_ITERATIONS = 2 # Reviewer 最大迭代次数
MULTI_MODAL_REVIEW = True # 是否启用多模态 Review 看图片
ALTERNATIVE_NUM = 2 # Consultant 生成的备选方案数量
MAX_COMMENTS = 5 # Reviewer 每轮最多评论数
ConsModule_ON = True # 是否启用 Consultant 模块
CUT_OUTEXT = 2000  # 截断输出文本的最大字符数 

# 设置工作区
WORKSPACE = os.environ.get("RE4_WORKSPACE", "./workspace")
os.makedirs(WORKSPACE, exist_ok=True)

# API 设置
Cons_name = os.environ.get("RE4_CONSULTANT_MODEL", os.environ.get("RE4_MODEL", "gpt-4.1-mini"))
Rev_name = os.environ.get("RE4_REVIEWER_MODEL", os.environ.get("RE4_MODEL", "gpt-4.1-mini"))

PROG_name = os.environ.get("RE4_PROGRAMMER_MODEL", os.environ.get("RE4_MODEL", "gpt-4.1-mini"))


# Supply OPENAI_API_KEY through your environment; never store credentials here.
OPENAI_API_BASE = os.environ.get("RE4_API_BASE", os.environ.get("OPENAI_BASE_URL", os.environ.get("OPENAI_API_BASE", "https://api.openai.com/v1")))

