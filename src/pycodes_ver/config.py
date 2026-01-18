import os

# Config 参数
MAX_ITERATIONS = 2 # Reviewer 最大迭代次数
MULTI_MODAL_REVIEW = True # 是否启用多模态 Review 看图片
ALTERNATIVE_NUM = 2 # Consultant 生成的备选方案数量
MAX_COMMENTS = 5 # Reviewer 每轮最多评论数
ConsModule_ON = True # 是否启用 Consultant 模块
CUT_OUTEXT = 2000  # 截断输出文本的最大字符数 

# 设置工作区
WORKSPACE = "./workspace"
os.makedirs(WORKSPACE, exist_ok=True)

# API 设置
Cons_name = "gpt-4.1-mini"  # Consultant
Rev_name = "gpt-4.1-mini"   # Reviewer

PROG_name = "gpt-4.1-mini"  # or deepseek-reasoner / gemini


os.environ['OPENAI_API_KEY'] = 'sk-chatanywhere-yourkeyhere'
os.environ['OPENAI_API_BASE'] = 'https://api.chatanywhere.tech/v1'

# os.environ['OPENAI_API_KEY'] = 'sk-deepseek-yourkeyhere'
# os.environ['OPENAI_API_BASE'] = 'https://api.deepseek.com'
