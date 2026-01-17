from openai import OpenAI
from typing import Optional


class LLMAgent:
    def __init__(self, model_name: str, api_key: str, base_url: Optional[str] = None,):
        """
        :param model_name: 模型名称（gpt-4.1-mini-2025-04-14、deepseek-reasoner）
        :param api_key: 模型API密钥
        :param base_url: 模型API基础地址（非OpenAI官方需指定）
        """
        self.model_name = model_name if model_name else "gpt-4.1-mini-2025-04-14"

        self.api_key = api_key 
        self.base_url = base_url

        self.system_prompt = f"""
        [ROLE: CODE EXPERT] You are a code expert specializing in scientific computing and numerical methods.
        Provide deep technical insights into algorithm implementation, code optimization, and numerical solution strategies.
        Use programming-specific terminology to analyze code structures, identify bugs, and offer advanced coding solutions.
        """.strip()

        self.client = self.create_client()


    def create_client(self):
        """根据model_name + api_key + url, 创建OpenAI客户端

           Default url is chatanywhere---  'https://api.chatanywhere.tech/v1'

           api_deepseek = "sk-"
           api_doubao = 'sk-'
           api_chatanywhere = "sk-"

        """
        if not self.api_key:
            raise ValueError("API key is required to create OpenAI client.")
        
        # For deepseek
        elif self.model_name in ['deepseek-chat','deepseek-reasoner']:
            client = OpenAI(api_key=self.api_key, base_url=self.base_url if self.base_url is not None else "https://api.deepseek.com/v1")

        # For doubao
        elif self.model_name in ["doubao-1.5-thinking-vision-pro",'doubao-1.5-vision-pro','doubao-1.5-thinking-pro']:
            client = OpenAI(api_key = self.api_key, base_url =self.base_url if self.base_url is not None else  'https://wcode.net/api/gpt/v1')

        # For gemini, claude, chatgpt, etc.
        else:
            client=OpenAI(api_key=self.api_key, base_url="https://api.chatanywhere.tech/v1")

        return client


    def execute(self, user_prompt: str, temperature: float = 0.0, ) -> str:
        """
        调用LLM生成响应
        :param user_prompt: 用户提示词（问题/需求）
        :return: LLM生成的文本
        """
        # 构建消息列表
        messages = [
            {"role": "system", "content": self.system_prompt},
            {"role": "user", "content": user_prompt}
        ]

        try:
            response = self.client.chat.completions.create(
                model=self.model_name,
                messages=messages,
                temperature=temperature,
                stream=False
            )
            return response.choices[0].message.content.strip()
        except Exception as e:
            raise RuntimeError(f"failed to use LLM api: {str(e)}") from e