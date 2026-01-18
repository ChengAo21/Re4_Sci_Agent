import os
from openai import OpenAI
from config import *

llm_client = OpenAI(base_url=os.environ['OPENAI_API_BASE'], api_key=os.environ['OPENAI_API_KEY'])

# LLM call 
def llm_call_(config_dict, structured_module):
    # .parse()
    response_ = llm_client.beta.chat.completions.parse(
        model=config_dict["llm_name"],
        messages=[
            {"role": "system", "content": config_dict["system_content"]},
            {"role": "user", "content": config_dict["human_content"]}
        ],
        response_format=structured_module,
        temperature=0.0,
        )

    return response_.choices[0].message.parsed