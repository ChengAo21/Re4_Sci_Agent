import os
from openai import OpenAI
from config import OPENAI_API_BASE


def get_llm_client():
    """Read credentials at call time, without importing or displaying a stored key."""
    api_key = os.environ.get("OPENAI_API_KEY", "").strip()
    if not api_key:
        raise RuntimeError("Set OPENAI_API_KEY in your environment before running the LLM workflow.")
    return OpenAI(base_url=OPENAI_API_BASE, api_key=api_key)

# LLM call 
def llm_call_(config_dict, structured_module):
    # .parse()
    response_ = get_llm_client().beta.chat.completions.parse(
        model=config_dict["llm_name"],
        messages=[
            {"role": "system", "content": config_dict["system_content"]},
            {"role": "user", "content": config_dict["human_content"]}
        ],
        response_format=structured_module,
        temperature=0.0,
        )

    return response_
