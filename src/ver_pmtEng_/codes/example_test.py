from re4_Agent import re4_SciComp_Agent

problem_description = r"""
give me a python code to calculate the eigenvalues and eigenvectors of a 3x3 matrix using numpy.
""".strip()


config1={
    "role": "consultant",
    "model": "gpt-4.1-mini-2025-04-14",
    "api_key": "your_apikey_here",
    "base_url": "https://api.chatanywhere.tech/v1"}

config2={
    "role": "programmer",
    "model": "gpt-4.1-mini-2025-04-14",
    "api_key": "your_apikey_here",
    "base_url": "https://api.chatanywhere.tech/v1"}

config3={
    "role": "reviewer",
    "model": "gpt-4.1-mini-2025-04-14",
    "api_key": "your_apikey_here",
    "base_url": "https://api.chatanywhere.tech/v1"}

agent1 = re4_SciComp_Agent(problem_description, 
                           consultant_config=config1, 
                           programmer_config=config2, 
                           reviewer_config=config3,
                           time_rev=1,)
agent1.run()