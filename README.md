# Re4-Scientific-Computing-Agent-with-Rewriting-Resolution-Review-and-Revision
The proposed agent, incorporating a "rewriting-resolution-review-revision" logical chain via three reasoning LLMs (functioning as the Consultant, Reviewer, and Programmer, respectively), is integrated in a collaborative and interactive manner.
![schematic](schematic_re4.png)

For more details, see our arXiv preprint [arXiv:2508.20729](https://img.shields.io/badge/arXiv-2508.20729-b31b1b.svg)(https://arxiv.org/abs/2508.20729)

## Future work
```
- release streamlit version and employ langchain
```

## Project Structure
```
src/
├── utils.py             # General tool functions
├── FileHandler.py       # File operation
├── code_executor.py     # Code processing & execution  
├── llm_agent.py         # Client configuration
├── prompts_template.py  # Prompt generation
└── re4_Agent.py         # Workflow
```

## Contact
If you have any questions on our work or implementation, feel free to reach out to [chengao23@mails.ucas.ac.cn](mailto:chengao23@mails.ucas.ac.cn)!

If you find this repository useful, please consider giving a star ⭐ and cite our paper.

```
@article{cheng2025re4,
  title={Re4: Scientific Computing Agent with Rewriting, Resolution, Review and Revision},
  author={Cheng, Ao and Zhang, Lei and He, Guowei},
  journal={arXiv preprint arXiv:2508.20729},
  year={2025}
}
```
